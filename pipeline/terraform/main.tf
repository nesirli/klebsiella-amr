terraform {
  required_version = ">= 1.5"
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
    cloudinit = {
      source  = "hashicorp/cloudinit"
      version = "~> 2.3"
    }
    random = {
      source  = "hashicorp/random"
      version = "~> 3.6"
    }
  }
}

provider "aws" {
  region = var.aws_region
}

# ---------------------------------------------------------------------------
# Network: reuse the default VPC's public subnets (no NAT gateway -> no
# ongoing hourly cost when idle; instances get internet access via the
# default VPC's internet gateway since default-VPC subnets auto-assign
# public IPs).
# ---------------------------------------------------------------------------

data "aws_vpc" "default" {
  default = true
}

data "aws_subnets" "default" {
  filter {
    name   = "vpc-id"
    values = [data.aws_vpc.default.id]
  }
}

resource "aws_security_group" "batch" {
  name        = "${var.name_prefix}-batch-sg"
  description = "Outbound-only SG for Batch compute environment instances"
  vpc_id      = data.aws_vpc.default.id

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

# ---------------------------------------------------------------------------
# IAM: Batch service role + EC2 instance role/profile
# ---------------------------------------------------------------------------

data "aws_iam_policy_document" "batch_assume" {
  statement {
    actions = ["sts:AssumeRole"]
    principals {
      type        = "Service"
      identifiers = ["batch.amazonaws.com"]
    }
  }
}

resource "aws_iam_role" "batch_service" {
  name               = "${var.name_prefix}-batch-service-role"
  assume_role_policy = data.aws_iam_policy_document.batch_assume.json
}

resource "aws_iam_role_policy_attachment" "batch_service" {
  role       = aws_iam_role.batch_service.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
}

data "aws_iam_policy_document" "ec2_assume" {
  statement {
    actions = ["sts:AssumeRole"]
    principals {
      type        = "Service"
      identifiers = ["ec2.amazonaws.com"]
    }
  }
}

resource "aws_iam_role" "ecs_instance" {
  name               = "${var.name_prefix}-ecs-instance-role"
  assume_role_policy = data.aws_iam_policy_document.ec2_assume.json
}

resource "aws_iam_role_policy_attachment" "ecs_instance" {
  role       = aws_iam_role.ecs_instance.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

# S3 access for Nextflow task file staging (scoped to the work bucket below).
resource "aws_iam_role_policy" "ecs_instance_s3" {
  name = "${var.name_prefix}-s3-workdir-access"
  role = aws_iam_role.ecs_instance.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect   = "Allow"
      Action   = ["s3:GetObject", "s3:PutObject", "s3:DeleteObject", "s3:ListBucket"]
      Resource = [aws_s3_bucket.work.arn, "${aws_s3_bucket.work.arn}/*"]
    }]
  })
}

resource "aws_iam_instance_profile" "ecs_instance" {
  name = "${var.name_prefix}-ecs-instance-profile"
  role = aws_iam_role.ecs_instance.name
}

# ---------------------------------------------------------------------------
# S3 bucket for Nextflow's remote work directory
# ---------------------------------------------------------------------------

resource "aws_s3_bucket" "work" {
  bucket        = "${var.name_prefix}-nf-work-${data.aws_caller_identity.current.account_id}"
  force_destroy = true
}

data "aws_caller_identity" "current" {}

# ---------------------------------------------------------------------------
# Launch template: standard ECS-optimized AMI + AWS CLI v2 installed via
# user data (Nextflow's aws.batch executor needs the CLI on the host to
# stage files to/from S3; only needed since we're not using Wave/Fusion).
# ---------------------------------------------------------------------------

data "aws_ssm_parameter" "ecs_ami" {
  name = "/aws/service/ecs/optimized-ami/amazon-linux-2/recommended/image_id"
}

data "cloudinit_config" "batch" {
  gzip          = false
  base64_encode = true

  part {
    content_type = "text/x-shellscript"
    content      = <<-EOF
      #!/bin/bash
      exec >>/var/log/aws-cli-install.log 2>&1
      set -x
      yum install -y unzip
      curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip
      unzip -q -o /tmp/awscliv2.zip -d /tmp
      /tmp/aws/install -i /opt/aws-cli -b /opt/aws-cli/bin

      # The bundled AWS CLI binary depends on libz.so.1 from its own
      # dist/ dir, but that's not on the default linker search path --
      # harmless on this host (which has a system libz), but breaks when
      # Nextflow mounts /opt into minimal biocontainers images that don't.
      # Wrap the binary so cliPath resolves to something self-contained.
      mv /opt/aws-cli/bin/aws /opt/aws-cli/bin/aws-real
      cat > /opt/aws-cli/bin/aws <<'WRAPPER'
      #!/bin/bash
      export LD_LIBRARY_PATH="/opt/aws-cli/v2/current/dist:$LD_LIBRARY_PATH"
      exec /opt/aws-cli/bin/aws-real "$@"
      WRAPPER
      chmod +x /opt/aws-cli/bin/aws
    EOF
  }
}

resource "aws_launch_template" "batch" {
  name_prefix = "${var.name_prefix}-batch-lt-"
  image_id    = data.aws_ssm_parameter.ecs_ami.value
  user_data   = data.cloudinit_config.batch.rendered

  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size = 30
      volume_type = "gp3"
    }
  }
}

# ---------------------------------------------------------------------------
# Batch compute environment + job queue
# ---------------------------------------------------------------------------

resource "random_id" "ce_suffix" {
  byte_length = 4
}

resource "aws_batch_compute_environment" "this" {
  compute_environment_name = "${var.name_prefix}-ce-${random_id.ce_suffix.hex}"
  type                     = "MANAGED"
  service_role             = aws_iam_role.batch_service.arn

  compute_resources {
    type                = "EC2"
    allocation_strategy = "BEST_FIT_PROGRESSIVE"
    min_vcpus           = 0
    max_vcpus           = var.max_vcpus
    # This account is restricted to Free Tier instance types only.
    # m7i-flex.large (2 vCPU / 8GB) is enough for smoke-testing the
    # Batch/S3 wiring but too small for the pipeline's heavier stages
    # (kraken2 needs 8GB just for its DB, assembly needs 16GB) -- lift
    # this restriction before running the real pipeline on Batch.
    instance_type       = ["m7i-flex.large"]
    subnets             = data.aws_subnets.default.ids
    security_group_ids  = [aws_security_group.batch.id]
    instance_role       = aws_iam_instance_profile.ecs_instance.arn

    launch_template {
      launch_template_id = aws_launch_template.batch.id
      version             = "$Latest"
    }
  }

  lifecycle {
    create_before_destroy = true
  }

  depends_on = [aws_iam_role_policy_attachment.batch_service]
}

resource "aws_batch_job_queue" "this" {
  name     = "${var.name_prefix}-queue"
  state    = "ENABLED"
  priority = 1

  compute_environment_order {
    order               = 1
    compute_environment = aws_batch_compute_environment.this.arn
  }
}
