variable "aws_region" {
  type    = string
  default = "us-east-1"
}

variable "name_prefix" {
  type    = string
  default = "kleb-amr-nf"
}

variable "max_vcpus" {
  type    = number
  default = 8
}
