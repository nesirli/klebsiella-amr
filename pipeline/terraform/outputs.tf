output "job_queue_name" {
  value = aws_batch_job_queue.this.name
}

output "job_queue_arn" {
  value = aws_batch_job_queue.this.arn
}

output "work_bucket" {
  value = aws_s3_bucket.work.bucket
}

output "aws_region" {
  value = var.aws_region
}
