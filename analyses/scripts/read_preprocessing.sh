#!/usr/bin/env bash

nextflow run main.nf \
   -resume \
   -profile singularity \
   --fastq_dir ../fetchngs/results/fastq/ \
   --fastq_pattern *_[12].fastq.*