#!/usr/bin/env bash

for file in *.fastq.gz; do
  if [[ ! "$file" =~ _[12]\.fastq\.gz$ ]]; then
    mv "$file" "${file%.fastq.gz}_1.fastq.gz"
  fi
done