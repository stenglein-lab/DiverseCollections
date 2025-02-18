#!/usr/bin/env bash

nextflow run nf-core/fetchngs \
   -profile <docker/singularity/.../institute> \
   --input ids.csv \
   --outdir <OUTDIR>