#!/usr/bin/env bash

# Use this script with the consensus sequences for galbut/chaq positive samples from geneious to capture all possible sequences


# Build bowtie indexes
bowtie2-build 1020-D3_RNA1_ref_supp.fasta 1020-D3_RNA1_ref_supp
bowtie2-build 1020-F1_RNA1_ref_supp.fasta 1020-F1_RNA1_ref_supp
bowtie2-build 1020-F3_Chaq_ref_supp.fasta 1020-F3_Chaq_ref_supp
bowtie2-build 1020-F3_RNA3_new_ref.fasta 1020-F3_RNA3_new_ref
bowtie2-build 1020-M-19_RNA3_new_ref.fasta 1020-M-19_RNA3_new_ref
bowtie2-build 1020-M-0818-A7_RNA1_ref_supp.fasta 1020-M-0818-A7_RNA1_ref_supp
bowtie2-build 1428-M-3_RNA3_new_ref.fasta 1428-M-3_RNA3_new_ref
bowtie2-build 1428-M-42_RNA3_new_ref.fasta 1428-M-42_RNA3_new_ref
bowtie2-build 1428-M-43_Chaq_new_ref.fasta 1428-M-43_Chaq_new_ref

# Run samples against new index
bowtie2 -x 1020-D3_RNA1_ref_supp -U ../results/host_filtered_fastq/1020-D3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-D3_RNA1_ref_supp_self_realign1.bt.sam
bowtie2 -x 1020-F1_RNA1_ref_supp -U ../results/host_filtered_fastq/1020-F1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F1_RNA1_ref_supp_self_realign.bt.sam
bowtie2 -x 1020-F3_Chaq_ref_supp -U ../results/host_filtered_fastq/1020-F3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_Chaq_ref_supp_self_realign.bt.sam
bowtie2 -x 1020-F3_RNA3_new_ref -U ../results/host_filtered_fastq/1020-F3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_RNA3_new_ref_self_realign.bt.sam
bowtie2 -x 1020-M-19_RNA3_new_ref -U ../results/host_filtered_fastq/1020-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-19_RNA3_new_ref_self_realign.bt.sam
bowtie2 -x 1020-M-0818-A7_RNA1_ref_supp -U ../results/host_filtered_fastq/1020-M-0818-A7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A7_RNA1_ref_supp_self_realign.bt.sam
bowtie2 -x 1428-M-3_RNA3_new_ref -U ../results/host_filtered_fastq/1428-M-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-3_RNA3_new_ref_self_realign.bt.sam
bowtie2 -x 1428-M-42_RNA3_new_ref -U ../results/host_filtered_fastq/1428-M-42_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-42_RNA3_new_ref_self_realign.bt.sam
bowtie2 -x 1428-M-43_Chaq_new_ref -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_Chaq_new_ref_self_realign.bt.sam

bowtie2 -x 1020-D3_RNA1_ref_supp -U ../results/host_filtered_fastq/1020-D3_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-D3_RNA1_ref_supp_self2_realign1.bt.sam
bowtie2 -x 1020-F1_RNA1_ref_supp -U ../results/host_filtered_fastq/1020-F1_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F1_RNA1_ref_supp_self2_realign.bt.sam
bowtie2 -x 1020-F3_Chaq_ref_supp -U ../results/host_filtered_fastq/1020-F3_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_Chaq_ref_supp_self2_realign.bt.sam
bowtie2 -x 1020-F3_RNA3_new_ref -U ../results/host_filtered_fastq/1020-F3_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_RNA3_new_ref_self2_realign.bt.sam
bowtie2 -x 1020-M-19_RNA3_new_ref -U ../results/host_filtered_fastq/1020-M-19_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-19_RNA3_new_ref_self2_realign.bt.sam
bowtie2 -x 1020-M-0818-A7_RNA1_ref_supp -U ../results/host_filtered_fastq/1020-M-0818-A7_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A7_RNA1_ref_supp_self2_realign.bt.sam
bowtie2 -x 1428-M-3_RNA3_new_ref -U ../results/host_filtered_fastq/1428-M-3_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-3_RNA3_new_ref_self2_realign.bt.sam
bowtie2 -x 1428-M-42_RNA3_new_ref -U ../results/host_filtered_fastq/1428-M-42_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-42_RNA3_new_ref_self2_realign.bt.sam
bowtie2 -x 1428-M-43_Chaq_new_ref -U ../results/host_filtered_fastq/1428-M-43_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_Chaq_new_ref_self2_realign.bt.sam