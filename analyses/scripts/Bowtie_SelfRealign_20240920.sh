#!/usr/bin/env bash

# Use this script with the consensus sequences for galbut/chaq positive samples from geneious to capture all possible sequences


# Build bowtie indexes
bowtie2-build 500-F-21_RNA2.fasta 500-F-21_RNA2
bowtie2-build ME-M-7_RNA1.fasta ME-M-7_RNA1
bowtie2-build ME-F-1_Chaq.fasta ME-F-1_Chaq
bowtie2-build Penn-F-4_RNA2_2.fasta Penn-F-4_RNA2_2
bowtie2-build Penn-F-4_RNA2_2.fasta Penn-F-4_RNA2_2
bowtie2-build Peach-13-H7_RNA3_2.fasta Peach-13-H7_RNA3_2
bowtie2-build ME-F-1_RNA1.fasta ME-F-1_RNA1
bowtie2-build mel-20-TD-3_RNA3_1.fasta mel-20-TD-3_RNA3_1
bowtie2-build 500-F-41_RNA3.fasta 500-F-41_RNA3

# Run samples against new index
bowtie2 -x 500-F-21_RNA2 -U ../results/host_filtered_fastq/500-F-21_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-21_RNA2_self_realign1.bt.sam
bowtie2 -x ME-M-7_RNA1 -U ../results/host_filtered_fastq/ME-M-7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-7_RNA1_self_realign.bt.sam
bowtie2 -x ME-F-1_Chaq -U ../results/host_filtered_fastq/ME-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_Chaq_self_realign.bt.sam
bowtie2 -x Penn-F-4_RNA2_2 -U ../results/host_filtered_fastq/Penn-F-4_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-4_RNA2_2_self_realign.bt.sam
bowtie2 -x Penn-F-4_RNA2_2 -U ../results/host_filtered_fastq/Penn-F-4_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-4_RNA2_2_self_realign.bt.sam
bowtie2 -x Peach-13-H7_RNA3_2 -U ../results/host_filtered_fastq/Peach-13-H7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-13-H7_RNA3_2_self_realign.bt.sam
bowtie2 -x ME-F-1_RNA1 -U ../results/host_filtered_fastq/ME-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_RNA1_self_realign.bt.sam
bowtie2 -x mel-20-TD-3_RNA3_1-U ../results/host_filtered_fastq/mel-20-TD-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S mel-20-TD-3_RNA3_1_self_realign.bt.sam
bowtie2 -x 500-F-41_RNA3 -U ../results/host_filtered_fastq/500-F-41_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-41_RNA3_self_realign.bt.sam

bowtie2 -x 500-F-21_RNA2 -U ../results/host_filtered_fastq/500-F-21_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-21_RNA2_self_realign1.bt.sam
bowtie2 -x ME-M-7_RNA1 -U ../results/host_filtered_fastq/ME-M-7_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-7_RNA1_self_realign.bt.sam
bowtie2 -x ME-F-1_Chaq -U ../results/host_filtered_fastq/ME-F-1_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_Chaq_self_realign.bt.sam
bowtie2 -x Penn-F-4_RNA2_2 -U ../results/host_filtered_fastq/Penn-F-4_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-4_RNA2_2_self_realign.bt.sam
bowtie2 -x Penn-F-4_RNA2_2 -U ../results/host_filtered_fastq/Penn-F-4_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-4_RNA2_2_self_realign.bt.sam
bowtie2 -x Peach-13-H7_RNA3_2 -U ../results/host_filtered_fastq/Peach-13-H7_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-13-H7_RNA3_2_self_realign.bt.sam
bowtie2 -x ME-F-1_RNA1 -U ../results/host_filtered_fastq/ME-F-1_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_RNA1_self_realign.bt.sam
bowtie2 -x mel-20-TD-3_RNA3_1-U ../results/host_filtered_fastq/mel-20-TD-3_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S mel-20-TD-3_RNA3_1_self_realign.bt.sam
bowtie2 -x 500-F-41_RNA3 -U ../results/host_filtered_fastq/500-F-41_R2_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-41_RNA3_self_realign.bt.sam