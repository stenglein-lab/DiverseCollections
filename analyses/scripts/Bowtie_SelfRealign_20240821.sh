#!/usr/bin/env bash

# Use this script with the consensus sequences for galbut/chaq positive samples from geneious to capture all possible sequences


# Build bowtie indexes
bowtie2-build 500-F-21_RNA1.fasta 500-F-21_RNA1
bowtie2-build 500-F-41_RNA3.fasta 500-F-41_RNA3
bowtie2-build 500-M-17_RNA3.fasta 500-M-17_RNA3
bowtie2-build 500-M-60_RNA1.fasta 500-M-60_RNA1
bowtie2-build 500-M-71-2_RNA1.fasta 500-M-71-2_RNA1
bowtie2-build 500-M-84_RNA1.fasta 500-M-84_RNA1
bowtie2-build 500-M-F11_RNA1.fasta 500-M-F11_RNA1
bowtie2-build 1020-A7_RNA2.fasta 1020-A7_RNA2
bowtie2-build 1020-D3_RNA1.fasta 1020-D3_RNA1
bowtie2-build 1020-F1_RNA1.fasta 1020-F1_RNA1
bowtie2-build 1020-F3_RNA2.fasta 1020-F3_RNA2
bowtie2-build 1020-F3_Chaq.fasta 1020-F3_Chaq
bowtie2-build 1020-F3_RNA3.fasta 1020-F3_RNA3
bowtie2-build 1020-M-19_RNA3.fasta 1020-M-19_RNA3
bowtie2-build 1020-M-0818-A7_RNA1.fasta 1020-M-0818-A7_RNA1
bowtie2-build 1428-M-3_RNA3.fasta 1428-M-3_RNA3
bowtie2-build 1428-M-6_RNA1.fasta 1428-M-6_RNA1
bowtie2-build 1428-M-6_RNA2.fasta 1428-M-6_RNA2
bowtie2-build 1428-M-42_RNA3.fasta 1428-M-42_RNA3
bowtie2-build 1428-M-43_Chaq.fasta 1428-M-43_Chaq
bowtie2-build 1428-M-43_RNA1.fasta 1428-M-43_RNA1
bowtie2-build 1428-M-43_RNA2.fasta 1428-M-43_RNA2
bowtie2-build 1428-M-43_RNA3.fasta 1428-M-43_RNA3
bowtie2-build CVID-0825-B1_RNA1.fasta CVID-0825-B1_RNA1
bowtie2-build ME-F-1_Chaq.fasta ME-F-1_Chaq
bowtie2-build ME-F-1_RNA1.fasta ME-F-1_RNA1
bowtie2-build ME-F-2_RNA1.fasta ME-F-2_RNA1
bowtie2-build ME-F-9_RNA1.fasta ME-F-9_RNA1
bowtie2-build ME-M-8_RNA1.fasta ME-M-8_RNA1
bowtie2-build ME-M-10_RNA2.fasta ME-M-10_RNA2
bowtie2-build ME-M-19_RNA1.fasta ME-M-19_RNA1
bowtie2-build ME-M-19_RNA3.fasta ME-M-19_RNA3
bowtie2-build Peach-3-B5_RNA1.fasta Peach-3-B5_RNA1
bowtie2-build Peach-3-B5_RNA2.fasta Peach-3-B5_RNA2
bowtie2-build Peach-3-B5_RNA3.fasta Peach-3-B5_RNA3
bowtie2-build Peach-4-B7_RNA1.fasta Peach-4-B7_RNA1
bowtie2-build Peach-4-B7_RNA2.fasta Peach-4-B7_RNA2
bowtie2-build Peach-4-B7_RNA3.fasta Peach-4-B7_RNA3
bowtie2-build Peach-5-B9_Chaq.fasta Peach-5-B9_Chaq
bowtie2-build Peach-5-B9_RNA1.fasta Peach-5-B9_RNA1
bowtie2-build Peach-5-B9_RNA2.fasta Peach-5-B9_RNA2
bowtie2-build Peach-5-B9_RNA3.fasta Peach-5-B9_RNA3
bowtie2-build Peach-9-F8_Chaq.fasta Peach-9-F8_Chaq
bowtie2-build Peach-9-F8_RNA2.fasta Peach-9-F8_RNA2
bowtie2-build Penn-F-1_RNA1.fasta Penn-F-1_RNA1
bowtie2-build Penn-F-1_RNA2.fasta Penn-F-1_RNA2
bowtie2-build Penn-F-6_RNA1.fasta Penn-F-6_RNA1
bowtie2-build Penn-F-6_RNA3.fasta Penn-F-6_RNA3
bowtie2-build Penn-F-9_RNA1.fasta Penn-F-9_RNA1
bowtie2-build Penn-M-11_RNA1.fasta Penn-M-11_RNA1
bowtie2-build Penn-M-14_RNA1.fasta Penn-M-14_RNA1


# Run samples against new index
bowtie2 -x 500-F-21_RNA1 -U ../results/host_filtered_fastq/500-F-21_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-21_RNA1_self_realign1.bt.sam
bowtie2 -x 500-F-41_RNA3 -U ../results/host_filtered_fastq/500-F-41_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-41_RNA3_self_realign.bt.sam
bowtie2 -x 500-M-17_RNA3 -U ../results/host_filtered_fastq/500-M-17_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-17_RNA3_self_realign.bt.sam
bowtie2 -x 500-M-60_RNA1 -U ../results/host_filtered_fastq/500-M-60_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-60_RNA1_self_realign.bt.sam
bowtie2 -x 500-M-71-2_RNA1 -U ../results/host_filtered_fastq/500-M-71-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-71-2_RNA1_self_realign.bt.sam
bowtie2 -x 500-M-84_RNA1 -U ../results/host_filtered_fastq/500-M-84_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-84_RNA1_self_realign.bt.sam
bowtie2 -x 500-M-F11_RNA1 -U ../results/host_filtered_fastq/500-M-F11_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-F11_RNA1_self_realign.bt.sam
bowtie2 -x 1020-A7_RNA2 -U ../results/host_filtered_fastq/1020-A7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-A7_RNA2_self_realign.bt.sam
bowtie2 -x 1020-D3_RNA1 -U ../results/host_filtered_fastq/1020-D3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-D3_RNA1_self_realign.bt.sam
bowtie2 -x 1020-F1_RNA1 -U ../results/host_filtered_fastq/1020-F1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F1_RNA1self_realign.bt.sam
bowtie2 -x 1020-F3_Chaq -U ../results/host_filtered_fastq/1020-F3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_Chaq_self_realign.bt.sam
bowtie2 -x 1020-F3_RNA2 -U ../results/host_filtered_fastq/1020-F3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_RNA2_self_realign.bt.sam
bowtie2 -x 1020-F3_RNA3 -U ../results/host_filtered_fastq/1020-F3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_RNA3_self_realign.bt.sam
bowtie2 -x 1020-M-19_RNA3 -U ../results/host_filtered_fastq/1020-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-19_RNA3_self_realign.bt.sam
bowtie2 -x 1020-M-0818-A7_RNA1 -U ../results/host_filtered_fastq/1020-M-0818-A7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A7_RNA1_self_realign.bt.sam
bowtie2 -x 1428-M-3_RNA3 -U ../results/host_filtered_fastq/1428-M-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-3_RNA3_self_realign.bt.sam
bowtie2 -x 1428-M-6_RNA1 -U ../results/host_filtered_fastq/1428-M-6_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-6_RNA1_self_realign.bt.sam
bowtie2 -x 1428-M-6_RNA2 -U ../results/host_filtered_fastq/1428-M-6_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-6_RNA2_self_realign.bt.sam
bowtie2 -x 1428-M-42_RNA3 -U ../results/host_filtered_fastq/1428-M-42_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-42_RNA3_self_realign.bt.sam
bowtie2 -x 1428-M-43_Chaq -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_Chaq_self_realign.bt.sam
bowtie2 -x 1428-M-43_RNA1 -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_RNA1_self_realign.bt.sam
bowtie2 -x 1428-M-43_RNA2 -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_RNA2_self_realign.bt.sam
bowtie2 -x 1428-M-43_RNA3 -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_RNA3_self_realign.bt.sam
bowtie2 -x CVID-0825-B1_RNA1 -U ../results/host_filtered_fastq/CVID-0825-B1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S CVID-0825-B1_RNA1_self_realign.bt.sam
bowtie2 -x ME-F-1_Chaq -U ../results/host_filtered_fastq/ME-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_Chaq_self_realign.bt.sam
bowtie2 -x ME-F-1_RNA1 -U ../results/host_filtered_fastq/ME-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_RNA1_self_realign.bt.sam
bowtie2 -x ME-F-2_RNA1 -U ../results/host_filtered_fastq/ME-F-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-2_RNA1_self_realign.bt.sam
bowtie2 -x ME-F-9_RNA1 -U ../results/host_filtered_fastq/ME-F-9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-9_RNA1_self_realign.bt.sam
bowtie2 -x ME-M-8_RNA1 -U ../results/host_filtered_fastq/ME-M-8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-8_RNA1_self_realign.bt.sam
bowtie2 -x ME-M-10_RNA2 -U ../results/host_filtered_fastq/ME-M-10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-10_RNA2_self_realign.bt.sam
bowtie2 -x ME-M-19_RNA1 -U ../results/host_filtered_fastq/ME-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-19_RNA1_self_realign.bt.sam
bowtie2 -x ME-M-19_RNA3 -U ../results/host_filtered_fastq/ME-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-19_RNA3_self_realign.bt.sam
bowtie2 -x Peach-3-B5_RNA1 -U ../results/host_filtered_fastq/Peach-3-B5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-3-B5_RNA1_self_realign.bt.sam
bowtie2 -x Peach-3-B5_RNA2 -U ../results/host_filtered_fastq/Peach-3-B5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-3-B5_RNA2_self_realign.bt.sam
bowtie2 -x Peach-3-B5_RNA3 -U ../results/host_filtered_fastq/Peach-3-B5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-3-B5_RNA3_self_realign.bt.sam
bowtie2 -x Peach-4-B7_RNA1 -U ../results/host_filtered_fastq/Peach-4-B7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-4-B7_RNA1_self_realign.bt.sam
bowtie2 -x Peach-4-B7_RNA2 -U ../results/host_filtered_fastq/Peach-4-B7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-4-B7_RNA2_self_realign.bt.sam
bowtie2 -x Peach-4-B7_RNA3 -U ../results/host_filtered_fastq/Peach-4-B7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-4-B7_RNA3_self_realign.bt.sam
bowtie2 -x Peach-5-B9_Chaq -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_Chaq_self_realign.bt.sam
bowtie2 -x Peach-5-B9_RNA1 -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_RNA1_self_realign.bt.sam
bowtie2 -x Peach-5-B9_RNA2 -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_RNA2_self_realign.bt.sam
bowtie2 -x Peach-5-B9_RNA3 -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_RNA3_self_realign.bt.sam
bowtie2 -x Peach-9-F8_Chaq -U ../results/host_filtered_fastq/Peach-9-F8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-9-F8_Chaq_self_realign.bt.sam
bowtie2 -x Peach-9-F8_RNA2 -U ../results/host_filtered_fastq/Peach-9-F8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-9-F8_RNA2_self_realign.bt.sam
bowtie2 -x Penn-F-1_RNA1 -U ../results/host_filtered_fastq/Penn-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-1_RNA1_self_realign.bt.sam
bowtie2 -x Penn-F-1_RNA2 -U ../results/host_filtered_fastq/Penn-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-1_RNA2_self_realign.bt.sam
bowtie2 -x Penn-F-6_RNA1 -U ../results/host_filtered_fastq/Penn-F-6_R1_fuh.fastq  --local --very-sensitive --no-unal --threads 18 -S Penn-F-6_RNA1_self_realign.bt.sam
bowtie2 -x Penn-F-6_RNA3 -U ../results/host_filtered_fastq/Penn-F-6_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-6_RNA3_self_realign.bt.sam
bowtie2 -x Penn-F-9_RNA1 -U ../results/host_filtered_fastq/Penn-F-9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-9_RNA1_self_realign.bt.sam
bowtie2 -x Penn-M-11_RNA1 -U ../results/host_filtered_fastq/Penn-M-11_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-11_RNA1_self_realign.bt.sam
bowtie2 -x Penn-M-14_RNA1 -U ../results/host_filtered_fastq/Penn-M-14_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-14_RNA1_self_realign.bt.sam