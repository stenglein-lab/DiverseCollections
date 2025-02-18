#!/usr/bin/env bash

# Use this script with the consensus sequences for galbut/chaq positive samples from geneious to capture all possible sequences


# Build bowtie indexes
bowtie2-build sim-20-TD-4_RNA3.fasta sim-20-TD-4_RNA3
bowtie2-build Penn-M-15_RNA3.fasta Penn-M-15_RNA3
bowtie2-build Penn-M-15_RNA2.fasta Penn-M-15_RNA2
bowtie2-build Penn-M-15_RNA1.fasta Penn-M-15_RNA1
bowtie2-build Penn-M-14_RNA1.fasta Penn-M-14_RNA1
bowtie2-build Penn-M-11_RNA1.fasta Penn-M-11_RNA1
bowtie2-build Penn-M-10_RNA3.fasta Penn-M-10_RNA3
bowtie2-build Penn-M-10_RNA2.fasta Penn-M-10_RNA2
bowtie2-build Penn-M-10_RNA1.fasta Penn-M-10_RNA1
bowtie2-build Penn-M-10_chaq.fasta Penn-M-10_chaq
bowtie2-build Penn-M-8_RNA2.fasta Penn-M-8_RNA2
bowtie2-build Penn-M-2_chaq.fasta Penn-M-2_chaq
bowtie2-build Penn-F-9_RNA3.fasta Penn-F-9_RNA3
bowtie2-build Penn-F-9_RNA2.fasta Penn-F-9_RNA2
bowtie2-build Penn-F-9_RNA1.fasta Penn-F-9_RNA1
bowtie2-build Penn-F-9_chaq.fasta Penn-F-9_chaq
bowtie2-build Penn-F-6_RNA3.fasta Penn-F-6_RNA3
bowtie2-build Penn-F-6_RNA1.fasta Penn-F-6_RNA1
bowtie2-build Penn-F-4_RNA3.fasta Penn-F-4_RNA3
bowtie2-build Penn-F-4_RNA2.fasta Penn-F-4_RNA2
bowtie2-build Penn-F-3_RNA3.fasta Penn-F-3_RNA3
bowtie2-build Penn-F-3_RNA2.fasta Penn-F-3_RNA2
bowtie2-build Penn-F-2_RNA3.fasta Penn-F-2_RNA3
bowtie2-build Penn-F-2_RNA2.fasta Penn-F-2_RNA2
bowtie2-build Penn-F-2_RNA1.fasta Penn-F-2_RNA1
bowtie2-build Penn-F-2_chaq_b.fasta Penn-F-2_chaq_b
bowtie2-build Penn-F-2_chaq_a.fasta Penn-F-2_chaq_a
bowtie2-build Penn-F-1_RNA2.fasta Penn-F-1_RNA2
bowtie2-build Penn-F-1_RNA1.fasta Penn-F-1_RNA1
bowtie2-build Peach-13-H7_RNA3_b.fasta Peach-13-H7_RNA3_b
bowtie2-build Peach-13-H7_RNA3_a.fasta Peach-13-H7_RNA3_a
bowtie2-build Peach-9-F8_RNA2.fasta Peach-9-F8_RNA2
bowtie2-build Peach-9-F8_chaq.fasta Peach-9-F8_chaq
bowtie2-build Peach-5-B9_RNA3.fasta Peach-5-B9_RNA3
bowtie2-build Peach-5-B9_RNA2.fasta Peach-5-B9_RNA2
bowtie2-build Peach-5-B9_RNA1.fasta Peach-5-B9_RNA1
bowtie2-build Peach-5-B9_chaq.fasta Peach-5-B9_chaq
bowtie2-build Peach-4-B7_RNA3.fasta Peach-4-B7_RNA3
bowtie2-build Peach-4-B7_RNA2.fasta Peach-4-B7_RNA2
bowtie2-build Peach-4-B7_RNA1.fasta Peach-4-B7_RNA1
bowtie2-build Peach-3-B5_RNA3.fasta Peach-3-B5_RNA3
bowtie2-build Peach-3-B5_RNA2.fasta Peach-3-B5_RNA2
bowtie2-build Peach-3-B5_RNA1.fasta Peach-3-B5_RNA1
bowtie2-build OH-M-5_RNA2.fasta OH-M-5_RNA2
bowtie2-build ME-M-19_RNA3.fasta ME-M-19_RNA3
bowtie2-build ME-M-19_RNA1.fasta ME-M-19_RNA1
bowtie2-build ME-M-13_RNA3.fasta ME-M-13_RNA3
bowtie2-build ME-M-13_RNA2.fasta ME-M-13_RNA2
bowtie2-build ME-M-13_RNA1.fasta ME-M-13_RNA1
bowtie2-build ME-M-13_chaq.fasta ME-M-13_chaq
bowtie2-build ME-M-11_RNA3.fasta ME-M-11_RNA3
bowtie2-build ME-M-11_RNA2.fasta ME-M-11_RNA2
bowtie2-build ME-M-8_RNA1.fasta ME-M-8_RNA1
bowtie2-build ME-M-2_RNA2.fasta ME-M-2_RNA2
bowtie2-build ME-M-1_RNA3.fasta ME-M-1_RNA3
bowtie2-build ME-M-1_RNA2.fasta ME-M-1_RNA2
bowtie2-build ME-M-1_RNA1.fasta ME-M-1_RNA1
bowtie2-build mel-20-TD-3_RNA3.fasta mel-20-TD-3_RNA3
bowtie2-build ME-F-9_RNA1.fasta ME-F-9_RNA1
bowtie2-build ME-F-3_RNA2.fasta ME-F-3_RNA2
bowtie2-build ME-F-1_RNA1.fasta ME-F-1_RNA1
bowtie2-build ME-F-1_chaq.fasta ME-F-1_chaq
bowtie2-build CVID-1006-H2_RNA2.fasta CVID-1006-H2_RNA2
bowtie2-build CVID-0825-C10_RNA2.fasta CVID-0825-C10_RNA2
bowtie2-build Adobe-B8_RNA3.fasta Adobe-B8_RNA3
bowtie2-build Adobe-B8_RNA2.fasta Adobe-B8_RNA2
bowtie2-build Adobe-B8_RNA1.fasta Adobe-B8_RNA1
bowtie2-build 1428-M-43_RNA3.fasta 1428-M-43_RNA3
bowtie2-build 1428-M-43_RNA2.fasta 1428-M-43_RNA2
bowtie2-build 1428-M-43_RNA1.fasta 1428-M-43_RNA1
bowtie2-build 1428-M-43_chaq.fasta 1428-M-43_chaq
bowtie2-build 1428-M-42_RNA3.fasta 1428-M-42_RNA3
bowtie2-build 1428-M-42_RNA2.fasta 1428-M-42_RNA2
bowtie2-build 1428-M-42_RNA1.fasta 1428-M-42_RNA1
bowtie2-build 1428-M-42_chaq.fasta 1428-M-42_chaq
bowtie2-build 1428-M-6_RNA2.fasta 1428-M-6_RNA2
bowtie2-build 1428-M-6_RNA1.fasta 1428-M-6_RNA1
bowtie2-build 1428-M-4_RNA2.fasta 1428-M-4_RNA2
bowtie2-build 1428-M-3_RNA3.fasta 1428-M-3_RNA3
bowtie2-build 1020-M-0818-A10_RNA1.fasta 1020-M-0818-A10_RNA1
bowtie2-build 1020-M-0818-A5_RNA3.fasta 1020-M-0818-A5_RNA3
bowtie2-build 1020-M-0818-A5_RNA2.fasta 1020-M-0818-A5_RNA2
bowtie2-build 1020-M-0818-A5_RNA1.fasta 1020-M-0818-A5_RNA1
bowtie2-build 1020-M-0818-A5_chaq.fasta 1020-M-0818-A5_chaq
bowtie2-build 1020-M-19_RNA3.fasta 1020-M-19_RNA3
bowtie2-build 1020-M-19_RNA2.fasta 1020-M-19_RNA2
bowtie2-build 1020-M-19_RNA1.fasta 1020-M-19_RNA1
bowtie2-build 1020-M-19_chaq.fasta 1020-M-19_chaq
bowtie2-build 1020-F-31_RNA3.fasta 1020-F-31_RNA3
bowtie2-build 1020-F-31_RNA2.fasta 1020-F-31_RNA2
bowtie2-build 1020-F-31_chaq.fasta 1020-F-31_chaq
bowtie2-build 1020-F3-chaq.fasta 1020-F3-chaq
bowtie2-build 1020-F3_RNA3.fasta 1020-F3_RNA3
bowtie2-build 1020-F1_RNA1.fasta 1020-F1_RNA1
bowtie2-build 1020-D3_RNA1.fasta 1020-D3_RNA1
bowtie2-build 1020-0818-A1_RNA1.fasta 1020-0818-A1_RNA1
bowtie2-build 500-M-F10_RNA2.fasta 500-M-F10_RNA2
bowtie2-build 500-M-85_RNA1.fasta 500-M-85_RNA1
bowtie2-build 500-M-84_RNA3.fasta 500-M-84_RNA3
bowtie2-build 500-M-84_RNA2.fasta 500-M-84_RNA2
bowtie2-build 500-M-84_chaq.fasta 500-M-84_chaq
bowtie2-build 500-M-83_RNA2.fasta 500-M-83_RNA2
bowtie2-build 500-M-71-2_RNA3.fasta 500-M-71-2_RNA3
bowtie2-build 500-M-71-2_RNA2.fasta 500-M-71-2_RNA2
bowtie2-build 500-M-71-2_chaq.fasta 500-M-71-2_chaq
bowtie2-build 500-M-57_RNA2.fasta 500-M-57_RNA2
bowtie2-build 500-M-57_RNA1.fasta 500-M-57_RNA1
bowtie2-build 500-M-17_RNA3.fasta 500-M-17_RNA3
bowtie2-build 500-M-17-2_RNA1.fasta 500-M-17-2_RNA1
bowtie2-build 500-F-41_RNA3.fasta 500-F-41_RNA3
bowtie2-build 500-F-41_RNA2.fasta 500-F-41_RNA2
bowtie2-build 500-F-41_RNA1.fasta 500-F-41_RNA1
bowtie2-build 500-F-41_chaq.fasta 500-F-41_chaq
bowtie2-build 500-F-21_RNA2.fasta 500-F-21_RNA2
bowtie2-build 500-F-21_RNA1.fasta 500-F-21_RNA1



# Run samples against new index
bowtie2 -x sim-20-TD-4_RNA3 -U ../results/host_filtered_fastq/sim-20-TD-4_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S sim-20-TD-4_RNA3_self_realign1.bt.sam
bowtie2 -x Penn-M-15_RNA3 -U ../results/host_filtered_fastq/Penn-M-15_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-15_RNA3_self_realign.bt.sam
bowtie2 -x Penn-M-15_RNA2 -U ../results/host_filtered_fastq/Penn-M-15_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-15_RNA2_self_realign.bt.sam
bowtie2 -x Penn-M-15_RNA1 -U ../results/host_filtered_fastq/Penn-M-15_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-15_RNA1_self_realign.bt.sam
bowtie2 -x Penn-M-14_RNA1 -U ../results/host_filtered_fastq/Penn-M-14_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-14_RNA1_self_realign.bt.sam
bowtie2 -x Penn-M-11_RNA1 -U ../results/host_filtered_fastq/Penn-M-11_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-11_RNA1_self_realign.bt.sam
bowtie2 -x Penn-M-10_RNA3 -U ../results/host_filtered_fastq/Penn-M-10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-10_RNA3_self_realign.bt.sam
bowtie2 -x Penn-M-10_RNA2 -U ../results/host_filtered_fastq/Penn-M-10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-10_RNA2_self_realign.bt.sam
bowtie2 -x Penn-M-10_RNA1 -U ../results/host_filtered_fastq/Penn-M-10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-10_RNA1_self_realign.bt.sam
bowtie2 -x Penn-M-10_chaq -U ../results/host_filtered_fastq/Penn-M-10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-10_chaq_self_realign.bt.sam
bowtie2 -x Penn-M-8_RNA2 -U ../results/host_filtered_fastq/Penn-M-8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-8_RNA2_self_realign.bt.sam
bowtie2 -x Penn-M-2_chaq -U ../results/host_filtered_fastq/Penn-M-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-M-2_chaq_self_realign.bt.sam
bowtie2 -x Penn-F-9_RNA3 -U ../results/host_filtered_fastq/Penn-F-9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-9_RNA3_self_realign.bt.sam
bowtie2 -x Penn-F-9_RNA2 -U ../results/host_filtered_fastq/Penn-F-9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-9_RNA2_self_realign.bt.sam
bowtie2 -x Penn-F-9_RNA1 -U ../results/host_filtered_fastq/Penn-F-9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-9_RNA1_self_realign.bt.sam
bowtie2 -x Penn-F-9_chaq -U ../results/host_filtered_fastq/Penn-F-9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-9_chaq_self_realign.bt.sam
bowtie2 -x Penn-F-6_RNA3 -U ../results/host_filtered_fastq/Penn-F-6_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-6_RNA3_self_realign.bt.sam
bowtie2 -x Penn-F-6_RNA1 -U ../results/host_filtered_fastq/Penn-F-6_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-6_RNA1_self_realign.bt.sam
bowtie2 -x Penn-F-4_RNA3 -U ../results/host_filtered_fastq/Penn-F-4_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-4_RNA3_self_realign.bt.sam
bowtie2 -x Penn-F-4_RNA2 -U ../results/host_filtered_fastq/Penn-F-4_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-4_RNA2_self_realign.bt.sam
bowtie2 -x Penn-F-3_RNA3 -U ../results/host_filtered_fastq/Penn-F-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-3_RNA3_self_realign.bt.sam
bowtie2 -x Penn-F-3_RNA2 -U ../results/host_filtered_fastq/Penn-F-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-3_RNA2_self_realign.bt.sam
bowtie2 -x Penn-F-2_RNA3 -U ../results/host_filtered_fastq/Penn-F-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-2_RNA3_self_realign.bt.sam
bowtie2 -x Penn-F-2_RNA2 -U ../results/host_filtered_fastq/Penn-F-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-2_RNA2_self_realign.bt.sam
bowtie2 -x Penn-F-2_RNA1 -U ../results/host_filtered_fastq/Penn-F-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-2_RNA1_self_realign.bt.sam
bowtie2 -x Penn-F-2_chaq_b -U ../results/host_filtered_fastq/Penn-F-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-2_chaq_b_self_realign.bt.sam
bowtie2 -x Penn-F-2_chaq_a -U ../results/host_filtered_fastq/Penn-F-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-2_chaq_a_self_realign.bt.sam
bowtie2 -x Penn-F-1_RNA2 -U ../results/host_filtered_fastq/Penn-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-1_RNA2_self_realign.bt.sam
bowtie2 -x Penn-F-1_RNA1 -U ../results/host_filtered_fastq/Penn-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Penn-F-1_RNA1_self_realign.bt.sam
bowtie2 -x Peach-13-H7_RNA3_b -U ../results/host_filtered_fastq/Peach-13-H7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-13-H7_RNA3_b_self_realign.bt.sam
bowtie2 -x Peach-13-H7_RNA3_a -U ../results/host_filtered_fastq/Peach-13-H7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-13-H7_RNA3_a_self_realign.bt.sam
bowtie2 -x Peach-9-F8_RNA2 -U ../results/host_filtered_fastq/Peach-9-F8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-9-F8_RNA2_self_realign.bt.sam
bowtie2 -x Peach-9-F8_chaq -U ../results/host_filtered_fastq/Peach-9-F8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-9-F8_chaq_self_realign.bt.sam
bowtie2 -x Peach-5-B9_RNA3 -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_RNA3_self_realign.bt.sam
bowtie2 -x Peach-5-B9_RNA2 -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_RNA2_self_realign.bt.sam
bowtie2 -x Peach-5-B9_RNA1 -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_RNA1_self_realign.bt.sam
bowtie2 -x Peach-5-B9_chaq -U ../results/host_filtered_fastq/Peach-5-B9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-5-B9_chaq_self_realign.bt.sam
bowtie2 -x Peach-4-B7_RNA3 -U ../results/host_filtered_fastq/Peach-4-B7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-4-B7_RNA3_self_realign.bt.sam
bowtie2 -x Peach-4-B7_RNA2 -U ../results/host_filtered_fastq/Peach-4-B7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-4-B7_RNA2_self_realign.bt.sam
bowtie2 -x Peach-4-B7_RNA1 -U ../results/host_filtered_fastq/Peach-4-B7_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-4-B7_RNA1_self_realign.bt.sam
bowtie2 -x Peach-3-B5_RNA3 -U ../results/host_filtered_fastq/Peach-3-B5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-3-B5_RNA3_self_realign.bt.sam
bowtie2 -x Peach-3-B5_RNA2 -U ../results/host_filtered_fastq/Peach-3-B5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-3-B5_RNA2_self_realign.bt.sam
bowtie2 -x Peach-3-B5_RNA1 -U ../results/host_filtered_fastq/Peach-3-B5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Peach-3-B5_RNA1_self_realign.bt.sam
bowtie2 -x OH-M-5_RNA2 -U ../results/host_filtered_fastq/OH-M-5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S OH-M-5_RNA2_self_realign.bt.sam
bowtie2 -x ME-M-19_RNA3 -U ../results/host_filtered_fastq/ME-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-19_RNA3_self_realign.bt.sam
bowtie2 -x ME-M-19_RNA1_a -U ../results/host_filtered_fastq/ME-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-19_RNA1_self_realign.bt.sam
bowtie2 -x ME-M-13_RNA3 -U ../results/host_filtered_fastq/ME-M-13_R1_fuh.fastq  --local --very-sensitive --no-unal --threads 18 -S ME-M-13_RNA3_self_realign.bt.sam
bowtie2 -x ME-M-13_RNA2 -U ../results/host_filtered_fastq/ME-M-13_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-13_RNA2_self_realign.bt.sam
bowtie2 -x ME-M-13_RNA1 -U ../results/host_filtered_fastq/ME-M-13_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-13_RNA1_self_realign.bt.sam
bowtie2 -x ME-M-13_chaq -U ../results/host_filtered_fastq/ME-M-13_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-13_chaq_self_realign.bt.sam
bowtie2 -x ME-M-11_RNA3 -U ../results/host_filtered_fastq/ME-M-11_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-11_RNA3_self_realign.bt.sam
bowtie2 -x ME-M-11_RNA2 -U ../results/host_filtered_fastq/ME-M-11_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-11_RNA2_self_realign.bt.sam
bowtie2 -x ME-M-8_RNA1 -U ../results/host_filtered_fastq/ME-M-8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-8_RNA1_self_realign.bt.sam
bowtie2 -x ME-M-2_RNA2 -U ../results/host_filtered_fastq/ME-M-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-2_RNA2_self_realign.bt.sam
bowtie2 -x ME-M-1_RNA3 -U ../results/host_filtered_fastq/ME-M-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-1_RNA3_self_realign.bt.sam
bowtie2 -x ME-M-1_RNA2 -U ../results/host_filtered_fastq/ME-M-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-1_RNA2_self_realign.bt.sam
bowtie2 -x ME-M-1_RNA1 -U ../results/host_filtered_fastq/ME-M-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-M-1_RNA1_self_realign.bt.sam
bowtie2 -x mel-20-TD-3_RNA3 -U ../results/host_filtered_fastq/mel-20-TD-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S mel-20-TD-3_RNA3_self_realign.bt.sam
bowtie2 -x ME-F-9_RNA1 -U ../results/host_filtered_fastq/ME-F-9_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-9_RNA1_self_realign.bt.sam
bowtie2 -x ME-F-3_RNA2 -U ../results/host_filtered_fastq/ME-F-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-3_RNA2_self_realign.bt.sam
bowtie2 -x ME-F-1_RNA1 -U ../results/host_filtered_fastq/ME-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_RNA1_self_realign.bt.sam
bowtie2 -x ME-F-1_chaq -U ../results/host_filtered_fastq/ME-F-1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S ME-F-1_chaq_self_realign.bt.sam
bowtie2 -x CVID-1006-H2_RNA2 -U ../results/host_filtered_fastq/CVID-1006-H2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S CVID-1006-H2_RNA2_self_realign.bt.sam
bowtie2 -x CVID-0825-C10_RNA2 -U ../results/host_filtered_fastq/CVID-0825-C10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S CVID-0825-C10_RNA2_self_realign.bt.sam
bowtie2 -x Adobe-B8_RNA3 -U ../results/host_filtered_fastq/Adobe-B8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Adobe-B8_RNA3_self_realign.bt.sam
bowtie2 -x Adobe-B8_RNA2 -U ../results/host_filtered_fastq/Adobe-B8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Adobe-B8_RNA2_self_realign.bt.sam
bowtie2 -x Adobe-B8_RNA1 -U ../results/host_filtered_fastq/Adobe-B8_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S Adobe-B8_RNA1_self_realign.bt.sam
bowtie2 -x 1428-M-43_RNA3 -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_RNA3_self_realign.bt.sam
bowtie2 -x 1428-M-43_RNA2 -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_RNA2_self_realign.bt.sam
bowtie2 -x 1428-M-43_RNA1 -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_RNA1_self_realign.bt.sam
bowtie2 -x 1428-M-43_chaq -U ../results/host_filtered_fastq/1428-M-43_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-43_chaq_self_realign.bt.sam
bowtie2 -x 1428-M-42_RNA3 -U ../results/host_filtered_fastq/1428-M-42_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-42_RNA3_self_realign.bt.sam
bowtie2 -x 1428-M-42_RNA2 -U ../results/host_filtered_fastq/1428-M-42_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-42_RNA2_self_realign.bt.sam
bowtie2 -x 1428-M-42_RNA1 -U ../results/host_filtered_fastq/1428-M-42_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-42_RNA1_self_realign.bt.sam
bowtie2 -x 1428-M-42_chaq -U ../results/host_filtered_fastq/1428-M-42_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-42_chaq_self_realign.bt.sam
bowtie2 -x 1428-M-6_RNA2 -U ../results/host_filtered_fastq/1428-M-6_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-6_RNA2_self_realign.bt.sam
bowtie2 -x 1428-M-6_RNA1 -U ../results/host_filtered_fastq/1428-M-6_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-6_RNA1_self_realign.bt.sam
bowtie2 -x 1428-M-4_RNA2 -U ../results/host_filtered_fastq/1428-M-4_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-4_RNA2_self_realign.bt.sam
bowtie2 -x 1428-M-3_RNA3 -U ../results/host_filtered_fastq/1428-M-3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1428-M-3_RNA3_self_realign.bt.sam
bowtie2 -x 1020-M-0818-A10_RNA1 -U ../results/host_filtered_fastq/1020-M-0818-A10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A10_RNA1_self_realign.bt.sam
bowtie2 -x 1020-M-0818-A5_RNA3 -U ../results/host_filtered_fastq/1020-M-0818-A5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A5_RNA3_self_realign.bt.sam
bowtie2 -x 1020-M-0818-A5_RNA2 -U ../results/host_filtered_fastq/1020-M-0818-A5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A5_RNA2_self_realign.bt.sam
bowtie2 -x 1020-M-0818-A5_RNA1 -U ../results/host_filtered_fastq/1020-M-0818-A5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A5_RNA1_self_realign.bt.sam
bowtie2 -x 1020-M-0818-A5_chaq -U ../results/host_filtered_fastq/1020-M-0818-A5_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-0818-A5_chaq_self_realign.bt.sam
bowtie2 -x 1020-M-19_RNA3 -U ../results/host_filtered_fastq/1020-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-19_RNA3_self_realign.bt.sam
bowtie2 -x 1020-M-19_RNA2 -U ../results/host_filtered_fastq/1020-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-19_RNA2_self_realign.bt.sam
bowtie2 -x 1020-M-19_RNA1 -U ../results/host_filtered_fastq/1020-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-19_RNA1_self_realign.bt.sam
bowtie2 -x 1020-M-19_chaq -U ../results/host_filtered_fastq/1020-M-19_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-M-19_chaq_self_realign.bt.sam
bowtie2 -x 1020-F-31_RNA3 -U ../results/host_filtered_fastq/1020-F-31_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F-31_RNA3_self_realign.bt.sam
bowtie2 -x 1020-F-31_RNA2 -U ../results/host_filtered_fastq/1020-F-31_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F-31_RNA2_self_realign.bt.sam
bowtie2 -x 1020-F-31_chaq -U ../results/host_filtered_fastq/1020-F-31_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F-31_chaq_self_realign.bt.sam
bowtie2 -x 1020-F3-chaq -U ../results/host_filtered_fastq/1020-F3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3-chaq_self_realign.bt.sam
bowtie2 -x 1020-F3_RNA3 -U ../results/host_filtered_fastq/1020-F3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F3_RNA3_self_realign.bt.sam
bowtie2 -x 1020-F1_RNA1 -U ../results/host_filtered_fastq/1020-F1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-F1_RNA1_self_realign.bt.sam
bowtie2 -x 1020-D3_RNA1 -U ../results/host_filtered_fastq/1020-D3_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-D3_RNA1_self_realign.bt.sam
bowtie2 -x 1020-0818-A1_RNA1 -U ../results/host_filtered_fastq/1020-0818-A1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 1020-0818-A1_RNA1_self_realign.bt.sam
bowtie2 -x 500-M-F10_RNA2 -U ../results/host_filtered_fastq/500-M-F10_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-F10_RNA2_self_realign.bt.sam
bowtie2 -x 500-M-85_RNA1 -U ../results/host_filtered_fastq/500-M-85_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-85_RNA1_self_realign.bt.sam
bowtie2 -x 500-M-84_RNA3 -U ../results/host_filtered_fastq/500-M-84_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-84_RNA3_self_realign.bt.sam
bowtie2 -x 500-M-84_RNA2 -U ../results/host_filtered_fastq/500-M-84_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-84_RNA2_self_realign.bt.sam
bowtie2 -x 500-M-84_chaq -U ../results/host_filtered_fastq/500-M-84_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-84_chaq_self_realign.bt.sam
bowtie2 -x 500-M-83_RNA2 -U ../results/host_filtered_fastq/500-M-83_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-83_RNA2_self_realign.bt.sam
bowtie2 -x 500-M-71-2_RNA3 -U ../results/host_filtered_fastq/500-M-71-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-71-2_RNA3_self_realign.bt.sam
bowtie2 -x 500-M-71-2_RNA2 -U ../results/host_filtered_fastq/500-M-71-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-71-2_RNA2_self_realign.bt.sam
bowtie2 -x 500-M-71-2_chaq -U ../results/host_filtered_fastq/500-M-71-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-71-2_chaq_self_realign.bt.sam
bowtie2 -x 500-M-57_RNA2 -U ../results/host_filtered_fastq/500-M-57_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-57_RNA2_self_realign.bt.sam
bowtie2 -x 500-M-57_RNA1 -U ../results/host_filtered_fastq/500-M-57_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-57_RNA1_self_realign.bt.sam
bowtie2 -x 500-M-17_RNA3 -U ../results/host_filtered_fastq/500-M-17_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-17_RNA3_self_realign.bt.sam
bowtie2 -x 500-M-17-2_RNA1 -U ../results/host_filtered_fastq/500-M-17-2_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-M-17-2_RNA11_self_realign.bt.sam
bowtie2 -x 500-F-41_RNA3 -U ../results/host_filtered_fastq/500-F-41_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-41_RNA3_self_realign.bt.sam
bowtie2 -x 500-F-41_RNA2 -U ../results/host_filtered_fastq/500-F-41_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-41_RNA2_self_realign.bt.sam
bowtie2 -x 500-F-41_RNA1 -U ../results/host_filtered_fastq/500-F-41_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-41_RNA1_self_realign.bt.sam
bowtie2 -x 500-F-41_chaq -U ../results/host_filtered_fastq/500-F-41_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-41_chaq_self_realign.bt.sam
bowtie2 -x 500-F-21_RNA2 -U ../results/host_filtered_fastq/500-F-21_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-21_RNA2_self_realign.bt.sam
bowtie2 -x 500-F-21_RNA1 -U ../results/host_filtered_fastq/1020-0818-A1_R1_fuh.fastq --local --very-sensitive --no-unal --threads 18 -S 500-F-21_RNA1_self_realign.bt.sam