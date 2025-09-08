library(phytools)
library(ape)
library(tidyverse)
library(readxl)
library(TreeTools)

# Read in trees
RNA1 <- read.tree("analyses/trees/RNA1/RNA1_nucleotide_alignment.fasta-1.contree")
RNA2 <- read.tree("analyses/trees/RNA2/RNA2_cont_removed_align.fasta.contree")
RNA3 <- read.tree("analyses/trees/RNA3/RNA3_cont_removed_align.fasta.contree")
Chaq <- read.tree("analyses/trees/Chaq/chaq_cont_removed_align.fasta.contree")

# Midpoint root trees
rna1_tree_rooted <- midpoint_root(RNA1)
rna2_tree_rooted <- midpoint_root(RNA2)
rna3_tree_rooted <- midpoint_root(RNA3)
Chaq_tree_rooted <- midpoint_root(Chaq)

# write out midpoint rooted trees in .newick format
write.tree(rna1_tree_rooted, file="analyses/trees/RNA1/rna1_tree_rooted.from_phytools.newick")
write.tree(rna2_tree_rooted, file="analyses/trees/RNA2/rna2_tree_rooted.from_phytools.newick")
write.tree(rna3_tree_rooted, file="analyses/trees/RNA3/rna3_tree_rooted.from_phytools.newick")
write.tree(Chaq_tree_rooted, file="analyses/trees/Chaq/chaq_tree_rooted.from_phytools.newick")
