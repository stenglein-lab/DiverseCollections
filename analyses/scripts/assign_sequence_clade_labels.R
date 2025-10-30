library(tidyverse)
library(ape)
library(phytools)

# this script assigns clade labels to all galbut virus acccessions 
# namely, major clade / genotype labels like melA or simA or melB
#
# It does this using trees of the sequences and, for each tree and each clade,
# it uses 2 accessions whose MRCA defines the clade.  In other words, 2 accessions
# whose MRCA is the node that is the common ancestor of all sequences in the clade.
# 
# Then it assigns the same clade label to all tips downstream of that MRCA.
#
# This is a data-driven way to map accessions to clades that doesn't rely on 
# possibly error-prone metadata sheets
#
# MDS 9/17/2025

# Read in iqtree trees
RNA1 <- read.tree("analyses/trees/RNA1/RNA1_nucleotide_alignment.fasta-1.contree")
RNA2 <- read.tree("analyses/trees/RNA2/RNA2_cont_removed_align.fasta.contree")
RNA3 <- read.tree("analyses/trees/RNA3/RNA3_cont_removed_align.fasta.contree")
Chaq <- read.tree("analyses/trees/Chaq/chaq_cont_removed_align.fasta.contree")

# handle input trees:
# midpoint root, and 
# create unique node labels (instead of support values, which is what iqtree trees uses for node labels)
handle_tree <- function(tree) {
  
  # Midpoint root tree
  tree_rooted <- midpoint_root(tree)
  
  # node labels in trees from iqtree are support values
  # replace these values with a unique label for each node
  tree_rooted$node.label <- paste0("node_", seq(1:length(tree_rooted$node.label)))
  
  tree_rooted
}

# process input trees
rna1_tree <- handle_tree(RNA1)
rna2_tree <- handle_tree(RNA2)
rna3_tree <- handle_tree(RNA3)
chaq_tree <- handle_tree(Chaq)

# get node label of a clade using the getMRCA function, which returns a node index not a label
get_mrca_node_label <- function(tree, tips) {
  # get the node label for the MRCA node of 2 tips
  # see: https://stackoverflow.com/questions/51696837/r-phylo-object-how-to-connect-node-label-and-node-number
  node_index <- getMRCA(tree, tip = tips)
  tree$node.label[node_index-Ntip(tree)]
}

# this function takes as input:
# - a tree
# - tip labels that will define a clade of interest MRCA node
# - a label to attach to the tips within that clade
# it returns a tibble with 
# - a column of tip labels
# - a column with a clade label to associate with those tips (1 value)
label_clade_tips <- function(tree, mrca_defining_tips, label) {
  mrca_node_label <- get_mrca_node_label(tree, mrca_defining_tips)
  clade           <- extract.clade(tree, mrca_node_label)
  clade_tips      <- clade$tip.label 
  df <- tibble(accession   = clade_tips, 
               clade_label = label)
  df
}

# RNA 1 
plot.phylo(rna1_tree, show.node.label = T, cex=0.5)
rna1_melA  <- label_clade_tips(rna1_tree, c("galbut_EVE_from_Wallace_supp_fig_S10", "PQ625069"), "melA")
rna1_simA  <- label_clade_tips(rna1_tree, c("MW976829", "MW976831"), "simA")
rna1_melB  <- label_clade_tips(rna1_tree, c("PQ625078", "PQ625064"), "melB")
# double check that we've accounted for all tips with our MRCA strategy
stopifnot(nrow(rna1_melA) + nrow(rna1_simA) + nrow(rna1_melB) == length(rna1_tree$tip.label)) 

# RNA 2
plot.phylo(rna2_tree, show.node.label = T, cex=0.5)
rna2_melA  <- label_clade_tips(rna2_tree, c("PQ624966", "OR729044"), "melA")
rna2_simA  <- label_clade_tips(rna2_tree, c("OR820595", "SRX10234780_RNA2"), "simA")
rna2_melB  <- label_clade_tips(rna2_tree, c("PQ624975", "PQ624997"), "melB")
# double check that we've accounted for all tips with our MRCA strategy
stopifnot(nrow(rna2_melA) + nrow(rna2_simA) + nrow(rna2_melB) == length(rna2_tree$tip.label)) 

# RNA3
plot.phylo(rna3_tree, show.node.label = T, cex=0.5)
rna3_melA  <- label_clade_tips(rna3_tree, c("PQ624895", "PQ624857"), "melA")
rna3_simA  <- label_clade_tips(rna3_tree, c("OR820601", "SRX10234787_RNA3"), "simA")
rna3_melB  <- label_clade_tips(rna3_tree, c("PQ624866", "PQ624913"), "melB")
# double check that we've accounted for all tips with our MRCA strategy
stopifnot(nrow(rna3_melA) + nrow(rna3_simA) + nrow(rna3_melB) == length(rna3_tree$tip.label)) 

# Chaq
plot.phylo(chaq_tree, show.node.label = T, cex=0.5)
chaq_melA  <- label_clade_tips(chaq_tree, c("PQ625188", "PQ625178"), "melA")
chaq_simA  <- label_clade_tips(chaq_tree, c("OR820594", "MW976841"), "simA")
# double check that we've accounted for all tips with our MRCA strategy
stopifnot(nrow(chaq_melA) + nrow(chaq_simA) == length(chaq_tree$tip.label)) 

# create a combined table of accessions -> clade label for all accessions
combined_clade_labels <- 
  rbind(
    rna1_melA, rna1_melB, rna1_simA,
    rna2_melA, rna2_melB, rna2_simA,
    rna3_melA, rna3_melB, rna3_simA,
    chaq_melA,            chaq_simA
  )

# this table maps accessions -> clades (melA, melB, etc)
write.table(combined_clade_labels, 
            "analyses/data/accession_clade_labels.txt",
            quote=F,
            sep="\t",
            row.names = F)

