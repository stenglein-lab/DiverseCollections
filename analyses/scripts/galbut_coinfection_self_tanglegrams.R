library(phytools)
library(ape)
library(tidyverse)
library(readxl)
# library(forcats)

metadata <- read_xlsx("analyses/data/Metadata_Table.xlsx")

# make a map of sample IDs -> accessions
id_acc_map <- metadata %>% select(sample_name, accession)
id_acc_map$accession <- str_replace_all(id_acc_map$accession, " ", "")
id_acc_map <- id_acc_map %>% separate_longer_delim(accession, delim=",") %>% filter(!is.na(accession))

# Read in trees
RNA1 <- read.tree("analyses/trees/RNA1/RNA1_nucleotide_alignment.fasta.contree")
RNA2 <- read.tree("analyses/trees/RNA2/RNA2_nucleotide_alignment.fasta.contree")
RNA3 <- read.tree("analyses/trees/RNA3/RNA3_nucelotide_alignment.fasta.contree")
Chaq <- read.tree("analyses/trees/Chaq/Chaq_nucleotide_alignment.fasta.contree")

# midpoint root both trees
rna1_tree_rooted <- midpoint_root(RNA1)
rna2_tree_rooted <- midpoint_root(RNA2)
rna3_tree_rooted <- midpoint_root(RNA3)
Chaq_tree_rooted <- midpoint_root(Chaq)

# Read in sequence ids
RNA1 <- read_xlsx("analyses/trees/RNA1/RNA1_Seq_IDs.xlsx")
RNA2 <- read_xlsx("analyses/trees/RNA2/RNA2_Seq_IDs.xlsx")
RNA3 <- read_xlsx("analyses/trees/RNA3/RNA3_SeqIDs.xlsx")
Chaq <- read_xlsx("analyses/trees/Chaq/Chaq_Seq_IDs.xlsx")

RNA1$sample_id <- paste(RNA1$location, RNA1$date, RNA1$accession, sep = "_")
RNA1_ids <- RNA1 %>% 
  select(accession, sample_id)

RNA2$sample_id <- paste(RNA2$location, RNA2$date, RNA2$accession, sep = "_")
RNA2_ids <- RNA2 %>% 
  select(accession, sample_id)

RNA3$sample_id <- paste(RNA3$location, RNA3$date, RNA3$accession, sep = "_")
RNA3_ids <- RNA3 %>% 
  select(accession, sample_id)

Chaq$sample_id <- paste(Chaq$location, Chaq$date, Chaq$accession, sep = "_")
Chaq_ids <- Chaq %>% 
  select(accession, sample_id)

# relabel tips on trees
old_tip_labels_1 <- RNA1_ids$accession
new_tip_labels_1 <- RNA1_ids$sample_id

old_tip_labels_2 <- RNA2_ids$accession
new_tip_labels_2 <- RNA2_ids$sample_id

old_tip_labels_3 <- RNA3_ids$accession
new_tip_labels_3 <- RNA3_ids$sample_id

old_tip_labels_c <- Chaq_ids$accession
new_tip_labels_c <- Chaq_ids$sample_id

# relabel tips
rna1_tree_relabeled <- updateLabel(rna1_tree_rooted, old_tip_labels_1, new_tip_labels_1)
rna2_tree_relabeled <- updateLabel(rna2_tree_rooted, old_tip_labels_2, new_tip_labels_2)
rna3_tree_relabeled <- updateLabel(rna3_tree_rooted, old_tip_labels_3, new_tip_labels_3)
Chaq_tree_relabeled <- updateLabel(Chaq_tree_rooted, old_tip_labels_c, new_tip_labels_c)

new_tip_labels_1

# -----------------------------------------
# make self-tanglegram co-infection figure 
# -----------------------------------------

# link lines should connect sequences from the same sample (but not the same sequence)
RNA1_ids <- left_join(RNA1_ids, id_acc_map)
RNA2_ids <- left_join(RNA2_ids, id_acc_map)
RNA3_ids <- left_join(RNA3_ids, id_acc_map)
Chaq_ids <- left_join(Chaq_ids, id_acc_map)

make_self_tanglegram <- function(tree, ids, label) {
  
  coinf_ids  <- ids %>% filter(!is.na(sample_name)) %>% group_by(sample_name) %>% mutate(n_seg = n(), tip_n = row_number()) %>% filter(n_seg > 1)
  coinf_link <- coinf_ids %>% select(sample_id, tip_n) %>% pivot_wider(names_from = tip_n, names_prefix = "tip_", values_from = sample_id) %>% select(-sample_name)
  assoc <- cbind(coinf_link$tip_1, coinf_link$tip_2)
  
  # make cophylogeny (this makes the object but doesn't plot it yet)
  cophy <- cophylo(tree, tree, assoc = assoc, rotate=F)
  
  # will save as PDF
  st_filename=paste0("analyses/trees/tanglegrams/", label, "_self_tanglegram.pdf")
  pdf(file = st_filename, width=8.5, height=11)
  
  # plot tanglegram
  plot.cophylo(cophy, 
               link.type = "curved", 
               link.lwd=1,
               link.lty="solid",
               link.col=make.transparent("coral4", 0.5),
               fsize=0.3, 
               pts=F,
               scale.bar = c(0.025, 0.025))
  
  # this doesn't plot node support values, but it could.
  # see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

  # turn off PDF
  dev.off()
}

make_self_tanglegram(rna1_tree_relabeled, RNA1_ids, "rna1")
make_self_tanglegram(rna2_tree_relabeled, RNA2_ids, "rna2")
make_self_tanglegram(rna3_tree_relabeled, RNA3_ids, "rna3")
make_self_tanglegram(Chaq_tree_relabeled, Chaq_ids, "chaq")



