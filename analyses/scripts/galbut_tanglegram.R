library(phytools)
library(ape)
library(tidyverse)
library(readxl)

# a script to create a galbut virus RNA1/RNA2 tanglegram
# MDS 12/19/2024

# RNA1 RNA2
# read in 2 individual trees
# use ape read.tree
RNA1 <- read.tree("analyses/trees/RNA1/RNA1_nucleotide_alignment.fasta.contree")
RNA2 <- read.tree("analyses/trees/RNA2/RNA2_nucleotide_alignment.fasta.contree")
RNA3 <- read.tree("analyses/trees/RNA3/RNA3_nucelotide_alignment.fasta.contree")
Chaq <- read.tree("analyses/trees/Chaq/Chaq_alignment_recomb_removed.fasta.contree")

# midpoint root both trees
rna1_tree_rooted <- midpoint_root(RNA1)
rna2_tree_rooted <- midpoint_root(RNA2)
rna3_tree_rooted <- midpoint_root(RNA3)
Chaq_tree_rooted <- midpoint_root(Chaq)

# Read in sequence ids
RNA1_id <- read_xlsx("analyses/trees/RNA1/RNA1_Seq_IDs.xlsx")
RNA2_id <- read_xlsx("analyses/trees/RNA2/RNA2_Seq_IDs.xlsx")
RNA3_id <- read_xlsx("analyses/trees/RNA3/RNA3_SeqIDs.xlsx")
Chaq_id <- read_xlsx("analyses/trees/Chaq/Chaq_Seq_IDs.xlsx")

# Make ID that is the same as individual trees
RNA1_id <- RNA1_id %>% 
  unite(id, c(location, date, accession), sep = "_", remove = FALSE) %>% 
  separate(isolate, c("D", "spp", "country", "year", "isolate", "virus", "virus_2", 
                      "segment", "CDS", "CDS_2"), sep = "_", remove = FALSE) %>% 
  mutate(isolate = coalesce(isolate, accession)) %>% 
  select(isolate, id, accession)

RNA2_id <- RNA2_id %>% 
  unite(id, c(location, date, accession), sep = "_", remove = FALSE)%>% 
  separate(isolate, c("D", "spp", "country", "year", "isolate", "virus", "virus_2", 
                      "segment", "CDS", "CDS_2"), sep = "_", remove = FALSE) %>% 
  mutate(isolate = coalesce(isolate, accession)) %>% 
  select(isolate, id, accession)

RNA3_id <- RNA3_id %>% 
  unite(id, c(location, date, accession), sep = "_", remove = FALSE)%>% 
  separate(isolate, c("D", "spp", "country", "year", "isolate", "virus", "virus_2", 
                      "segment", "CDS", "CDS_2"), sep = "_", remove = FALSE) %>% 
  mutate(isolate = coalesce(isolate, accession)) %>% 
  select(isolate, id, accession)

Chaq_id <- Chaq_id %>% 
  unite(id, c(location, date, accession), sep = "_", remove = FALSE) %>% 
  separate(isolate, c("D", "spp", "country", "year", "isolate", "virus", "virus_2", 
                      "CDS", "CDS_2"), sep = "_", remove = FALSE) %>% 
  mutate(isolate = coalesce(isolate, accession)) %>% 
  select(isolate, id, accession)

# Manually fix a few names (EVE, SRA seqs, available seqs)
write.csv(RNA1_id, file = "analyses/trees/RNA1_newnames.csv")
RNA1_id <- read_csv("analyses/trees/RNA1_newnames.csv") %>% 
  rename("RNA1" = id)

write.csv(RNA2_id, file = "analyses/trees/RNA2_newnames.csv")
RNA2_id <- read_csv("analyses/trees/RNA2_newnames.csv") %>% 
  rename("RNA2" = id)


write.csv(RNA3_id, file = "analyses/trees/RNA3_newnames.csv")
RNA3_id <- read_csv("analyses/trees/RNA3_newnames.csv") %>% 
  rename("RNA3" = id)

write.csv(Chaq_id, file = "analyses/trees/Chaq_newnames.csv")
Chaq_id <- read_csv("analyses/trees/Chaq_newnames.csv") %>% 
  rename("Chaq" = id)

# Combine label dfs

# relabel tips on trees
old_tip_labels_1 <- RNA1_id$accession
new_tip_labels_1 <- RNA1_id$RNA1

old_tip_labels_2 <- RNA2_id$accession
new_tip_labels_2 <- RNA2_id$RNA2

old_tip_labels_3 <- RNA3_id$accession
new_tip_labels_3 <- RNA3_id$RNA3

old_tip_labels_c <- Chaq_id$accession
new_tip_labels_c <- Chaq_id$Chaq

# relabel tips
rna1_tree_relabeled <- updateLabel(rna1_tree_rooted, old_tip_labels_1, new_tip_labels_1)
rna2_tree_relabeled <- updateLabel(rna2_tree_rooted, old_tip_labels_2, new_tip_labels_2)
rna3_tree_relabeled <- updateLabel(rna3_tree_rooted, old_tip_labels_3, new_tip_labels_3)
Chaq_tree_relabeled <- updateLabel(Chaq_tree_rooted, old_tip_labels_c, new_tip_labels_c)

# collapse 0 branch lengths
rna1_tree_relabeled <- di2multi(rna1_tree_relabeled, tol = 0.001)
rna2_tree_relabeled <- di2multi(rna2_tree_relabeled, tol = 0.001)
rna3_tree_relabeled <- di2multi(rna3_tree_relabeled, tol = 0.001)
Chaq_tree_relabeled <- di2multi(Chaq_tree_relabeled, tol = 0.001)

# RNA 1 and 2
# Make association for labels
assoc_1_2 <- full_join(RNA1_id, RNA2_id, by = "isolate") %>% 
  select(RNA1, RNA2)

assoc_1 <- assoc_1_2$RNA1
assoc_2 <- assoc_1_2$RNA2

assoc <- cbind(assoc_1, assoc_2)
  
# RNA 1 and 2
# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, rna2_tree_relabeled, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna1_rna2_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("navyblue", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()

# RNA 1 and 3
# Make association for labels
assoc_1_3 <- full_join(RNA1_id, RNA3_id, by = "isolate") %>% 
  select(RNA1, RNA3)

assoc_1 <- assoc_1_3$RNA1
assoc_3 <- assoc_1_3$RNA3

assoc <- cbind(assoc_1, assoc_3)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, rna3_tree_relabeled, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna1_rna3_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("navyblue", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()

# RNA2 RNA3
# Make association for labels
assoc_2_3 <- full_join(RNA2_id, RNA3_id, by = "isolate") %>% 
  select(RNA2, RNA3)

assoc_2 <- assoc_2_3$RNA2
assoc_3 <- assoc_2_3$RNA3

assoc <- cbind(assoc_2, assoc_3)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna2_tree_relabeled, rna3_tree_relabeled, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna2_rna3_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("navyblue", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()


# RNA1 Chaq
# Make association for labels
assoc_1_c <- full_join(RNA1_id, Chaq_id, by = "isolate") %>% 
  select(RNA1, Chaq)

assoc_1 <- assoc_1_c$RNA1
assoc_c <- assoc_1_c$Chaq

assoc <- cbind(assoc_1, assoc_c)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, Chaq_tree_relabeled, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna1_Chaq_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("navyblue", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()

# RNA2 Chaq
# Make association for labels
assoc_2_c <- full_join(RNA2_id, Chaq_id, by = "isolate") %>% 
  select(RNA2, Chaq)

assoc_2 <- assoc_2_c$RNA2
assoc_c <- assoc_2_c$Chaq

assoc <- cbind(assoc_2, assoc_c)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna2_tree_relabeled, Chaq_tree_relabeled, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna2_Chaq_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("navyblue", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()

# RNA3 Chaq
# Make association for labels
assoc_3_c <- full_join(RNA3_id, Chaq_id, by = "isolate") %>% 
  select(RNA3, Chaq)

assoc_3 <- assoc_3_c$RNA3
assoc_c <- assoc_3_c$Chaq

assoc <- cbind(assoc_3, assoc_c)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna3_tree_relabeled, Chaq_tree_relabeled, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna3_Chaq_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("navyblue", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()



