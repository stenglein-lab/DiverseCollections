library(phytools)
library(ape)
library(tidyverse)
library(readxl)
library(TreeTools)

## Remake trees with coinfection samples removed
# Read in new trees
RNA1 <- read.tree("analyses/trees/tanglegrams/RNA1_no_coinfection.fasta.contree")
RNA2 <- read.tree("analyses/trees/tanglegrams/RNA2_updated_no_coinf.fasta.contree")
RNA3 <- read.tree("analyses/trees/tanglegrams/RNA3_updated_no_coinfect.fasta.contree")
Chaq <- read.tree("analyses/trees/tanglegrams/Chaq_updated_no_coinf.fasta.contree")

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

# Manually fix a few names (EVE, SRA seqs, available seqs, remove coinfections)
#write.csv(RNA1_id, file = "analyses/trees/tanglegrams/RNA1_newnames.csv")
RNA1_id <- read_csv("analyses/trees/tanglegrams/RNA1_newnames.csv") %>% 
  rename("RNA1" = id)

#write.csv(RNA2_id, file = "analyses/trees/tanglegrams/RNA2_newnames.csv")
RNA2_id <- read_csv("analyses/trees/tanglegrams/RNA2_newnames.csv") %>% 
  rename("RNA2" = id)


#write.csv(RNA3_id, file = "analyses/trees/tanglegrams/RNA3_newnames.csv")
RNA3_id <- read_csv("analyses/trees/tanglegrams/RNA3_newnames.csv") %>% 
  rename("RNA3" = id)

#write.csv(Chaq_id, file = "analyses/trees/tanglegrams/Chaq_newnames.csv")
Chaq_id <- read_csv("analyses/trees/tanglegrams/Chaq_newnames.csv") %>% 
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
             link.col=make.transparent("grey", 0.25),
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
             link.col=make.transparent("grey", 0.25),
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
             link.col=make.transparent("grey", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()

# ----------------------------------------------------------------------------------------------
# *** Chaq Tanglegrams **
#
# There is an issue with chaq tanglegrams: because chaq is never present in clade B infections,
# there is a "missing" clade in galbut-chaq tanglegrams. 
#
# This causes tanglegrams to not look great: spacing is off
# 
# To fix this, add in a stand-in subtree to make chaq-containing tanglegrams align better. 
#
# This will be an artifical clade in the trees that we will delete manually 
# after tanglegram creation and PDF export
# ----------------------------------------------------------------------------------------------


# add a new branch of length 0.2 to root of chaq tree to start filler clade
chaq_tree_with_filler <- AddTip(Chaq_tree_relabeled, where = 0, edgeLength=0.2, label="filler_clade")

# add 18 new filler tips on 0-length branches to simulate "missing" clade B chaq sequences
for (new_tip in 1:17) {
  # add as polytomy
  chaq_tree_with_filler <- AddTip(chaq_tree_with_filler, 
                                  where = "filler_clade", 
                                  edgeLength=0,
                                  lengthBelow = 0,
                                  label=paste0("filler_tip_", new_tip))
}
# ggtree(chaq_tree_with_filler) + geom_tiplab(size=2)

# RNA1 Chaq
# Make association for labels
assoc_1_c <- full_join(RNA1_id, Chaq_id, by = "isolate") %>% 
  select(RNA1, Chaq)

assoc_1 <- assoc_1_c$RNA1
assoc_c <- assoc_1_c$Chaq

assoc <- cbind(assoc_1, assoc_c)


# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, chaq_tree_with_filler, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna1_Chaq_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("grey", 0.25),
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
cophy <- cophylo(rna2_tree_relabeled, chaq_tree_with_filler, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna2_Chaq_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("grey", 0.25),
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
cophy <- cophylo(rna3_tree_relabeled, chaq_tree_with_filler, assoc = assoc)

# will save as PDF
pdf(file = "analyses/trees/tanglegrams/rna3_Chaq_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("grey", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()

