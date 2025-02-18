library(phytools)
library(ape)

# a script to create a galbut virus RNA1/RNA2 tanglegram
# MDS 12/19/2024

# read in 2 individual trees
# use ape read.tree
RNA1 <- read.tree("analyses/trees/RNA1/RNA1_nucleotide_alignment.fasta.contree")
RNA2 <- read.tree("analyses/trees/RNA2/RNA2_nucleotide_alignment.fasta.contree")

# midpoint root both trees
rna1_tree_rooted <- midpoint_root(RNA1)
rna2_tree_rooted <- midpoint_root(RNA2)

# need to relabel tips so they have shared identifiers (here: strain/isolate name)
# https://rdrr.io/cran/ape/man/updateLabel.html

# read in map of sequence name (in tree) -> isolate identifier 
isolate_map1 <- read.delim("analyses/trees/tanglegrams/RNA1_RNA2_labels.txt", header=F, sep="\t")
colnames(isolate_map1) <- c("sequence_id", "isolate", "segment")

old_tip_labels <- isolate_map1$sequence_id
new_tip_labels <- isolate_map1$isolate


# relabel tips
rna1_tree_relabeled <- updateLabel(rna1_tree_rooted, old_tip_labels, new_tip_labels)
rna2_tree_relabeled <- updateLabel(rna2_tree_rooted, old_tip_labels, new_tip_labels)


# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, rna2_tree_relabeled)

# will save as PDF
pdf(file = "rna1_rna2_tanglegram.pdf", width=8.5, height=11)

# plot tanglegram
plot.cophylo(cophy, 
             link.type = "curved", 
             link.lwd=1,
             link.lty="solid",
             link.col=make.transparent("slateblue", 0.25),
             fsize=0.3, 
             pts=F,
             scale.bar = c(0.025, 0.025))

# this doesn't plot node support values, but it could.
# see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html

# turn off PDF
dev.off()



