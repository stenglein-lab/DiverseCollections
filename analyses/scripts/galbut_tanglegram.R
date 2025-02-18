library(phytools)
library(ape)

# a script to create a galbut virus RNA1/RNA2 tanglegram
# MDS 12/19/2024

# RNA1 RNA2
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

# collapse 0 branch lengths
rna1_tree_relabeled <- di2multi(rna1_tree_relabeled, tol = 0.001)
rna2_tree_relabeled <- di2multi(rna2_tree_relabeled, tol = 0.001)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, rna2_tree_relabeled)

# will save as PDF
pdf(file = "rna1_rna2_tanglegram.pdf", width=8.5, height=11)

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


# RNA1 RNA3
# use ape read.tree
RNA3 <- read.tree("analyses/trees/RNA3/RNA3_nucelotide_alignment.fasta.contree")

# midpoint root both trees
rna3_tree_rooted <- midpoint_root(RNA3)

# need to relabel tips so they have shared identifiers (here: strain/isolate name)
# https://rdrr.io/cran/ape/man/updateLabel.html

# read in map of sequence name (in tree) -> isolate identifier 
isolate_map1 <- read.delim("analyses/trees/tanglegrams/RNA1_RNA3_labels.txt", header=F, sep="\t")
colnames(isolate_map1) <- c("sequence_id", "isolate", "segment")

old_tip_labels <- isolate_map1$sequence_id
new_tip_labels <- isolate_map1$isolate


# relabel tips
rna1_tree_relabeled <- updateLabel(rna1_tree_rooted, old_tip_labels, new_tip_labels)
rna3_tree_relabeled <- updateLabel(rna3_tree_rooted, old_tip_labels, new_tip_labels)

# collapse 0 branch lengths
rna1_tree_relabeled <- di2multi(rna1_tree_relabeled, tol = 0.001)
rna3_tree_relabeled <- di2multi(rna3_tree_relabeled, tol = 0.001)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, rna3_tree_relabeled)

# will save as PDF
pdf(file = "rna1_rna3_tanglegram.pdf", width=8.5, height=11)

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
# need to relabel tips so they have shared identifiers (here: strain/isolate name)
# https://rdrr.io/cran/ape/man/updateLabel.html

# read in map of sequence name (in tree) -> isolate identifier 
isolate_map1 <- read.delim("analyses/trees/tanglegrams/RNA2_RNA3_labels.txt", header=F, sep="\t")
colnames(isolate_map1) <- c("sequence_id", "isolate", "segment")

old_tip_labels <- isolate_map1$sequence_id
new_tip_labels <- isolate_map1$isolate


# relabel tips
rna2_tree_relabeled <- updateLabel(rna2_tree_rooted, old_tip_labels, new_tip_labels)
rna3_tree_relabeled <- updateLabel(rna3_tree_rooted, old_tip_labels, new_tip_labels)

# collapse 0 branch lengths
rna1_tree_relabeled <- di2multi(rna1_tree_relabeled, tol = 0.001)
rna2_tree_relabeled <- di2multi(rna2_tree_relabeled, tol = 0.001)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna2_tree_relabeled, rna3_tree_relabeled)

# will save as PDF
pdf(file = "rna2_rna3_tanglegram.pdf", width=8.5, height=11)

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
# read in 2 individual trees
# use ape read.tree
Chaq <- read.tree("analyses/trees/Chaq/Chaq_alignment_recomb_removed.fasta.contree")

# midpoint root both trees
Chaq_tree_rooted <- midpoint_root(Chaq)

# need to relabel tips so they have shared identifiers (here: strain/isolate name)
# https://rdrr.io/cran/ape/man/updateLabel.html

# read in map of sequence name (in tree) -> isolate identifier 
isolate_map1 <- read.delim("analyses/trees/tanglegrams/RNA1_Chaq_labels.txt", header=F, sep="\t")
colnames(isolate_map1) <- c("sequence_id", "isolate", "segment")

old_tip_labels <- isolate_map1$sequence_id
new_tip_labels <- isolate_map1$isolate


# relabel tips
rna1_tree_relabeled <- updateLabel(rna1_tree_rooted, old_tip_labels, new_tip_labels)
Chaq_tree_relabeled <- updateLabel(Chaq_tree_rooted, old_tip_labels, new_tip_labels)

# collapse 0 branch lengths
rna1_tree_relabeled <- di2multi(rna1_tree_relabeled, tol = 0.001)
rna2_tree_relabeled <- di2multi(rna2_tree_relabeled, tol = 0.001)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna1_tree_relabeled, Chaq_tree_relabeled)

# will save as PDF
pdf(file = "rna1_Chaq_tanglegram.pdf", width=8.5, height=11)

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

# need to relabel tips so they have shared identifiers (here: strain/isolate name)
# https://rdrr.io/cran/ape/man/updateLabel.html

# read in map of sequence name (in tree) -> isolate identifier 
isolate_map1 <- read.delim("analyses/trees/tanglegrams/RNA2_Chaq_labels.txt", header=F, sep="\t")
colnames(isolate_map1) <- c("sequence_id", "isolate", "segment")

old_tip_labels <- isolate_map1$sequence_id
new_tip_labels <- isolate_map1$isolate


# relabel tips
rna2_tree_relabeled <- updateLabel(rna2_tree_rooted, old_tip_labels, new_tip_labels)
Chaq_tree_relabeled <- updateLabel(Chaq_tree_rooted, old_tip_labels, new_tip_labels)

# collapse 0 branch lengths
rna1_tree_relabeled <- di2multi(rna1_tree_relabeled, tol = 0.001)
rna2_tree_relabeled <- di2multi(rna2_tree_relabeled, tol = 0.001)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna2_tree_relabeled, Chaq_tree_relabeled)

# will save as PDF
pdf(file = "rna2_chaq_tanglegram.pdf", width=8.5, height=11)

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
# need to relabel tips so they have shared identifiers (here: strain/isolate name)
# https://rdrr.io/cran/ape/man/updateLabel.html

# read in map of sequence name (in tree) -> isolate identifier 
isolate_map1 <- read.delim("analyses/trees/tanglegrams/RNA3_Chaq_labels.txt", header=F, sep="\t")
colnames(isolate_map1) <- c("sequence_id", "isolate", "segment")

old_tip_labels <- isolate_map1$sequence_id
new_tip_labels <- isolate_map1$isolate


# relabel tips
rna3_tree_relabeled <- updateLabel(rna3_tree_rooted, old_tip_labels, new_tip_labels)
Chaq_tree_relabeled <- updateLabel(Chaq_tree_rooted, old_tip_labels, new_tip_labels)

# collapse 0 branch lengths
rna1_tree_relabeled <- di2multi(rna1_tree_relabeled, tol = 0.001)
rna2_tree_relabeled <- di2multi(rna2_tree_relabeled, tol = 0.001)

# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(rna3_tree_relabeled, Chaq_tree_relabeled)

# will save as PDF
pdf(file = "rna3_chaq_tanglegram.pdf", width=8.5, height=11)

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




