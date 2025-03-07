library(ape)
library(adephylo)
library(phytools)
library(tidyverse)
library(ggridges)
library(ggtree)
library(readxl)

# Read in trees
RNA1 <- read.tree("analyses/trees/RNA1/RNA1_nucleotide_alignment.fasta.contree")
RNA2 <- read.tree("analyses/trees/RNA2/RNA2_nucleotide_alignment.fasta.contree")
RNA3 <- read.tree("analyses/trees/RNA3/RNA3_nucelotide_alignment.fasta.contree")
Chaq <- read.tree("analyses/trees/Chaq/Chaq_alignment_recomb_removed.fasta.contree")

# Midpoint root trees
RNA1_tree_rooted <- midpoint.root(RNA1)
RNA2_tree_rooted <- midpoint.root(RNA2)
RNA3_tree_rooted <- midpoint.root(RNA3)
Chaq_tree_rooted <- midpoint.root(Chaq)

# Get root to tip distances
RNA1_root_dists <- as.data.frame(distRoot(RNA1_tree_rooted))
RNA2_root_dists <- as.data.frame(distRoot(RNA2_tree_rooted))
RNA3_root_dists <- as.data.frame(distRoot(RNA3_tree_rooted))
Chaq_root_dists <- as.data.frame(distRoot(Chaq_tree_rooted))

# Rename column
colnames(RNA1_root_dists) <- c("dist")
colnames(RNA2_root_dists) <- c("dist")
colnames(RNA3_root_dists) <- c("dist")
colnames(Chaq_root_dists) <- c("dist")

# Add segment column
RNA1_root_dists$segment="RNA1"
RNA2_root_dists$segment="RNA2"
RNA3_root_dists$segment="RNA3"
Chaq_root_dists$segment="Chaq"

# Combine dataframes
root_dists <- rbind(RNA1_root_dists, RNA2_root_dists, RNA3_root_dists, Chaq_root_dists)
root_dists <- root_dists %>% 
  mutate(segment = factor(segment, levels = c("RNA1", "RNA2", "RNA3", "Chaq")))

# Plot
plot <- ggplot(root_dists) +
  geom_boxplot(aes(y = dist, x = segment), color="black", size=0.25, alpha=0.8) +
  geom_point(aes(x = segment, y = dist, color = segment), alpha = 0.75, size = 2) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(y = "Tip to Root Distances", x = "", fill = "Segment")

plot
ggsave("analyses/plots/root_dists.pdf", units = "in", width = 14, height = 8)

# Read in clade data
RNA1_clade <- read_xlsx("analyses/trees/RNA1/RNA1_clade_info.xlsx") %>% 
  mutate(accession = str_replace(accession, "EVE", "galbut_EVE_from_Wallace_supp_fig_S10"))
RNA2_clade <- read_xlsx("analyses/trees/RNA2/RNA2_clade_info.xlsx")
RNA3_clade <- read_xlsx("analyses/trees/RNA3/RNA3_clade_info.xlsx")
Chaq_clade <- read_xlsx("analyses/trees/Chaq/Chaq_clade_info.xlsx")

# Combine dfs 
all_clade <- rbind(RNA1_clade, RNA2_clade, RNA3_clade, Chaq_clade) %>% 
  select(accession, clade)

# Add first column name in root_dist
root_dists <- rownames_to_column(root_dists, var = "accession")

# Join root_dists and clade data
all_data <- full_join(all_clade, root_dists, by = "accession")

# Plot
plot2 <- ggplot(all_data) +
  geom_boxplot(aes(y = dist, x = segment), color="black", size=0.25, alpha=0.8) +
  geom_point(aes(x = segment, y = dist, color = segment), alpha = 0.75, size = 2) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~clade, nrow = 3) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(y = "Tip to Root Distances", x = "", fill = "Segment")

plot2
ggsave("analyses/plots/root_dists_clade.pdf", units = "in", width = 14, height = 8)
