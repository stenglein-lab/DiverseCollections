library(tidyverse)
library(seqinr) # for dist.mat

# read in MSAs 
RNA1_msa <- read.alignment("analyses/trees/alignments/RNA1_nucleotide_alignment.fasta", format="fasta")
RNA2_msa <- read.alignment("analyses/trees/alignments/RNA2_nucleotide_alignment.fasta", format="fasta")
RNA3_msa <- read.alignment("analyses/trees/alignments/RNA3_nucelotide_alignment.fasta", format="fasta")
Chaq_msa <- read.alignment("analyses/trees/alignments/Chaq_nucleotide_alignment.fasta", format="fasta")

# calculate pairwise distances using seqinr dist.alignment function
RNA1_dist <- tibble(segment = "RNA1", distance = as.numeric(dist.alignment(RNA1_msa, matrix="identity", gap=T)))
RNA2_dist <- tibble(segment = "RNA2", distance = as.numeric(dist.alignment(RNA2_msa, matrix="identity", gap=T)))
RNA3_dist <- tibble(segment = "RNA3", distance = as.numeric(dist.alignment(RNA3_msa, matrix="identity", gap=T)))
Chaq_dist <- tibble(segment = "Chaq", distance = as.numeric(dist.alignment(Chaq_msa, matrix="identity", gap=T)))

# combine 
distances <- rbind(RNA1_dist, RNA2_dist, RNA3_dist, Chaq_dist)

# the distance reported by dist.alignment is the square root of the pairwise distances.
# from the dist.alignment help text:
# "The resulting matrix contains the squared root of the pairwise distances. "
# "For example, if identity between 2 sequences is 80 the squared root of (1.0 - 0.8) i.e. 0.4472136"
pairwise_percent_identities <- distances %>% mutate(pct_identity = 100 * (1-(distance^2)))

# plot histogram of pairwise percent identities
ggplot(pairwise_percent_identities) +
  geom_histogram(aes(x=pct_identity), bins=101) +
  facet_wrap(~segment, ncol = 1) +
  theme_bw(base_size = 12) +
  xlab("Pairwise percent nucleotide identity") +
  ylab("Number pairwise comparisons")

# calculate averages and stats

average_percent_identities <- pairwise_percent_identities %>% 
  group_by(segment) %>% 
  summarize(median_pct_identity = median(pct_identity)) %>% 
  pivot_wider(names_from = segment, values_from = median_pct_identity)

# Kruskal-Wallis test
kw <- kruskal.test(pct_identity ~ segment, data = pairwise_percent_identities)
kw$p.value

# if calculated p-value is "0", that means
# it is <2.2e-16, the minimum floating point value (double)
if (kw$p.value == 0) {
  p_value = paste0("p < ", 
                   sprintf("%0.1e", .Machine$double.eps))
} else {
  p_value = paste0("p = ", 
                   sprintf("%0.1e", kw$p.value))
}

# generate output text for paper
output_text <- paste0(
  "Galbut virus segments exhibited varying levels of sequence diversity (Kruskal-Wallis test",
  p_value, 
  "). ",
  "RNA3 sequences exhibited the most diversity and RNA2 the least (Figs. 4, 5, 6, and 8). ",
  "The median percent identity between RNA1, RNA2, RNA3, and chaq sequences was ",
  sprintf("%0.1f", average_percent_identities %>% pull(RNA1) * 100),
  "%, ",
  sprintf("%0.1f", average_percent_identities %>% pull(RNA2) * 100),
  "%, ",
  sprintf("%0.1f", average_percent_identities %>% pull(RNA3) * 100),
  "%, ",
  sprintf("%0.1f", average_percent_identities %>% pull(Chaq) * 100),
  "%, respectively (Fig. 8). ",
  "Because there are no clade B chaq virus sequences, the diversity of chaq sequences is underestimated by this measure."
)

output_text
