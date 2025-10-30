# This script was used to plot pair-wise distances to visualize which segments were more diverse than others.
library(tidyverse)
library(seqinr) # for dist.mat
library(MKinfer)

# read in MSAs 
RNA1_msa <- read.alignment("analyses/trees/alignments/RNA1_nucleotide_alignment.fasta", format="fasta")
RNA2_msa <- read.alignment("analyses/trees/alignments/RNA2_nucleotide_alignment.fasta", format="fasta")
RNA3_msa <- read.alignment("analyses/trees/alignments/RNA3_nucleotide_alignment.fasta", format="fasta")
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

pairwise_percent_identities$segment <- fct_relevel(pairwise_percent_identities$segment, "Chaq", after=3)

# plot histogram of pairwise percent identities
ggplot(pairwise_percent_identities) +

    geom_histogram(aes(x=pct_identity), bins=101) +
  facet_wrap(~segment, ncol = 1, scales="free_y") +
  theme_bw(base_size = 12) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  xlab("Pairwise nucleotide identity (%)") +
  ylab("Number pairwise alignments")

ggsave("analyses/plots/Supplemental_figure_8_pairwise_pct_identity_by_segment.pdf", width=10, height=7.5, units="in")


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


# do permutation/bootstrap testing for significance of differences in pairwise differences
# this is a homebrew bootstrap testing function
assess_diversity_difference <- function(seg_1, seg_2) {
  
  debug <- 0
  if (debug) {
    seg_1 = "RNA1"
    seg_2 = "RNA1"
  }
  
  seg_1_ids <- filter(pairwise_percent_identities, segment == seg_1) %>% pull(pct_identity)
  seg_2_ids <- filter(pairwise_percent_identities, segment == seg_2) %>% pull(pct_identity)
  
  num_seg_1 <- length(seg_1_ids) 
  num_seg_2 <- length(seg_2_ids)
  
  combined_ids <- c(seg_1_ids, seg_2_ids)
  
  observed_diff <- abs(median(seg_1_ids) - median(seg_2_ids))
  
  number_permutations <- 5000
  number_permuted_exceeds <- 0
  
  # do bootstrapping
  for (i in 1:number_permutations) {
     seg_1_perm_ids <- sample(combined_ids, num_seg_1, replace = T)
     seg_2_perm_ids <- sample(combined_ids, num_seg_2, replace = T)
     permuted_diff <- abs(median(seg_1_perm_ids) - median(seg_2_perm_ids))
     
     if (permuted_diff > observed_diff) {
       number_permuted_exceeds <- number_permuted_exceeds + 1
     }
  }
  
  # return p-value from permutation testing
  if (number_permuted_exceeds > 0) {
    number_permuted_exceeds / number_permutations
  } else {
    paste0("< 1 / ", number_permutations)
  }
}

# do similar bootstrap testing using the MKinfer package
assess_diversity_difference_mkinfer <- function(seg_1, seg_2) {
  
  debug <- 0
  if (debug) {
    seg_1 = "RNA1"
    seg_2 = "RNA2"
  }
  
  seg_1_ids <- filter(pairwise_percent_identities, segment == seg_1) %>% pull(pct_identity)
  seg_2_ids <- filter(pairwise_percent_identities, segment == seg_2) %>% pull(pct_identity)
  
  btt <- boot.t.test(seg_1_ids, seg_2_ids, alternative = "two.sided", 
                     paired = F, R = 5000)
  
  # btt$p.value
  btt$boot.p.value
}

# homebrew bootstrap testing 
# a control: the same distribution should not be significantly different
assess_diversity_difference("RNA1", "RNA1")

# test for differences between segments
assess_diversity_difference("RNA1", "RNA2")
assess_diversity_difference("RNA1", "RNA3")
assess_diversity_difference("RNA1", "Chaq")
assess_diversity_difference("RNA2", "RNA3")
assess_diversity_difference("RNA2", "Chaq")
assess_diversity_difference("RNA3", "Chaq")

# bootstrap testing using MKinfer boot.t.test to double check
# a control: the same distribution should not be significantly different
assess_diversity_difference_mkinfer("RNA1", "RNA1")

# test for differences between segments
assess_diversity_difference_mkinfer("RNA1", "RNA2")
assess_diversity_difference_mkinfer("RNA1", "RNA3")
assess_diversity_difference_mkinfer("RNA1", "Chaq")
assess_diversity_difference_mkinfer("RNA2", "RNA3")
assess_diversity_difference_mkinfer("RNA2", "Chaq")
assess_diversity_difference_mkinfer("RNA3", "Chaq")


# generate output text for paper
output_text <- paste0(
  "Galbut virus segments exhibited varying levels of sequence diversity (Kruskal-Wallis test",
  p_value, 
  "). ",
  "RNA3 sequences exhibited the most diversity and RNA2 the least (Figs. 4, 5, 6, and 8). ",
  "The median percent identity between RNA1, RNA2, RNA3, and chaq sequences was ",
  sprintf("%0.1f", average_percent_identities %>% pull(RNA1)),
  "%, ",
  sprintf("%0.1f", average_percent_identities %>% pull(RNA2)),
  "%, ",
  sprintf("%0.1f", average_percent_identities %>% pull(RNA3)),
  "%, ",
  sprintf("%0.1f", average_percent_identities %>% pull(Chaq)),
  "%, respectively (Fig. 8). ",
  "Because there are no clade B chaq virus sequences, the diversity of chaq sequences is underestimated by this measure."
)

output_text



