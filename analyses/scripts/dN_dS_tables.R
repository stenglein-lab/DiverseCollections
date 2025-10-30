library(tidyverse)

chaq <- read.csv("analyses/data/Chaq_FEL.csv") %>% 
  mutate(segment = "Chaq")

RNA1 <- read.csv("analyses/data/RNA1_FEL.csv") %>% 
  mutate(segment = "RNA1")

RNA2 <- read.csv("analyses/data/RNA2_FEL.csv") %>% 
  mutate(segment = "RNA2")

RNA3 <- read.csv("analyses/data/RNA3_FEL.csv") %>% 
  mutate(segment = "RNA3")

all_invariable <- rbind(chaq, RNA1, RNA2, RNA3) %>% 
  filter(class == "Invariable") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_invariable, file = "analyses/tables/invariable_selection.csv")

all_neutral <- rbind(chaq, RNA1, RNA2, RNA3) %>% 
  filter(class == "Neutral") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_neutral, file = "analyses/tables/invariable_neutral.csv")

all_purifying <- rbind(chaq, RNA1, RNA2, RNA3) %>% 
  filter(class == "Purifying") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_purifying, file = "analyses/tables/Supp_table_6_purifying_selection.csv")

all_diversifying <- rbind(chaq, RNA1, RNA2, RNA3) %>% 
  filter(class == "Diversifying") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_diversifying, file = "analyses/tables/Supp_table_7_diversifying_selection.csv")

significant_diversifying <- filter(all_diversifying, p.value < 0.1)
significant_purifying    <- filter(all_purifying,    p.value < 0.1)

significant_diversifying$selection_type <- "diversifying"
significant_purifying$selection_type    <- "purifying"

all_significant <- rbind(significant_diversifying, significant_purifying)

# how long are coding sequence?
# these #s based on FoCo17 sequences, MT742160 - MT742163

RNA1_codons <- 541
RNA2_codons <- 495
RNA3_codons <- 449
Chaq_codons <- 320

codon_lengths <- tibble(
  segment    = c("RNA1",     "RNA2",       "RNA3",      "Chaq"),
  num_codons = c(RNA1_codons, RNA2_codons, RNA3_codons, Chaq_codons)
)

# reorder for plotting: chaq last
all_significant$segment <- fct_relevel(all_significant$segment, "Chaq", after=3)
codon_lengths$segment <- fct_relevel(codon_lengths$segment, "Chaq", after=3)

significant_counts <- 
  all_significant %>% 
  group_by(segment, selection_type) %>% 
  summarize(n=n(), .groups="drop")

# pull out specific values for paper text
rna1_pur <- significant_counts %>% filter(segment == "RNA1" & selection_type == "purifying") %>% pull(n)
rna2_pur <- significant_counts %>% filter(segment == "RNA2" & selection_type == "purifying") %>% pull(n)
rna3_pur <- significant_counts %>% filter(segment == "RNA3" & selection_type == "purifying") %>% pull(n)
chaq_pur <- significant_counts %>% filter(segment == "Chaq" & selection_type == "purifying") %>% pull(n)

rna1_div <- significant_counts %>% filter(segment == "RNA1" & selection_type == "diversifying") %>% pull(n)
rna2_div <- significant_counts %>% filter(segment == "RNA2" & selection_type == "diversifying") %>% pull(n)
rna3_div <- significant_counts %>% filter(segment == "RNA3" & selection_type == "diversifying") %>% pull(n)
chaq_div <- significant_counts %>% filter(segment == "Chaq" & selection_type == "diversifying") %>% pull(n)

# output text for paper
paste0(
  "Consistent with this, we found evidence of purifying selection in all segments, with ",
  sprintf("%0.0f, ", rna1_pur),
  sprintf("%0.0f, ", rna2_pur),
  sprintf("%0.0f, ", rna3_pur),
  sprintf("%0.0f ", chaq_pur),
  "sites experiencing purifying selection in galbut virus RNA 1, 2, 3, and chaq virus respectively ")
  
paste0(
  "There was little evidence of diversifying selection in galbut virus RNA 1, RNA 2, or chaq virus with ",
  sprintf("%0.0f, ", rna1_div),
  sprintf("%0.0f, ", rna2_div),
  sprintf("%0.0f, ", chaq_div),
  "diversifying sites identified, respectively ")

paste0(
  "In contrast, galbut virus RNA 3 showed more evidence of diversifying selection, with ",
  sprintf("%0.0f, ", rna3_div),
  " sites identified " )

ggplot(all_significant) +
  geom_rect(data=codon_lengths, 
            aes(xmin = 0, xmax = num_codons, ymin=0, ymax=1),
            linewidth=0.25, color="black", fill="grey95") +
  geom_segment(data=filter(all_significant, selection_type == "purifying"), 
               aes(x=codon, xend=codon, y=0, yend=1, color=selection_type),
               linewidth=0.25, alpha=0.75) +
  geom_segment(data=filter(all_significant, selection_type == "diversifying"), 
               aes(x=codon, xend=codon, y=0, yend=1, color=selection_type),
               linewidth=0.6) +
  facet_wrap(~segment, ncol = 1) + 
  theme_minimal(base_size = 12) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        # text = element_text(size = 20),
        axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  scale_color_manual(values=c("coral3", "slateblue")) + 
  labs(x = "Codon", y = "", color = "Selection Type") 

ggsave("analyses/plots/Figure_5_sites_under_selection.pdf", width=5, height=4, units="in")

