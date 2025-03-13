library(tidyverse)

chaq <- read.csv("analyses/trees/Chaq/Chaq_FEL.csv") %>% 
  mutate(segment = "Chaq")

RNA1 <- read.csv("analyses/trees/RNA1/RNA1_FEL.csv") %>% 
  mutate(segment = "RNA1")

RNA2 <- read.csv("analyses/trees/RNA2/RNA2_FEL.csv") %>% 
  mutate(segment = "RNA2")

RNA3_sim <- read.csv("analyses/trees/RNA3/RNA3_FEL_sim.csv") %>% 
  mutate(segment = "RNA3")

RNA3_else <- read.csv("analyses/trees/RNA3/RNA3_FEL.csv") %>% 
  mutate(segment = "RNA3")

all_invariable <- rbind(chaq, RNA1, RNA2, RNA3_else, RNA3_sim) %>% 
  filter(class == "Invariable") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_invariable, file = "analyses/trees/invariable_selection.csv")

all_neutral <- rbind(chaq, RNA1, RNA2, RNA3_else, RNA3_sim) %>% 
  filter(class == "Neutral") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_neutral, file = "analyses/trees/invariable_neutral.csv")

all_purifying <- rbind(chaq, RNA1, RNA2, RNA3_else, RNA3_sim) %>% 
  filter(class == "Purifying") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_purifying, file = "analyses/trees/purifying_selection.csv")

all_diversifying <- rbind(chaq, RNA1, RNA2, RNA3_else, RNA3_sim) %>% 
  filter(class == "Diversifying") %>% 
  select(segment, codon, alpha, beta, alpha.beta, LRT, p.value, Total.branch.length)

write.csv(all_diversifying, file = "analyses/trees/diversifying_selection.csv")

significant_diversifying <- filter(all_diversifying, p.value < 0.05)
significant_purifying    <- filter(all_purifying,    p.value < 0.05)

significant_diversifying$selection_type <- "diversifying"
significant_purifying$selection_type    <- "purifying"

all_significant <- rbind(significant_diversifying, significant_purifying)

# how long are coding sequence?
# these #s based on FoCo17 sequences, MT742160 - MT742163

# *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
# NOTE!  The codons in the table above extend beyond these apparent sizes.  Need to check!
# *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!

RNA1_codons <- 541
RNA2_codons <- 495
RNA3_codons <- 443
Chaq_codons <- 320

codon_lengths <- tibble(
  segment    = c("RNA1",     "RNA2",       "RNA3",      "Chaq"),
  num_codons = c(RNA1_codons, RNA2_codons, RNA3_codons, Chaq_codons)
)

# reorder for plotting: chaq last
all_significant$segment <- fct_relevel(all_significant$segment, "Chaq", after=3)
codon_lengths$segment <- fct_relevel(codon_lengths$segment, "Chaq", after=3)

ggplot(all_significant) +
  geom_rect(data=codon_lengths, 
            aes(xmin = 0, xmax = num_codons, ymin=0, ymax=1),
            linewidth=0.25, color="black", fill="grey95") +
  geom_segment(data=all_significant, 
               aes(x=codon, xend=codon, y=0, yend=1, color=selection_type),
               linewidth=1) +
  facet_wrap(~segment, ncol = 1) + 
  theme_bw(base_size = 12) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        strip.text = element_text(hjust=0.03),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()) +
  scale_color_manual(values=c("coral3", "slateblue")) + 
  xlab("Codon") +
  ylab("")

ggsave("analyses/plots/Fig_X_sites_under_selection.pdf", width=7.5, height=6, units="in")
