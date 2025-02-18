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
