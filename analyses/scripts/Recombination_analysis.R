library(tidyverse)

# Read in the data
RNA3_maj <- read_csv("analyses/trees/recombination/RNA3_recomb_major-2.csv")
RNA3_min <- read_csv("analyses/trees/recombination/RNA3_recomb_minor-2.csv")

Chaq_maj <- read_csv("analyses/trees/recombination/Chaq_recomb_major.csv")
Chaq_min <- read_csv("analyses/trees/recombination/Chaq_recomb_minor.csv")

# Make plots
RNA3_plot <- ggplot(NULL, aes(x = `Window Start`, y = `Percent Identity`)) +
  geom_line(data = RNA3_maj, color = 'blue') +
  geom_line(data = RNA3_min, color = "orange") +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

RNA3_plot
ggsave("analyses/plots/RNA3_recombination.pdf", units = "in", width = 10, height = 8)  

Chaq_plot <- ggplot(NULL, aes(x = `Window Start`, y = `Percent Identity`)) +
  geom_line(data = Chaq_maj, color = 'blue') +
  geom_line(data = Chaq_min, color = "orange") +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

Chaq_plot
ggsave("analyses/plots/Chaq_recombination.pdf", units = "in", width = 10, height = 8)  
