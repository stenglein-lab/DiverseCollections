library(tidyverse)
library(readxl)

metadata <- read_xlsx("analyses/data/Metadata_Table.xlsx")

metadata <- metadata %>% 
  mutate(present = str_replace(present, "y", "yes"),
         present = str_replace(present, "n", "no"),
         complete = str_replace(complete, "y", "yes"),
         complete = str_replace(complete, "n", "no"),
         low_pos = str_replace(low_pos, "y", "yes"),
         low_pos = str_replace(low_pos, "n", "no"),
         expected_result = str_replace(expected_result, "y", "yes"),
         expected_result = str_replace(expected_result, "n", "no"),
         chaq_present = str_replace(chaq_present, "n", "no"),
         chaq_present = str_replace(chaq_present, "y", "yes"),
         chaq_complete = str_replace(chaq_complete, "n", "no"),
         chaq_complete = str_replace(chaq_complete, "y", "yes"),
         galbut_present = str_replace(galbut_present, "n", "no"),
         galbut_present = str_replace(galbut_present, "y", "yes"),
         galbut_complete = str_replace(galbut_complete, "n", "no"),
         galbut_complete = str_replace(galbut_complete, "y", "yes"))

location_v_positive <- metadata %>% 
  group_by(location, present) %>% 
  summarise(n())

location_v_complete <- metadata %>% 
  filter(present == "yes") %>% 
  group_by(location, complete) %>% 
  summarise(n())

location_v_chaq_positive <- metadata %>% 
  group_by(location, chaq_present) %>% 
  summarise(n())

location_v_chaq_complete <- metadata %>% 
  filter(chaq_present == "yes") %>% 
  group_by(location, chaq_complete) %>% 
  summarise(n())

location_v_galbut_positive <- metadata %>% 
  group_by(location, galbut_present) %>% 
  summarise(n())

location_v_galbut_complete <- metadata %>% 
  filter(galbut_present == "yes") %>% 
  group_by(location, galbut_complete) %>% 
  summarise(n())

# Normalize
metadata2 <- metadata %>% 
  mutate(ct_use = case_when(genotype == "A" & ct_1600_1601 < ct_2165_2170 ~ ct_1600_1601,
                            genotype == "A" & ct_1600_1601 > ct_2165_2170 ~ ct_2165_2170,
                            genotype == "B" ~ ct_2165_2170),
         delta_ct = ct_use - ct_rpl) %>% 
  select(sample_name, location, species, low_pos, ct_1600_1601, ct_2165_2170, 
         genotype, ct_rpl, chaq_present, ct_use, delta_ct) %>% 
# remove ohio samples that were definitely not positive for galbut virus  
  filter(!row_number() %in% c(119, 120, 122, 123, 124, 125, 126))

ct_genotype <- ggplot(filter(metadata2)) +
  geom_jitter(aes(x = genotype, y = delta_ct, color = chaq_present, 
                  shape = low_pos), size = 4, alpha = 0.7, width = 0.25, 
              height = 0) +
  scale_colour_manual(values = c("darkblue", "indianred4")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(~location) +
  labs(x = "Galbut virus genotype (based on RT-qPCR)",
       y = "Number of cycles to detection", 
       color = "Chaq virus present\n by sequencing", 
       shape = "RT-qPCR low positve")

ct_genotype
ggsave("analyses/plots/galbut_present.pdf", units = "in", width = 14, height = 8)


