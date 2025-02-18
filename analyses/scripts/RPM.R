library(tidyverse)
library(readxl)

all_counts <- read_xlsx("analyses/RPM/all_read_counts.xlsx") %>% 
  filter(count_type == "post_trimming") %>% 
  rename(id = sample_id) %>% 
  select(id, count)

remapping1 <- read_xlsx("analyses/RPM/remapping_num_aligned.xlsx")
remapping3 <- read_xlsx("analyses/RPM/remapping3_num_aligned.xlsx")
remappingM <- read_xlsx("analyses/RPM/remapping_mark.xlsx")

#parse to get sample ids
remapping1 <- remapping1 %>% 
  mutate(sample_id = str_replace(sample_id, "sim-20-TD-4.filt", "D_mel_USA_2020_sim-20-TD-4_Galbut_virus_RNA3_Complete_CDS_2")) %>% 
  separate(sample_id, into = c("D", "mel_sim", "loc", "year", "id", "galbut", "virus", "segment"), sep = "_") %>% 
  mutate(segment = str_replace(segment, "Complete", "Chaq"),
         segment = str_replace(segment, "Capsid", "RNA2"),
         segment = str_replace(segment, "RNA2s", "RNA2")) %>% 
  select(id, segment, num_reads) 

remapping3 <- remapping3 %>% 
  mutate(sample_id = str_replace(sample_id, "Penn-M-2_Chaq_virus_Complete_CDS_1", 
                                 "D_mel_USA_2023_Penn-M-2_Chaq_virus_Complete_CDS_1")) %>% 
  mutate(sample_id = str_replace(sample_id, "Penn-M-2_Chaq_virus_Complete_CDS_2", 
                                 "D_mel_USA_2023_Penn-M-2_Chaq_virus_Complete_CDS_2")) %>% 
  separate(sample_id, into = c("D", "mel_sim", "loc", "year", "id", "galbut", "virus", "segment"), sep = "_") %>% 
  mutate(segment = str_replace(segment, "Complete", "Chaq"),
         segment = str_replace(segment, "Capsid", "RNA2"),
         segment = str_replace(segment, "Partial", "Chaq")) %>% 
  select(id, segment, num_reads)

remappingM <- remappingM %>% 
  separate(id, into = c("D", "mel_sim", "loc", "year", "id", "galbut", "virus", "segment"), sep = "_") %>% 
  mutate(segment = str_replace(segment, "Complete", "Chaq")) %>% 
  mutate(num_reads = num_seqs) %>% 
  select(id, segment, num_reads)

# Get segments into own df
remap1_RNA1 <- remapping1 %>% 
  filter(segment == "RNA1")

remap1_RNA2 <- remapping1 %>% 
  filter(segment == "RNA2")

remap1_RNA3 <- remapping1 %>% 
  filter(segment == "RNA3")

remap1_Chaq <- remapping1 %>% 
  filter(segment == "Chaq")

remap3_RNA1 <- remapping3 %>% 
  filter(segment == "RNA1")

remap3_RNA2 <- remapping3 %>% 
  filter(segment == "RNA2")

remap3_RNA3 <- remapping3 %>% 
  filter(segment == "RNA3")

remap3_Chaq <- remapping3 %>% 
  filter(segment == "Chaq")

remapM_RNA2 <- remappingM %>% 
  filter(segment == "RNA2") 

remapM_RNA3 <- remappingM %>% 
  filter(segment == "RNA3")

remapM_Chaq <- remappingM %>% 
  filter(segment == "Chaq")

#join the df's - keeping all the remapping 3 and all other mapping in remapping 1
RNA1_remap1 <- anti_join(remap1_RNA1, remap3_RNA1, by = "id")
all_RNA1 <- full_join(remap3_RNA1, RNA1_remap1, by = c("segment", "id", "num_reads"))
all_RNA1 <- all_RNA1[-c(35), ]

RNA2_remap1 <- anti_join(remap1_RNA2, remap3_RNA2, by = "id")                    
all_RNA2 <- full_join(remap3_RNA2, RNA2_remap1, by = c("segment", "id", "num_reads"))
all_RNA2 <- full_join(all_RNA2, remapM_RNA2, by = c("segment", "id", "num_reads"))
all_RNA2 <- all_RNA2[-c(14, 31, 108, 107), ]

RNA3_remap1 <- anti_join(remap1_RNA3, remap3_RNA3, by = "id")                    
all_RNA3 <- full_join(remap3_RNA3, RNA3_remap1, by = c("segment", "id", "num_reads"))
all_RNA3 <- full_join(all_RNA3, remapM_RNA3, by = c("segment", "id", "num_reads"))
all_RNA3 <- all_RNA3[-c(30, 40, 57), ]

Chaq_remap1 <- anti_join(remap1_Chaq, remap3_Chaq, by = "id")                    
all_Chaq <- full_join(remap3_Chaq, Chaq_remap1, by = c("segment", "id", "num_reads"))
all_Chaq <- full_join(all_Chaq, remapM_Chaq, by = c("segment", "id", "num_reads"))
all_Chaq <- all_Chaq[-c(16, 23, 37), ]

# Read in metadata table
metadata <- read_csv("analyses/data/Metadata_Table.csv") %>% 
  mutate(id = sample_name) %>% 
  select(id, genotype, chaq_present) 

#combine all the data
all_RNA1 <- left_join(all_RNA1, metadata, by = "id")
all_RNA1 <- left_join(all_RNA1, all_counts, by = "id")

all_RNA2 <- left_join(all_RNA2, metadata, by = "id")
all_RNA2 <- left_join(all_RNA2, all_counts, by = "id")

all_RNA3 <- left_join(all_RNA3, metadata, by = "id")
all_RNA3 <- left_join(all_RNA3, all_counts, by = "id")

all_Chaq <- left_join(all_Chaq, metadata, by = "id")
all_Chaq <- left_join(all_Chaq, all_counts, by = "id")

all_data <- rbind(all_RNA1, all_RNA2, all_RNA3, all_Chaq) %>% 
  mutate(genotype = str_replace(genotype, "A", "Genotype A"),
         genotype = str_replace(genotype, "B", "Genotype B"),
         chaq_present = str_replace(chaq_present, "y", "Yes"),
         chaq_present = str_replace(chaq_present, "n", "No"))

# Calculate RPM
all_data <- all_data %>% 
  mutate(rpm = (num_reads/count)*1000000)

# Make plot
plot <- ggplot(all_data, aes(x = segment, y = rpm)) +
  geom_jitter(aes(colour = chaq_present, shape = chaq_present), alpha = 0.75, size = 3, width = 0.25, height = 0) +
  scale_colour_manual(values = c("#a6611a", "#01665e")) +
  scale_shape_manual(values = c(15, 16)) +
  scale_y_log10() +
  facet_wrap(~genotype) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "", y = "Reads Per Million", shape = "Chaq virus \npresent", 
       colour = "Chaq virus \npresent")

plot 
ggsave("analyses/plots/RPM.pdf", units = "in", width = 10, height = 8)  
