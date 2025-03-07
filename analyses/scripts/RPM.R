library(tidyverse)
library(readxl)
library(ggpmisc)

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

# Combine rows with two RNA1s
remap1_RNA1[75, 3] <- remap1_RNA1[75, 3] + remap1_RNA1[76, 3]
remap1_RNA1 <- remap1_RNA1[-c(76), ]

remap1_RNA1[88, 3] <- remap1_RNA1[88, 3] + remap1_RNA1[89, 3]
remap1_RNA1 <- remap1_RNA1[-c(89), ]

remap3_RNA1[44, 3] <- remap3_RNA1[44, 3] + remap3_RNA1[45, 3]
remap3_RNA1 <- remap3_RNA1[-c(45), ]

remap3_RNA1[54, 3] <- remap3_RNA1[54, 3] + remap3_RNA1[55, 3]
remap3_RNA1 <- remap3_RNA1[-c(55), ]


#Combine rows with two RNA2s
remap1_RNA2[78, 3] <- remap1_RNA2[78, 3] + remap1_RNA2[79, 3]
remap1_RNA2 <- remap1_RNA2[-c(79), ]

remap1_RNA2[92, 3] <- remap1_RNA2[92, 3] + remap1_RNA2[93, 3]
remap1_RNA2 <- remap1_RNA2[-c(93), ]

remap3_RNA2[46, 3] <- remap3_RNA2[46, 3] + remap3_RNA2[47, 3]
remap3_RNA2 <- remap3_RNA2[-c(47), ]

remap3_RNA2[55, 3] <- remap3_RNA2[55, 3] + remap1_RNA2[56, 3]
remap3_RNA2 <- remap3_RNA2[-c(56), ]

remap3_RNA2[56, 3] <- remap3_RNA2[56, 3] + remap1_RNA2[57, 3]
remap3_RNA2 <- remap3_RNA2[-c(57), ]

remap3_RNA2[31, 3] <- remap3_RNA2[31, 3] + remap1_RNA2[41, 3]
remap3_RNA2 <- remap3_RNA2[-c(41), ]

#Combine rows with two RNA3s
remap1_RNA3[1, 3] <- remap1_RNA3[1, 3] + remap1_RNA3[1, 3]
remap1_RNA3 <- remap1_RNA3[-c(2), ]

remap1_RNA3[62, 3] <- remap1_RNA3[62, 3] + remap1_RNA3[63, 3]
remap1_RNA3 <- remap1_RNA3[-c(63), ]

remap1_RNA3[76, 3] <- remap1_RNA3[76, 3] + remap1_RNA3[77, 3]
remap1_RNA3 <- remap1_RNA3[-c(77), ]

remap1_RNA3[83, 3] <- remap1_RNA3[83, 3] + remap1_RNA3[84, 3]
remap1_RNA3 <- remap1_RNA3[-c(84), ]

remap1_RNA3[90, 3] <- remap1_RNA3[90, 3] + remap1_RNA3[91, 3]
remap1_RNA3 <- remap1_RNA3[-c(91), ]

remap1_RNA3[96, 3] <- remap1_RNA3[96, 3] + remap1_RNA3[97, 3]
remap1_RNA3 <- remap1_RNA3[-c(97), ]

remap1_RNA3[98, 3] <- remap1_RNA3[98, 3] + remap1_RNA3[99, 3]
remap1_RNA3 <- remap1_RNA3[-c(99), ]

remap3_RNA3[1, 3] <- remap3_RNA3[1, 3] + remap3_RNA3[2, 3]
remap3_RNA3 <- remap3_RNA3[-c(2), ]

remap3_RNA3[7, 3] <- remap3_RNA3[7, 3] + remap3_RNA3[8, 3]
remap3_RNA3 <- remap3_RNA3[-c(8), ]

remap3_RNA3[11, 3] <- remap3_RNA3[11, 3] + remap3_RNA3[12, 3]
remap3_RNA3 <- remap3_RNA3[-c(12), ]

remap3_RNA3[16, 3] <- remap3_RNA3[16, 3] + remap3_RNA3[17, 3]
remap3_RNA3 <- remap3_RNA3[-c(17), ]

remap3_RNA3[42, 3] <- remap3_RNA3[42, 3] + remap3_RNA3[43, 3]
remap3_RNA3 <- remap3_RNA3[-c(43), ]

remap3_RNA3[46, 3] <- remap3_RNA3[46, 3] + remap3_RNA3[47, 3]
remap3_RNA3 <- remap3_RNA3[-c(47), ]

remap3_RNA3[50, 3] <- remap3_RNA3[50, 3] + remap3_RNA3[51, 3]
remap3_RNA3 <- remap3_RNA3[-c(51), ]

remap3_RNA3[37, 3] <- remap3_RNA3[37, 3] + remap3_RNA3[50, 3]
remap3_RNA3 <- remap3_RNA3[-c(50), ]

#Combine rows with two Chaqs

remap1_Chaq[39, 3] <- remap1_Chaq[39, 3] + remap1_Chaq[40, 3]
remap1_Chaq <- remap1_Chaq[-c(40), ]

remap1_Chaq[48, 3] <- remap1_Chaq[48, 3] + remap1_Chaq[49, 3]
remap1_Chaq <- remap1_Chaq[-c(49), ]

remap1_Chaq[56, 3] <- remap1_Chaq[56, 3] + remap1_Chaq[57, 3]
remap1_Chaq <- remap1_Chaq[-c(57), ]

remap3_Chaq[20, 3] <- remap3_Chaq[20, 3] + remap3_Chaq[21, 3]
remap3_Chaq <- remap3_Chaq[-c(21), ]

remap3_Chaq[21, 3] <- remap3_Chaq[21, 3] + remap3_Chaq[22, 3]
remap3_Chaq <- remap3_Chaq[-c(22), ]

remap3_Chaq[27, 3] <- remap3_Chaq[27, 3] + remap3_Chaq[28, 3]
remap3_Chaq <- remap3_Chaq[-c(28), ]

remap3_Chaq[33, 3] <- remap3_Chaq[33, 3] + remap3_Chaq[34, 3]
remap3_Chaq <- remap3_Chaq[-c(34), ]

remap3_Chaq[34, 3] <- remap3_Chaq[34, 3] + remap3_Chaq[35, 3]
remap3_Chaq <- remap3_Chaq[-c(35), ]

#join the df's - keeping all the remapping 3 and all other mapping in remapping 1
RNA1_remap1 <- anti_join(remap1_RNA1, remap3_RNA1, by = c("id", "segment"))
all_RNA1 <- full_join(remap3_RNA1, RNA1_remap1, by = c("id", "num_reads", "segment"))
all_RNA1 <- all_RNA1[-c(35), ]

RNA2_remap1 <- anti_join(remap1_RNA2, remap3_RNA2, by = c("id", "segment"))                    
all_RNA2 <- full_join(remap3_RNA2, RNA2_remap1, by = c("segment", "id", "num_reads"))

RNA3_remap1 <- anti_join(remap1_RNA3, remap3_RNA3, by = c("id", "segment"))                    
all_RNA3 <- full_join(remap3_RNA3, RNA3_remap1, by = c("segment", "id", "num_reads"))

Chaq_remap1 <- anti_join(remap1_Chaq, remap3_Chaq, by = c("id", "segment"))                    
all_Chaq <- full_join(remap3_Chaq, Chaq_remap1, by = c("segment", "id", "num_reads"))
all_Chaq <- all_Chaq[-c(57), ]

# Read in metadata table
metadata <- read_csv("analyses/data/Metadata_Table.csv") %>% 
  mutate(id = sample_name) %>% 
  select(id, genotype, chaq_present, infection_status, clade) %>% 
  mutate(id = str_replace_all(id, "_", "-"))

#combine all the data
all_RNA1 <- left_join(all_RNA1, metadata, by = "id")
all_RNA1 <- left_join(all_RNA1, all_counts, by = "id")

all_RNA2 <- left_join(all_RNA2, metadata, by = "id")
all_RNA2 <- left_join(all_RNA2, all_counts, by = "id")

all_RNA3 <- left_join(all_RNA3, metadata, by = "id")
all_RNA3 <- left_join(all_RNA3, all_counts, by = "id")

all_Chaq <- left_join(all_Chaq, metadata, by = "id")
all_Chaq <- left_join(all_Chaq, all_counts, by = "id")

# Calculate RPM
all_RNA1 <- all_RNA1 %>% 
  mutate(rpm = (num_reads/count)*1000000)

all_RNA2 <- all_RNA2 %>% 
  mutate(rpm = (num_reads/count)*1000000)

all_RNA3 <- all_RNA3 %>% 
  mutate(rpm = (num_reads/count)*1000000)

all_Chaq <- all_Chaq %>% 
  mutate(rpm = (num_reads/count)*1000000)

## Compare qPCR w/ sequencing abundance
# Read in qPCR data
qPCR_metadata <- read_csv("analyses/data/Metadata_Table.csv") %>% 
  select(sample_name, ct_2165_2170, ct_rpl) %>% 
  mutate(delta_ct = ct_2165_2170 - ct_rpl,
         relative = 2^-delta_ct,
         id = sample_name,
         id = str_replace_all(id, "_", "-")) 

# Combine dfs
all_RNA1_q <- left_join(all_RNA1, qPCR_metadata, by = "id" )
  
# Regression
reg <- glm(relative ~ num_reads, data = all_RNA1_q)
summary(reg)

# Make plot
plot <- ggplot(all_RNA1_q, aes(x = num_reads, y = relative)) +
  geom_point(alpha = 0.75, size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  stat_poly_line(aes(x = num_reads, relative), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", linewidth = 1) +
  stat_poly_eq(aes(x = num_reads, relative, label = paste(after_stat(eq.label),
                                                                  after_stat(rr.label), sep = "*\", \"*"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(x = "log10 Reads Per Million \n(RNA1)", y = "log10 Relative Galbut Virus Levels 
      \n (normalized to RpL32 mRNA)", color = "Chaq virus present")

plot 
ggsave("analyses/plots/RPM_vs_qPCR.pdf", units = "in", width = 10, height = 8)  

## Combine all segments

all <- rbind(all_RNA1, all_RNA2)
all <- rbind(all, all_RNA3)
all <- rbind(all, all_Chaq)

all <- all %>% 
  filter(infection_status == "sin")

A_vs_B <- ggplot(all, aes(x = clade, y = rpm)) +
  geom_boxplot() +
  geom_jitter(aes(color = clade), height = 0, width = 0.15, size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("darkolivegreen", "slateblue4")) +
  scale_y_log10() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(y = "Total Galbut Virus Reads Per Million \n(log10)", x = "Galbut Virus Clade",
       color = "Galbut Virus Clade")

A_vs_B  
ggsave("analyses/plots/clade_vs_rpm.pdf", units = "in", width = 10, height = 8) 

all_A <- all %>% 
  filter(clade == "A")

Chaq_vs_noChaq <- ggplot(all_A, aes(x = chaq_present, y = rpm)) +
  geom_boxplot() +
  geom_jitter(aes(color = chaq_present), height = 0, width = 0.15, size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("steelblue4", "orchid4")) +
  scale_y_log10() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(y = "Total Galbut Virus Reads Per Million \n(log10)", x = "Chaq Presence",
       color = "Chaq presence")

Chaq_vs_noChaq 
ggsave("analyses/plots/chaq_vs_noChaq.pdf", units = "in", width = 10, height = 8) 

# bi/trimodal infection phenotype?
all_wide <- pivot_wider(all, names_from = segment, values_from = c(num_reads, count, rpm)) %>% 
  mutate(num_reads_Chaq = if_else(is.na(num_reads_Chaq), 0, num_reads_Chaq)) %>% 
  mutate(num_reads_RNA1 = if_else(is.na(num_reads_RNA1), 0, num_reads_RNA1)) %>% 
  mutate(num_reads_RNA3 = if_else(is.na(num_reads_RNA3), 0, num_reads_RNA3)) %>% 
  mutate(num_reads_RNA2 = if_else(is.na(num_reads_RNA2), 0, num_reads_RNA2)) %>% 
  mutate(count_RNA3 = if_else(is.na(count_RNA3), count_RNA2, count_RNA3)) %>% 
  mutate(num_total = (num_reads_Chaq + num_reads_RNA1 + num_reads_RNA2 + num_reads_RNA3)) %>% 
  mutate(rpm_total = ((num_total/count_RNA3)*1000000))

phenotype <- ggplot(all_wide, aes(x = rpm_total)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap(~clade) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(y = "Number of samples", x = "Total Number of Galbut Virus Reads per Million (log10)")

phenotype            

phenotype_RNA <- ggplot(all, aes(x = rpm)) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(clade~segment) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(y = "Number of samples", x = "Number of Galbut Virus Reads per Million (log10)")

phenotype_RNA
ggsave("analyses/plots/phenotype_RNA.pdf", units = "in", width = 10, height = 8) 

# Total RMP for all segments
all_rpm <- rbind(all_RNA1, all_RNA2)
all_rpm <- rbind(all_rpm, all_RNA3)
all_rpm <- rbind(all_rpm, all_Chaq)

all_wide_rpm <- pivot_wider(all_rpm, names_from = segment, values_from = c(num_reads, count, rpm)) %>% 
  mutate(num_reads_Chaq = if_else(is.na(num_reads_Chaq), 0, num_reads_Chaq)) %>% 
  mutate(num_reads_RNA1 = if_else(is.na(num_reads_RNA1), 0, num_reads_RNA1)) %>% 
  mutate(num_reads_RNA3 = if_else(is.na(num_reads_RNA3), 0, num_reads_RNA3)) %>% 
  mutate(num_reads_RNA2 = if_else(is.na(num_reads_RNA2), 0, num_reads_RNA2)) %>% 
  mutate(count_RNA3 = if_else(is.na(count_RNA3), count_RNA2, count_RNA3)) %>% 
  mutate(num_total = (num_reads_Chaq + num_reads_RNA1 + num_reads_RNA2 + num_reads_RNA3)) %>% 
  mutate(rpm_total = ((num_total/count_RNA3)*1000000))

write_csv(all_wide_rpm, file = "analyses/data/RPM.csv")
