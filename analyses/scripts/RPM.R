library(tidyverse)
library(readxl)
library(ggpmisc)
library(car)

all_counts <- read_xlsx("analyses/RPM/all_read_counts.xlsx") %>% 
  filter(count_type == "post_trimming") %>% 
  rename(id = sample_id) %>% 
  select(id, count)

remapping <- read_xlsx("analyses/remapping_seqs/validation/num_aligned.xlsx")

#parse to get sample ids
remapping <- remapping %>% 
  mutate(sample_id = str_replace(sample_id, "sim-20-TD-4.filt", "D_mel_USA_2020_sim-20-TD-4_Galbut_virus_RNA3_Complete_CDS_2")) %>% 
  separate(sample_id, into = c("D", "mel_sim", "loc", "year", "id", "galbut", "virus", "segment"), sep = "_") %>% 
  mutate(segment = str_replace(segment, "Complete", "Chaq"),
         segment = str_replace(segment, "Capsid", "RNA2"),
         segment = str_replace(segment, "RNA2s", "RNA2"),
         segment = str_replace(segment, "Partial", "Chaq")) %>% 
  select(id, segment, num_reads) 

remap_comp_all <- left_join(remapping, all_counts, by = "id") %>% 
  mutate(rpm = ((num_reads/count)*1000000))

remap_part <- read_xlsx("analyses/remapping_seqs/partials/num_aligned_partials.xlsx") %>% 
  filter(align_seq != "MW976848")

remap_part_wide <- pivot_wider(remap_part, names_from = c(align_seq, segment), values_from = num_reads) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(chaq_total = (`D_mel_USA_2023_1428-M-42_Chaq_virus_Complete_CDS_Chaq` + 
                         `D_mel_USA_2023_ME-M-8_Chaq_virus_Complete_CDS_Chaq` +
                         `D_mel_USA_2023_1020-M-0818-A7_Chaq_virus_Complete_CDS_Chaq` +
                         `D_mel_USA_2023_OH-M-5_Chaq_virus_Complete_CDS_Chaq` +
                         `D_mel_USA_2023_Penn-M-2_Chaq_virus_Complete_CDS_2_Chaq`),
         RNA1_total = (MT742164_1 + `D_mel_USA_2023_ME-M-7_Galbut_virus_RNA1_Complete_CDS_1_1` +
                         `D_mel_USA_2023_ME-F-1_Galbut_virus_RNA1_Complete_CDS_1`),
         RNA2_total = (`D_mel_USA_2023_CVID-1006-G7_Galbut_virus_RNA2_Complete_CDS_2` + 
                         `D_mel_USA_2023_ME-M-4_Galbut_virus_RNA2_Complete_CDS_2` +
                         `D_mel_USA_2023_500-M-G2_Galbut_virus_RNA2_Complete_CDS_2`),
         RNA3_total = (`D_mel_USA_2023_500-F-41_Galbut_virus_RNA3_Complete_CDS_1_3` +
                         `D_mel_USA_2023_1428-M-42_Galbut_virus_RNA3_Complete_CDS_3` +
                         `D_mel_USA_2023_Adobe-B4_Galbut_virus_RNA3_Complete_CDS_3` +
                         `D_mel_USA_2023_500-M-G2_Galbut_virus_RNA3_Complete_CDS_3` +
                         `MH384362_3` +
                         `D_mel_USA_2023_500-M-61-2_Galbut_virus_RNA3_Complete_CDS_3`)) %>% 
  select(id, chaq_total, RNA1_total, RNA2_total, RNA3_total) %>% 
  mutate(total = (chaq_total + RNA1_total + RNA2_total + RNA3_total)) 

remap_part_all <- left_join(remap_part_wide, all_counts, by = "id") 

remap_part_all <- pivot_longer(remap_part_all, cols = chaq_total:RNA3_total,
                               names_to = "segment", values_to = "num_reads") %>% 
  filter(num_reads != 0) %>% 
  mutate(segment = str_replace(segment, "chaq_total", "Chaq"),
         segment = str_replace(segment, "RNA1_total", "RNA1"),
         segment = str_replace(segment, "RNA2_total", "RNA2"),
         segment = str_replace(segment, "RNA3_total", "RNA3")) %>% 
  mutate(rpm = ((num_reads/count)*1000000)) %>% 
  mutate(rpm_all = ((total/count)*1000000))

remap_all <- full_join(remap_comp_all, remap_part_all, 
                       by = c("id", "count", "num_reads", "rpm", "segment")) %>% 
  select(id, count, num_reads, segment)

# Combine rows with two RNA1s
remap_all[284, 3] <- remap_all[284, 3] + remap_all[287, 3]
remap_all <- remap_all[-c(287), ]

remap_all[340, 3] <- remap_all[340, 3] + remap_all[342, 3]
remap_all <- remap_all[-c(342), ]


#Combine rows with two RNA2s
remap_all[100, 3] <- remap_all[100, 3] + remap_all[101, 3]
remap_all <- remap_all[-c(101), ]

remap_all[282, 3] <- remap_all[282, 3] + remap_all[284, 3]
remap_all <- remap_all[-c(284), ]

remap_all[337, 3] <- remap_all[337, 3] + remap_all[339, 3]
remap_all <- remap_all[-c(339), ]

remap_all[341, 3] <- remap_all[341, 3] + remap_all[345, 3]
remap_all <- remap_all[-c(345), ]

#Combine rows with two RNA3s
remap_all[3, 3] <- remap_all[3, 3] + remap_all[4, 3]
remap_all <- remap_all[-c(4), ]

remap_all[38, 3] <- remap_all[38, 3] + remap_all[40, 3]
remap_all <- remap_all[-c(40), ]

remap_all[79, 3] <- remap_all[79, 3] + remap_all[80, 3]
remap_all <- remap_all[-c(80), ]

remap_all[100, 3] <- remap_all[100, 3] + remap_all[103, 3]
remap_all <- remap_all[-c(103), ]

remap_all[120, 3] <- remap_all[120, 3] + remap_all[122, 3]
remap_all <- remap_all[-c(122), ]

remap_all[276, 3] <- remap_all[276, 3] + remap_all[279, 3]
remap_all <- remap_all[-c(279), ]

remap_all[281, 3] <- remap_all[281, 3] + remap_all[284, 3]
remap_all <- remap_all[-c(284), ]

remap_all[227, 3] <- remap_all[227, 3] + remap_all[229, 3]
remap_all <- remap_all[-c(229), ]

remap_all[301, 3] <- remap_all[301, 3] + remap_all[304, 3]
remap_all <- remap_all[-c(304), ]

remap_all[327, 3] <- remap_all[327, 3] + remap_all[330, 3]
remap_all <- remap_all[-c(330), ]

remap_all[356, 3] <- remap_all[356, 3] + remap_all[359, 3]
remap_all <- remap_all[-c(359), ]

remap_all[363, 3] <- remap_all[363, 3] + remap_all[365, 3]
remap_all <- remap_all[-c(365), ]

#Combine rows with two Chaqs
remap_all[213, 3] <- remap_all[213, 3] + remap_all[214, 3]
remap_all <- remap_all[-c(214), ]

remap_all[261, 3] <- remap_all[261, 3] + remap_all[263, 3]
remap_all <- remap_all[-c(263), ]

remap_all[317, 3] <- remap_all[317, 3] + remap_all[318, 3]
remap_all <- remap_all[-c(318), ]

remap_all[351, 3] <- remap_all[351, 3] + remap_all[355, 3]
remap_all <- remap_all[-c(355), ]

#calculate rpm
remap_all <- remap_all %>% 
  mutate(rpm = ((num_reads/count)*1000000))

## Compare qPCR w/ sequencing abundance
# Read in qPCR data
qPCR_metadata <- read_csv("analyses/data/Metadata_Table.csv") %>% 
  select(sample_name, ct_2165_2170, ct_rpl, infection_status, clade, chaq_present) %>% 
  mutate(id = sample_name,
         id = str_replace_all(id, "_", "-")) 

#combine all the data
remap_all <- left_join(remap_all, qPCR_metadata, by = "id") %>% 
  mutate(delta_ct = ct_2165_2170 - ct_rpl,
         relative = 2^-delta_ct)

# Regression
reg <- glm(relative ~ num_reads, data = filter(remap_all, segment == "RNA1"))
summary(reg)

# Make plot
plot <- ggplot(filter(remap_all, segment == "RNA1")) +
  geom_point(aes(x = num_reads, y = relative), alpha = 0.75, size = 3) +
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
  labs(x = "Number of galbut viurs RNA1-mapping reads per million total reads",
  y = "Galbut virus RNA levels relative to RpL32 mRNA", color = "Chaq virus present")

plot 
ggsave("analyses/plots/RPM_vs_qPCR.pdf", units = "in", width = 10, height = 8)  

# Get total RPM
remap_all_wide <- remap_all %>% 
  select(id, count, num_reads, segment, rpm, infection_status, clade, chaq_present) %>% 
  pivot_wider(names_from = segment, values_from = c(num_reads, count, rpm )) %>% 
  mutate(rpm_Chaq = if_else(is.na(rpm_Chaq), 0, rpm_Chaq)) %>% 
  mutate(rpm_RNA1 = if_else(is.na(rpm_RNA1), 0, rpm_RNA1)) %>% 
  mutate(rpm_RNA3 = if_else(is.na(rpm_RNA3), 0, rpm_RNA3)) %>% 
  mutate(rpm_RNA2 = if_else(is.na(rpm_RNA2), 0, rpm_RNA2)) %>% 
  mutate(rpm_total = (rpm_Chaq + rpm_RNA1 + rpm_RNA2 + rpm_RNA3))

# Filter coinfections out
remap_no_co <- remap_all_wide %>% 
  filter(infection_status == "sin")

A_vs_B <- ggplot(remap_no_co, aes(x = clade, y = rpm_total)) +
  geom_boxplot() +
  geom_jitter(aes(color = clade), height = 0, width = 0.15, size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("darkolivegreen", "slateblue4", "pink")) +
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

all_A <- remap_all_wide %>% 
  filter(clade == "A")

Chaq_vs_noChaq <- ggplot(all_A, aes(x = chaq_present, y = rpm_total)) +
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
phenotype <- ggplot(remap_all, aes(x = rpm)) +
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

phenotype_RNA <- ggplot(remap_all, aes(x = rpm)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap(~segment) +
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

# Numbers
avg_a_vs_b <- remap_all_wide %>% 
  filter(infection_status == "sin") %>% 
  group_by(clade) %>% 
  summarise(mean = mean(rpm_total),
            median = median(rpm_total))

# STATS
chaq_p <- all_A %>% 
  filter(chaq_present == "y") 
chaq_a <- all_A %>% 
  filter(chaq_present == "n")

# Wilcoxon test to compare sample rpms
wilco <- wilcox.test(chaq_a$rpm_total, chaq_p$rpm_total)
wilco

# Levene's Test to compare variance
lev <- leveneTest(rpm_total ~ chaq_present, all_A)
lev
