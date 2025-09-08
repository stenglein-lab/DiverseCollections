library(tidyverse)
library(readxl)
library(ggpmisc)
library(car)
library(rstatix)
library(ggpubr)

all_counts <- read_xlsx("analyses/RPM/all_read_counts.xlsx") %>% 
  filter(count_type == "post_trimming") %>% 
  rename(id = sample_id) %>% 
  select(id, count)

remapping_in <- read_xlsx("analyses/remapping_seqs/validation/num_aligned.xlsx")

#parse to get sample ids
remapping <- remapping_in %>% 
  mutate(sample_id = str_replace(sample_id, "sim-20-TD-4.filt", "D_mel_USA_2020_sim-20-TD-4_Galbut_virus_RNA3_Complete_CDS_2")) %>% 
  separate(sample_id, into = c("D", "mel_sim", "loc", "year", "id", "galbut", "virus", "segment"), sep = "_") %>% 
  mutate(segment = str_replace(segment, "Complete", "Chaq"),
         segment = str_replace(segment, "Capsid", "RNA2"),
         segment = str_replace(segment, "RNA2s", "RNA2"),
         segment = str_replace(segment, "Partial", "Chaq")) %>% 
  select(id, segment, num_reads) 

# check <- cbind(remapping,  remapping_in)

remap_comp_all <- left_join(remapping, all_counts, by = "id") %>% 
  mutate(rpm = ((num_reads/count)*1000000))

remap_part <- read_xlsx("analyses/remapping_seqs/partials/num_aligned_partials.xlsx") %>% 
  filter(align_seq != "MW976848",
         num_reads >= 10)

remap_part_wide <- pivot_wider(remap_part, names_from = c(align_seq, segment), values_from = num_reads) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(chaq_total = (`D_mel_USA_2023_1428-M-42_Chaq_virus_Complete_CDS_Chaq` + 
                         `D_mel_USA_2023_1020-M-0818-A7_Chaq_virus_Complete_CDS_Chaq` +
                         `D_mel_USA_2023_OH-M-5_Chaq_virus_Complete_CDS_Chaq`),
         RNA1_total = (MT742164_1),
         RNA2_total = (`D_mel_USA_2023_CVID-1006-G7_Galbut_virus_RNA2_Complete_CDS_2`),
         RNA3_total = (`D_mel_USA_2023_500-F-41_Galbut_virus_RNA3_Complete_CDS_1_3` +
                         `D_mel_USA_2023_1428-M-42_Galbut_virus_RNA3_Complete_CDS_3` +
                         `D_mel_USA_2023_Adobe-B4_Galbut_virus_RNA3_Complete_CDS_3`)) %>% 
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

# tally counts in cases of co-infection
remap_all <- remap_all %>% group_by(id, count, segment) %>% summarize(num_reads = sum(num_reads))

#calculate rpm
remap_all <- remap_all %>% 
  mutate(rpm = ((num_reads/count)*1000000))

## Compare qPCR w/ sequencing abundance
# Read in qPCR data
qPCR_metadata <- read_xlsx("analyses/data/Metadata_Table.xlsx") %>% 
  select(sample_name, ct_2165_2170, ct_rpl, infection_status, clade, chaq_present) %>% 
  mutate(id = sample_name,
         id = str_replace_all(id, "_", "-")) 

#combine all the data
remap_all <- left_join(remap_all, qPCR_metadata, by = "id") %>% 
  mutate(delta_ct = ct_2165_2170 - ct_rpl,
         relative = 2^-delta_ct)

# Regression
reg <- glm(log10(relative) ~ log10(num_reads), data = filter(remap_all, segment == "RNA1"))
summary(reg)

# Make plot
rpm_plot <- ggplot(filter(remap_all, segment == "RNA1")) +
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

rpm_plot 
ggsave("analyses/plots/RPM_vs_qPCR.pdf", units = "in", width = 10, height = 8)  

# Get total RPM
remap_all_wide <- remap_all %>% 
  select(id, segment, rpm, infection_status, clade, chaq_present) %>% 
  pivot_wider(names_from = segment, values_from = rpm, values_fill = 0) %>% 
  mutate(rpm_total = (RNA2 + RNA1 + RNA3 + Chaq))

all_rpm <- remap_all_wide %>% 
  select(id, rpm_total) %>% 
  mutate(id = str_replace(id, "-", "_")) %>% 
  rename(sample_name = id)

write_csv(all_rpm, "analyses/data/all_rpm_counts.csv")

# Numbers
avg_a_vs_b <- remap_all_wide %>% 
  filter(infection_status == "sin") %>% 
  group_by(clade) %>% 
  summarise(mean = mean(rpm_total),
            median = median(rpm_total))

# Filter coinfections out
remap_no_co <- remap_all_wide %>% 
  filter(infection_status == "sin") %>% 
  filter(clade != "-")

#Wilcoxon test to compare genotype rpms
geno_wilcox <- remap_no_co %>% 
  ungroup() %>% 
  wilcox_test(rpm_total ~ clade)
geno_wilcox_plot <- geno_wilcox %>% add_xy_position(x = "clade")


A_vs_B <- ggplot(remap_no_co, aes(x = clade, y = rpm_total)) +
  geom_boxplot() +
  geom_jitter(aes(color = clade), height = 0, width = 0.15, size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("darkolivegreen", "slateblue4")) +
  scale_y_log10() +
  stat_pvalue_manual(geno_wilcox_plot, tip.length = 0.01, y.position = 5,
                     bracket.shorten=0.1, size=3, bracket.nudget.y=0.1) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1)) +
  labs(y = "Total Galbut Virus Reads Per Million", x = "Galbut Virus Clade",
       color = "Galbut Virus Clade")

A_vs_B  
ggsave("analyses/plots/clade_vs_rpm.pdf", units = "in", width = 10, height = 8) 

# Wilcoxon test to compare sample rpms
all_A <- remap_all_wide %>% 
  filter(clade == "A") %>% 
  filter(infection_status == "sin")

chaq_wilcox <- all_A %>% 
  ungroup() %>% 
  wilcox_test(rpm_total ~ chaq_present)
chaq_wilcox_plot <- chaq_wilcox %>% add_xy_position(x = "chaq_present")

# Levene's Test to compare variance
lev <- leveneTest(rpm_total ~ chaq_present, all_A)
lev

Chaq_vs_noChaq <- ggplot(all_A, aes(x = chaq_present, y = rpm_total)) +
  geom_boxplot() +
  geom_jitter(aes(color = chaq_present), height = 0, width = 0.15, size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("steelblue4", "orchid4")) +
  scale_y_log10() +
  stat_pvalue_manual(chaq_wilcox_plot, tip.length = 0.01, y.position = 5,
                     bracket.shorten=0.1, size=3, bracket.nudget.y=0.1) +
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

# do RPM counts differ by segment?

# first, reorder so chaq will be plotted after galbut virus segments
remap_comp_all$segment <- fct_relevel(remap_comp_all$segment, "Chaq", after=3)

# this Shapiro tests indicates fractional counts are not normally distributed
df_shapiro     <- remap_comp_all %>% group_by(segment) %>% shapiro_test(rpm)
filter(df_shapiro, p<0.05)

# not normally distributed, so use Wilcoxon test
df_wilcox <- remap_comp_all %>% wilcox_test(rpm ~ segment)
df_wilcox_plot <- df_wilcox %>% add_xy_position(x = "segment")

phenotype_RNA_jitter <- ggplot(remap_comp_all) + 
  # these lines will connect points from the same sample
  # geom_line(aes(x = segment, y = rpm, group=id), 
  # color="slateblue", linewidth=0.25, alpha=0.25) +
  geom_jitter(aes(x = segment, y = rpm),
              shape=21, fill="darkslateblue", color="black", size=2, stroke=0.25, 
              alpha=0.5,
              width=0.20, height=0) +
  geom_boxplot(aes(x = segment, y = rpm), 
               color="darkslateblue", outlier.shape=NA, fill=NA,
               width=0.5) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
  scale_y_log10() +
  labs(y = "galbut virus mapping reads per million reads", x = "")

phenotype_RNA_jitter
# ggsave("analyses/plots/RPM_per_segment.pdf", units = "in", width = 10, height = 8)

y_positions <- rep(c(5,5.2),1)
phenotype_RNA_jitter_stats <- phenotype_RNA_jitter + 
  stat_pvalue_manual(filter(df_wilcox_plot, p.adj.signif != "ns"), 
                     y.position = y_positions, tip.length = 0.01, 
                     bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)
phenotype_RNA_jitter_stats 
ggsave("analyses/plots/RPM_per_segment.pdf", units = "in", width = 10, height = 8)

# create text for paper

rpm_avgs <- remap_comp_all %>% group_by(segment) %>% summarize(median_rpm = median(rpm))
RNA1_median <- rpm_avgs %>% filter(segment == "RNA1") %>% pull(median_rpm)
RNA2_median <- rpm_avgs %>% filter(segment == "RNA2") %>% pull(median_rpm)
RNA3_median <- rpm_avgs %>% filter(segment == "RNA3") %>% pull(median_rpm)
Chaq_median <- rpm_avgs %>% filter(segment == "Chaq") %>% pull(median_rpm)

RNA1_2_p <- filter(df_wilcox_plot, group1 == "RNA1" & group2 == "RNA2") %>% pull(p.adj)
RNA1_3_p <- filter(df_wilcox_plot, group1 == "RNA1" & group2 == "RNA3") %>% pull(p.adj)

paper_text <- paste0(
  "This was consistent with lower average coverage levels for RNA1 (Fig. X). ",
  "The median level of RNA1 mapping reads per million total reads (RPM) was ",
  sprintf("%0.0f", RNA1_median),
  ". ",
  "This was significantly lower than the median RPM values for RNA2 and RNA3 of ",
  sprintf("%0.0f", RNA2_median),
  " and ",
  sprintf("%0.0f", RNA3_median),
  " respectively (Wilcoxon p adjusted = ",
  sprintf("%0.3f", RNA1_2_p),
  " and ",
  sprintf("%0.3f", RNA1_3_p),
  "; Fig. X). "
)
paper_text
