library(tidyverse)
library(readxl)
library(ggpmisc)
library(rstatix)
library(ggpubr)
library(scales)
library(patchwork)

# read in host mapping read counts
host_mapping_counts <- read.delim("analyses/data/host_mapping_read_counts.txt", sep="\t", header=F)
colnames(host_mapping_counts) <- c("id", "host_mapping_reads")

# read in RpL32 mRNA mapping reads counts
rpl32_mapping_counts <- read.delim("analyses/data/RpL32_mapping_reads_per_dataset.txt", sep="\t", header=F)
# output format from samtools coverage
colnames(rpl32_mapping_counts) <- c("id", "accession", "startpos", "endpos", 
                                    "numreads", "covbases", "coverage", 
                                    "meandepth",  "meanbaseq", "meanmapq")

# simplify
rpl32_mapping_counts <- rpl32_mapping_counts %>% select(id, numreads) %>% rename(rpl32_mapping_counts = numreads)

# RpL32 mRNA length (NM_079843)
rpl32_mRNA_length <- 578

# galbut virus RNA 1 length (e.g., MT742160) 
galbut_RNA1_length <- 1685

# per refseq read mapping counts
refseq_mapping_counts <- read.delim("analyses/data/collected_per_refseq_coverage.tsv", header=F, sep="\t")
# output format from samtools coverage
colnames(refseq_mapping_counts) <- c("id", "accession", "startpos", "endpos", 
                                     "numreads", "covbases", "coverage", 
                                     "meandepth",  "meanbaseq", "meanmapq")

# join in reference sequence metadata
refseq_metadata <- read.delim("analyses/data/sample_accession_map.txt", sep="\t", header=F)
colnames(refseq_metadata) <- c("id", "segment", "accession", "accession_description")

# read in clade labels for all acccessions (melA, simA, etc)
clade_labels <- read.delim("analyses/data/accession_clade_labels.txt", sep="\t", header=T)

refseq_mapping_counts <- left_join(refseq_mapping_counts, refseq_metadata, by = join_by(id, accession))
refseq_mapping_counts <- left_join(refseq_mapping_counts, host_mapping_counts, by = join_by(id))
refseq_mapping_counts <- left_join(refseq_mapping_counts, rpl32_mapping_counts, by = join_by(id))
refseq_mapping_counts <- left_join(refseq_mapping_counts, clade_labels , by = join_by(accession))

# calculate RPM: galbut virus mapping reads per million host-mapping reads
# and RPR: galbut virus mapping reads per RpL32 mapping read
# and rprk = galbut virus mapping reads per kbp of galbut virus RNA1 length per RPL32 mapping read per kbp of RpL32 length
refseq_mapping_counts <- refseq_mapping_counts %>% mutate(rpm = 1e6 * numreads / host_mapping_reads,
                                                          rpr = numreads / rpl32_mapping_counts,
                                                          rprk = (numreads / galbut_RNA1_length) / 
                                                            (rpl32_mapping_counts / rpl32_mRNA_length))
# calculate per sample sums
refseq_mapping_counts_sums <- 
  refseq_mapping_counts %>% 
  group_by(id) %>% 
  summarize(summed_rpm      = sum(rpm))

# pull out just RNA1 
rna1_mapping_counts <- refseq_mapping_counts %>% filter(segment == "RNA1") 

# sum RPM and RPR values to handle co-infected samples
rna1_mapping_counts <- rna1_mapping_counts %>% group_by(id) %>% summarize(rpm = sum(rpm),
                                                                          rpr = sum(rpr),
                                                                          rprk = sum(rprk))
                                                                          
# Read in qPCR data
# and calculate relative galbut virus RNA from Cts
qPCR_metadata <- read_xlsx("analyses/data/Metadata_Table.xlsx") %>% 
  select(sample_name, ct_2165_2170, ct_rpl, infection_status, clade, chaq_present) %>% 
  mutate( id = str_replace_all(sample_name, "_", "-"),
          delta_ct = ct_2165_2170 - ct_rpl,
          relative = 2^-delta_ct)
 
# join RPM and qPCR data
remap_all <- left_join(rna1_mapping_counts, qPCR_metadata, by = "id") 

# Linear regression: qPCR relative RNA vs NGS RPM 
# reg <- glm(relative ~ num_reads, data = filter(remap_all, segment == "RNA1"))
reg <- glm(relative ~ rpm, data = remap_all)
summary(reg)

# vs. just RpL32 mapping reads?
reg_rpr <- glm(relative ~ rpr, data = filter(remap_all, is.finite(rpr)))
summary(reg_rpr)

# vs. just RpL32 mapping reads?
reg_rpr <- glm(relative ~ rpr, data = filter(remap_all, is.finite(rpr)))
summary(reg_rpr)

# vs. just RpL32 mapping reads normalized to length of RNAs?
reg_rprk <- glm(relative ~ rprk, data = filter(remap_all, is.finite(rprk)))
summary(reg_rprk)

# Make RPM plot
rpm_plot <- ggplot(remap_all) +
  geom_point(aes(x = rpm, y = relative), alpha = 0.75, size = 3) +
  scale_x_log10(labels = trans_format("log10", label_math())) +
  scale_y_log10(labels = trans_format("log10", label_math())) +
  stat_poly_line(aes(x = rpm, relative), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", linewidth = 1) +
  stat_poly_eq(aes(x = rpm, relative, 
                   label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_bw(base_size = 11) +
  labs( x = "galbut virus RNA1-mapping reads per million host-mapping reads (NGS)",
        y = "Galbut virus RNA levels relative to RpL32 mRNA\n(RT-qPCR)")

rpm_plot 
# ggsave("analyses/plots/RPM_vs_qPCR.pdf", units = "in", width = 10, height = 8)  

# plot vs. RpL32 normalized read
rpr_plot <- ggplot(filter(remap_all, is.finite(rprk))) +
  geom_point(aes(x = rprk, y = relative), alpha = 0.75, size = 3) +
  scale_x_log10(labels = trans_format("log10", label_math())) +
  scale_y_log10(labels = trans_format("log10", label_math())) +
  stat_poly_line(aes(x = rprk, relative), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", linewidth = 1) +
  stat_poly_eq(aes(x = rprk, relative, 
                   label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_bw(base_size = 12) +
  labs( x = "Galbut virus RNA1-mapping reads per RpL32-mapping read (NGS)",
        y = "Galbut virus RNA levels relative to RpL32 mRNA\n(RT-qPCR)")

rpr_plot 
# ggsave("analyses/plots/RPR_vs_qPCR.pdf", units = "in", width = 10, height = 8)  

# plot RpL32 mapping reads vs total host mapping reads

# save a combined supplemental figure
rpr_plot + rpm_plot + plot_layout(ncol = 1) + plot_annotation(tag_levels = 'A')
ggsave("analyses/plots/Supplemental_Figure_2_qPCR_vs_NGS_quantification.pdf", units="in", width=7.5, height=10)

# use data to determine:
# - whether co-infection or not (>1 of any segment)

coinfected_samples <- 
  refseq_mapping_counts %>% 
  group_by(id, segment) %>%
  summarize(seqs_per_segment = n(), 
            accessions = paste0(accession, collapse=","),
            .groups="drop") %>%
  filter(seqs_per_segment > 1) 
  
coinfected_sample_ids <-
  coinfected_samples %>%
  group_by(id) %>%
  summarize() %>%
  pull(id)

# use data to determine:
# - whether chaq is present
has_chaq_sample_ids <- 
  refseq_mapping_counts %>%
  filter(segment == "Chaq") %>%
  group_by(id) %>%
  summarize() %>%
  pull(id)


# Filter coinfections out
refseq_mapping_counts_no_coinfection <-
  refseq_mapping_counts %>%
  filter (! id %in% coinfected_sample_ids)

# double check all sequences from single infections assigned to same clades
discordant_clade_assignments <- 
  refseq_mapping_counts_no_coinfection %>% 
  group_by(id, clade_label) %>%
  summarize(.groups="drop") %>%
  group_by(id) %>%
  summarize(number_clades = n()) %>%
  filter(number_clades > 1)

# if any issue will error here
stopifnot(nrow(discordant_clade_assignments) == 0)

# sum RPM values for each sample
refseq_mapping_counts_no_coinfection_sums <- 
  refseq_mapping_counts_no_coinfection %>% 
  group_by(id, clade_label) %>%
  summarize(summed_rpm = sum(rpm), .groups="drop")
  
# how many of each clade?
refseq_mapping_counts_no_coinfection_sums %>%
  group_by(clade_label) %>% 
  summarize(n = n())

# Wilcoxon test to compare genotype rpms
geno_wilcox <- refseq_mapping_counts_no_coinfection_sums %>% 
  wilcox_test(summed_rpm ~ clade_label)

geno_wilcox_plot <- geno_wilcox %>% add_xy_position(x = "clade_label") %>%
  mutate(p.signif = case_when(
    p > 0.05 ~ "ns",
    p > 0.01 ~ "*",
    p > 0.001 ~ "**",
    p > 0.0001 ~ "***",
    .default = "****",
  ))

A_vs_B <- ggplot(refseq_mapping_counts_no_coinfection_sums, 
                 aes(x = clade_label, y = summed_rpm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = clade_label), shape=21, color="black", stroke=0.25, height = 0, width = 0.15, size = 2, alpha = 0.85) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue4")) +
  scale_color_manual(values = c("darkolivegreen", "slateblue4")) +
  scale_y_log10() +
  stat_pvalue_manual(geno_wilcox_plot, tip.length = 0.01, y.position = 6, label="p.signif",
                     bracket.shorten=0.1, size=3, bracket.nudget.y=0.1) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        # axis.text = element_text(size = 12),
        legend.position = "none"
  ) +
  labs(y = "Galbut virus-mapping reads\nper million host-mapping reads", x = "Galbut virus clade")

A_vs_B  
# ggsave("analyses/plots/clade_vs_rpm.pdf", units = "in", width = 10, height = 8) 

# calculate medians per group
clade_median_rpm <-
  refseq_mapping_counts_no_coinfection_sums %>% 
  group_by(clade_label) %>%
  summarize(median_rpm = median(summed_rpm))

melA_median <- clade_median_rpm %>% filter(clade_label == "melA") %>% pull(median_rpm)
melB_median <- clade_median_rpm %>% filter(clade_label == "melB") %>% pull(median_rpm)

# output some text in paper re: levels of A vs B
paste0(
  "Galbut virus genotype melA median RPM values were ",
  sprintf("%0.1fx ", melA_median / melB_median),
  "times higher than in genotype melB infections (Fig. 8A, Wilcoxon p-value = ",
  sprintf("%0.1e", geno_wilcox_plot$p),
  "). "
)

# is presence of chaq associated with differences in galbut virus levels?  
refseq_mapping_counts_no_coinfection_galbut_only_sums <-
  refseq_mapping_counts_no_coinfection %>%
  filter(clade_label == "melA") %>%
  filter(segment != "Chaq") %>%
  group_by(id) %>%
  summarize(summed_galbut_rpm = sum(rpm))

# do the samples have a chaq sequence?
refseq_mapping_counts_no_coinfection_galbut_only_sums <- 
  refseq_mapping_counts_no_coinfection_galbut_only_sums %>%
  mutate(chaq_present = if_else(id %in% has_chaq_sample_ids, "Yes", "No"))

# Wilcoxon test to compare sample rpms
chaq_wilcox <- refseq_mapping_counts_no_coinfection_galbut_only_sums %>% 
  wilcox_test(summed_galbut_rpm ~ chaq_present)
chaq_wilcox_plot <- chaq_wilcox %>% add_xy_position(x = "chaq_present") %>%
  mutate(p.signif = case_when(
    p > 0.05 ~ "ns",
    p > 0.01 ~ "*",
    p > 0.001 ~ "**",
    p > 0.0001 ~ "***",
    .default = "****",
  ))

Chaq_vs_noChaq <- ggplot(refseq_mapping_counts_no_coinfection_galbut_only_sums, 
                         aes(x = chaq_present, y = summed_galbut_rpm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = chaq_present), shape=21, stroke=0.25, color="black", height = 0, width = 0.15, size = 2, alpha = 0.85) +
  scale_fill_manual(values = c("steelblue4", "orchid4")) +
  scale_color_manual(values = c("steelblue4", "orchid4")) +
  scale_y_log10(breaks=c(1e2,1e3,1e4,1e5,1e6)) + 
  stat_pvalue_manual(chaq_wilcox_plot, tip.length = 0.01, y.position = 6, label = "p.signif",
                     bracket.shorten=0.1, size=3, bracket.nudget.y=0.1) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        legend.position = "none",
        # axis.text = element_text(size=12)
        ) +
  labs(y = "Galbut virus-mapping reads\nper million host-mapping reads", x = "Chaq virus present")

Chaq_vs_noChaq 

# ggsave("analyses/plots/chaq_vs_noChaq.pdf", units = "in", width = 10, height = 8) 

# output text for paper re: chaq presence

paste0("There was no difference in median galbut virus RPM levels between infections ",
       "with or without chaq virus (Fig. 6B; Wilcoxon p-value = ",
       sprintf("%0.2f", chaq_wilcox$p)
       ) 

# create and output a combined plot
A_vs_B + Chaq_vs_noChaq + plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")

ggsave("analyses/plots/Figure_6_galbut_levels_by_genotype_and_chaq.pdf", width=5, height=4, units="in")

# total RPMs in co-infected samples vs. singly infected?
refseq_mapping_counts_sums <- 
  refseq_mapping_counts_sums %>%
  mutate(is_coinfected = if_else(id %in% coinfected_sample_ids, "Yes", "No")) 

coinf_wilcox <- refseq_mapping_counts_sums %>% 
  wilcox_test(summed_rpm ~ is_coinfected ) %>%
  add_xy_position(x = "is_coinfected")

coinfected_vs_single <- 
  ggplot(refseq_mapping_counts_sums, 
                 aes(x = is_coinfected, y = summed_rpm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = is_coinfected), height = 0, width = 0.15, size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("darkolivegreen", "slateblue4")) +
  scale_y_log10() +
  stat_pvalue_manual(coinf_wilcox, tip.length = 0.01, y.position = 6,
                     bracket.shorten=0.1, size=3, bracket.nudget.y=0.1) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 12),
        legend.position = "none"
  ) +
  labs(y = "Galbut virus-mapping reads\nper million host-mapping reads", x = "Coinfected sample")

coinfected_vs_single  
ggsave("analyses/plots/Supplemental_Figure_X_singly_vs_coinfected_rpm.pdf", width=7.5, height=5, units="in")
