all_seg_chart <- ggplot(metadata, aes(num_seqs_total, fill = location, color = location)) +
  geom_histogram(binwidth = 1, bins = 8, alpha = 0.75) +
  scale_fill_manual(values = c("indianred3", "orangered3", "plum4", "darkseagreen", "seagreen")) +
  scale_color_manual(values = c("indianred4", "orangered4", "pink4", "darkseagreen4", "seagreen4")) +
  theme_bw() +
  labs(x = "Number of sequences in a sample", y = "", fill = "Location", color = "Location")

all_seg_chart
ggsave("analyses/plots/galbut_segment_histogram_all.pdf", units = "in", width = 14, height = 8)

all_seg_chart2 <- ggplot(filter(metadata, num_seqs > 0), aes(num_seqs, fill = location, color = location)) +
  geom_histogram(binwidth = 1, bins = 8, alpha = 0.75) +
  scale_fill_manual(values = c("indianred3", "orangered3", "plum4", "darkseagreen", "seagreen")) +
  scale_color_manual(values = c("indianred4", "orangered4", "pink4", "darkseagreen4", "seagreen4")) +
  theme_bw() +
  labs(x = "Number of sequences in a sample", y = "", fill = "Location", color = "Location")

all_seg_chart2
ggsave("analyses/plots/galbut_segment_histogram.pdf", units = "in", width = 14, height = 8)

RNA1_chart <- ggplot(metadata, aes(RNA1, fill = location, color = location)) +
  geom_histogram(binwidth = 1, bins = 3, alpha = 0.75) +
  scale_fill_manual(values = c("indianred3", "orangered3", "plum4", "darkseagreen", "seagreen")) +
  scale_color_manual(values = c("indianred4", "orangered4", "pink4", "darkseagreen4", "seagreen4")) +
  theme_bw() +
  labs(x = "Number of RNA1 sequences in a sample", y = "", fill = "Location", color = "Location")

RNA1_chart
ggsave("analyses/plots/RNA1_histogram.pdf", units = "in", width = 14, height = 8)

RNA2_chart <- ggplot(metadata, aes(RNA2, fill = location, color = location)) +
  geom_histogram(binwidth = 1, bins = 3, alpha = 0.75) +
  scale_fill_manual(values = c("indianred3", "orangered3", "plum4", "darkseagreen", "seagreen")) +
  scale_color_manual(values = c("indianred4", "orangered4", "pink4", "darkseagreen4", "seagreen4")) +
  theme_bw() +
  labs(x = "Number of RNA2 sequences in a sample", y = "", fill = "Location", color = "Location")

RNA2_chart
ggsave("analyses/plots/RNA2_histogram.pdf", units = "in", width = 14, height = 8)

RNA3_chart <- ggplot(metadata, aes(RNA3, fill = location, color = location)) +
  geom_histogram(binwidth = 1, bins = 3, alpha = 0.75) +
  scale_fill_manual(values = c("indianred3", "orangered3", "plum4", "darkseagreen", "seagreen")) +
  scale_color_manual(values = c("indianred4", "orangered4", "pink4", "darkseagreen4", "seagreen4")) +
  theme_bw() +
  labs(x = "Number of RNA3 sequences in a sample", y = "", fill = "Location", color = "Location")

RNA3_chart
ggsave("analyses/plots/RNA3_histogram.pdf", units = "in", width = 14, height = 8)

Chaq_chart <- ggplot(metadata, aes(Chaq, fill = location, color = location)) +
  geom_histogram(binwidth = 1, bins = 3, alpha = 0.75) +
  scale_fill_manual(values = c("indianred3", "orangered3", "plum4", "darkseagreen", "seagreen")) +
  scale_color_manual(values = c("indianred4", "orangered4", "pink4", "darkseagreen4", "seagreen4")) +
  theme_bw() +
  labs(x = "Number of Chaq sequences in a sample", y = "", fill = "Location", color = "Location")

Chaq_chart
ggsave("analyses/plots/Chaq_histogram.pdf", units = "in", width = 14, height = 8)

# sequecing qPCR data
galbut_complete <- ggplot(filter(metadata, present == "yes")) +
  geom_jitter(aes(x = complete, y = ct_1600_1601, colour = low_pos, 
                  shape = low_pos), size = 3, alpha = 0.75, width = 0.2, 
              height = 0) +
  scale_colour_manual(values = c("darkblue", "orangered3")) +
  theme_bw() +
  facet_wrap(~location) +
  labs(x = "Complete galbut virus sequence",
       y = "Ct using 1600/1601 primers (detect galbut virus clade A)", 
       colour = "RT-qPCR low positve", 
       shape = "RT-qPCR low positve")

galbut_complete
ggsave("analyses/plots/galbut_complete_1600_1601.pdf", units = "in", width = 14, height = 8)

galbut_present_2 <- ggplot(metadata) +
  geom_jitter(aes(x = present, y = ct_2165_2170, colour = low_pos, 
                  shape = low_pos), size = 3, alpha = 0.75, width = 0.2, 
              height = 0) +
  scale_colour_manual(values = c("darkslateblue", "indianred3")) +
  theme_bw() +
  facet_wrap(~location) +
  labs(x = "Positive for galbut virus by sequencing",
       y = "Ct using 2165/2170 primers (detect galbut virus clade A or B)", 
       colour = "RT-qPCR low positve", 
       shape = "RT-qPCR low positve")

galbut_present_2
ggsave("analyses/plots/galbut_present_2165_2170.pdf", units = "in", width = 14, height = 8)

galbut_complete_2 <- ggplot(filter(metadata, present == "yes")) +
  geom_jitter(aes(x = complete, y = ct_2165_2170, colour = low_pos, 
                  shape = low_pos), size = 3, alpha = 0.75, width = 0.2, 
              height = 0) +
  scale_colour_manual(values = c("darkblue", "orangered3")) +
  theme_bw() +
  facet_wrap(~location) +
  labs(x = "Positive for galbut virus by sequencing",
       y = "Ct using 2165/2170 primers (detect galbut virus clade A or B)", 
       colour = "RT-qPCR low positve", 
       shape = "RT-qPCR low positve")

galbut_complete_2
ggsave("analyses/plots/galbut_complete_2165_2170.pdf", units = "in", width = 14, height = 8)


### Heatmaps

heatmap_co <- ggplot(filter(heatmap_data, location == "Colorado", 
                            segment %in% c("Chaq_total", "RNA3_total",
                                           "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  #  scale_fill_grey() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_co
ggsave("analyses/plots/co_all_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_pe <- ggplot(filter(heatmap_data, location == "Pennsylvania",
                            segment %in% c("Chaq_total", "RNA3_total",
                                           "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  #  scale_fill_grey() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_pe
ggsave("analyses/plots/pe_all_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_me <- ggplot(filter(heatmap_data, location == "Maine", 
                            segment %in% c("Chaq_total", "RNA3_total",
                                           "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  #  scale_fill_grey() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_me
ggsave("analyses/plots/me_all_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_un <- ggplot(filter(heatmap_data, location == "Unknown", 
                            segment %in% c("Chaq_total", "RNA3_total",
                                           "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  #  scale_fill_grey() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_un
ggsave("analyses/plots/unk_all_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_oh <- ggplot(filter(heatmap_data, location == "Ohio", 
                            segment %in% c("Chaq_total", "RNA3_total",
                                           "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  #  scale_fill_grey() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_oh
ggsave("analyses/plots/oh_all_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

ggplot(filter(qPCR, location %in% c("Adobe", "Briarwood", "CVID", "Elm", "James",
                                    "Linden", "Myrtle", "Wabash"))) +
  geom_jitter(aes(x = target, y = ct, color = positive, shape = positive), 
              alpha = 0.7, size = 3, width = 0.25) +
  scale_colour_manual(values = c("grey50", "blueviolet")) +
  facet_wrap(~factor(location, levels = c("Wabash", "Linden", "Elm", "CVID", 
                                          "Adobe", "Myrtle", "James", 
                                          "Briarwood")), nrow = 2) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "qPCR Target", y = "Ct", color = "True Positive?", shape = "True Positive?")
ggsave("local.pdf", units = "in", width = 10, height = 8)  


ggplot(filter(qPCR, location %in% c("Unknown", "Maine", "Pennsylvania", "Ohio"))) +
  geom_jitter(aes(x = target, y = ct, color = positive, shape = positive), 
              alpha = 0.7, size = 3, width = 0.25) +
  scale_colour_manual(values = c("grey50", "blueviolet")) +
  facet_wrap(~factor(location, levels = c("Pennsylvania", "Maine",
                                          "Ohio", "Unknown"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "qPCR Target", y = "Ct", color = "True Positive?", shape = "True Positive?")
ggsave("national.pdf", units = "in", width = 10, height = 8)  

#Normalizing

# Rpl for normalizing
Rpl <- qPCR %>% 
  filter(target == "RpL32",
         positive == "Yes") %>% 
  select(name, target, ct) 

qPCR_new <- qPCR %>% 
  filter(positive == "Yes",
         target != "RpL32")

qPCR_new <- left_join(qPCR_new, Rpl, by = join_by(name)) %>% 
  mutate(delta_ct = ct.x - ct.y)

local_normalized <- ggplot(filter(qPCR_new, location %in% c("Adobe", "Briarwood", 
                                                            "CVID", "Elm", "James", 
                                                            "Linden", "Myrtle", 
                                                            "Wabash"))) +
  geom_jitter(aes(x = target.x, y = delta_ct, color = target.x, shape = target.x), 
              alpha = 0.7, size = 3, width = 0.25) +
  scale_colour_manual(values = c("slateblue", "steelblue4")) +
  facet_wrap(~factor(location, levels = c("Wabash", "Linden", "Elm", "CVID", 
                                          "Adobe", "Myrtle", "James", 
                                          "Briarwood")), nrow = 2) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "RT-qPCR Target", y = "Delta Ct (Normalized to RpL32)", 
       color = "", shape = "")

local_normalized
ggsave("local_normalized.pdf", units = "in", width = 10, height = 8)  

national_normalized <- ggplot(filter(qPCR_new, location %in% c("Unknown", "Maine", 
                                                               "Pennsylvania", 
                                                               "Ohio"))) +
  geom_jitter(aes(x = target.x, y = delta_ct, color = target.x, shape = target.x), 
              alpha = 0.7, size = 3, width = 0.25) +
  scale_colour_manual(values = c("slateblue", "steelblue4")) +
  facet_wrap(~factor(location, levels = c("Pennsylvania", "Maine",
                                          "Ohio", "Unknown"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "RT-qPCR Target", y = "Delta Ct (Normalized to RpL32)", 
       color = "", shape = "")

national_normalized
ggsave("national_nromalized.pdf", units = "in", width = 10, height = 8)  

over_time <- ggplot(filter(qPCR_new, location %in% c("CVID", "Linden", "Wabash"))) +
  geom_jitter(aes(x = factor(date, levels = c("July", "August", "September", 
                                              "October")), y = delta_ct, 
                  color = target.x, shape = target.x), 
              alpha = 0.6, size = 2.75, width = 0.4) +
  scale_colour_manual(values = c("darkslategrey", "tomato4")) +
  facet_wrap(~factor(location), nrow = 3) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Month Collected", y = "Delta Ct (Normalized to RpL32)", 
       color = "RT-qPCR Target", shape = "RT-qPCR Target")

over_time
ggsave("over_time.pdf", units = "in", width = 10, height = 8)  

over_time_ct <- ggplot(filter(qPCR_new, location %in% c("CVID", "Linden", "Wabash"))) +
  geom_jitter(aes(x = factor(date, levels = c("July", "August", "September", 
                                              "October")), y = ct.x, 
                  color = target.x, shape = target.x), 
              alpha = 0.6, size = 2.75, width = 0.4) +
  scale_colour_manual(values = c("darkslategrey", "tomato4")) +
  facet_wrap(~factor(location), nrow = 3) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Month Collected", y = "Cycle Threshold", 
       color = "RT-qPCR Target", shape = "RT-qPCR Target")

over_time_ct
ggsave("over_time_ct.pdf", units = "in", width = 10, height = 8)  
genotype_ct <- ggplot(filter(qPCR_wide, location %in% c("Adobe", "Briarwood", 
                                                        "", "Elm", "James", 
                                                        "Linden", "Myrtle", 
                                                        "Wabash"))) +
  geom_jitter(aes(x = genotype, y = ct_use, color = genotype, shape = genotype), 
              alpha = 0.6, size = 3, width = 0.25) +
  scale_colour_manual(values = c("royalblue", "midnightblue")) +
  facet_wrap(~factor(location, levels = c("Wabash", "Linden", "Elm", "CVID", 
                                          "Adobe", "Myrtle", "James", 
                                          "Briarwood")), nrow = 2) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Galbut virus genotype", y = "Number of cycles to detection", 
       color = "Galbut virus\n genotype", shape = "Galbut virus\n genotype")

genotype_ct
ggsave("local_genotype.pdf", units = "in", width = 10, height = 8)  

genotype_ct_nat <- ggplot(filter(qPCR_wide, location %in% c("Unknown", "Maine", 
                                                            "Pennsylvania", 
                                                            "Ohio"))) +
  geom_jitter(aes(x = genotype, y = ct_use, color = genotype, shape = genotype), 
              alpha = 0.7, size = 3, width = 0.25) +
  scale_colour_manual(values = c("royalblue", "midnightblue")) +
  facet_wrap(~factor(location, levels = c("Pennsylvania", "Maine",
                                          "Ohio", "Unknown"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Galbut virus genotype", y = "Delta Ct (Normalized to RpL32)", 
       color = "Galbut virus\n genotype", shape = "Galbut virus\n genotype")

genotype_ct_nat
ggsave("national_genotype.pdf", units = "in", width = 10)

       