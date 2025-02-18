library(tidyverse)
library(readxl)
library(viridis)

metadata <- read_xlsx("analyses/data/Metadata_Table.xlsx")

heatmap_data <- metadata %>% 
  select(sample_name, location, RNA1_submit, RNA2_submit, RNA3_submit, Chaq_submit,
         RNA1_total, RNA2_total, RNA3_total, Chaq_total) %>% 
  pivot_longer(cols = c("RNA1_submit", "RNA2_submit", "RNA3_submit", "Chaq_submit",
                        "RNA1_total", "RNA2_total", "RNA3_total", "Chaq_total"),
               names_to = "segment", 
               values_to = "number") %>% 
  mutate(segment = factor(segment, levels = c("Chaq_submit", "Chaq_total", 
                                              "RNA3_submit", "RNA3_total", 
                                              "RNA2_submit", "RNA2_total",
                                              "RNA1_submit", "RNA1_total")))

heatmap_all <- ggplot(filter(heatmap_data, segment %in% c("Chaq_total", "RNA3_total",
                                                          "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd", "#08519c")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_all
ggsave("analyses/plots/allll_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_pres <- metadata %>% 
  filter(num_seqs_total > 0) %>% 
  select(sample_name, location, RNA1_submit, RNA2_submit, RNA3_submit, Chaq_submit,
         RNA1_total, RNA2_total, RNA3_total, Chaq_total) %>% 
  pivot_longer(cols = c("RNA1_submit", "RNA2_submit", "RNA3_submit", "Chaq_submit",
                        "RNA1_total", "RNA2_total", "RNA3_total", "Chaq_total"),
               names_to = "segment", 
               values_to = "number") %>% 
  mutate(segment = factor(segment, levels = c("Chaq_submit", "Chaq_total", 
                                              "RNA3_submit", "RNA3_total", 
                                              "RNA2_submit", "RNA2_total",
                                              "RNA1_submit", "RNA1_total")))

heatmap_all_pres <- ggplot(filter(heatmap_pres, segment %in% c("Chaq_total", "RNA3_total",
                                                               "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd", "#08519c")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6.5)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_all_pres
ggsave("analyses/plots/all_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_co2 <- ggplot(filter(heatmap_pres, location == "Colorado", 
                             segment %in% c("Chaq_total", "RNA3_total",
                                            "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd")) +
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

heatmap_co2
ggsave("analyses/plots/co_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_pe2 <- ggplot(filter(heatmap_pres, location == "Pennsylvania",
                             segment %in% c("Chaq_total", "RNA3_total",
                                            "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
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

heatmap_pe2
ggsave("analyses/plots/pe_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_me2 <- ggplot(filter(heatmap_pres, location == "Maine", 
                             segment %in% c("Chaq_total", "RNA3_total",
                                            "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd", "#08519c")) +
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

heatmap_me2
ggsave("analyses/plots/me_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_un2 <- ggplot(filter(heatmap_pres, location == "Unknown", 
                             segment %in% c("Chaq_total", "RNA3_total",
                                            "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd")) +
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

heatmap_un2
ggsave("analyses/plots/un_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_oh2 <- ggplot(filter(heatmap_pres, location == "Ohio", 
                             segment %in% c("Chaq_total", "RNA3_total",
                                            "RNA2_total", "RNA1_total"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.85)) +
  scale_fill_manual(values = c("#9ecae1")) +
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

heatmap_oh2
ggsave("analyses/plots/oh_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)


# submitted
heatmap_sub <- metadata %>% 
  filter(num_seqs_submit > 0) %>% 
  select(sample_name, location, RNA1_submit, RNA2_submit, RNA3_submit, Chaq_submit,
         RNA1_total, RNA2_total, RNA3_total, Chaq_total) %>% 
  pivot_longer(cols = c("RNA1_submit", "RNA2_submit", "RNA3_submit", "Chaq_submit",
                        "RNA1_total", "RNA2_total", "RNA3_total", "Chaq_total"),
               names_to = "segment", 
               values_to = "number") %>% 
  mutate(segment = factor(segment, levels = c("Chaq_submit", "Chaq_total", 
                                              "RNA3_submit", "RNA3_total", 
                                              "RNA2_submit", "RNA2_total",
                                              "RNA1_submit", "RNA1_total")))

heatmap_all_sub <- ggplot(filter(heatmap_sub, segment %in% c("Chaq_submit", "RNA3_submit",
                                                               "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6.5)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_all_sub
ggsave("analyses/plots/all_sub_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_co3 <- ggplot(filter(heatmap_sub, location == "Colorado", 
                             segment %in% c("Chaq_submit", "RNA3_submit",
                                            "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd")) +
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

heatmap_co3
ggsave("analyses/plots/co_sub_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_pe3 <- ggplot(filter(heatmap_sub, location == "Pennsylvania",
                             segment %in% c("Chaq_submit", "RNA3_submit",
                                            "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd")) +
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

heatmap_pe3
ggsave("analyses/plots/pe_sub_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_me3 <- ggplot(filter(heatmap_sub, location == "Maine", 
                             segment %in% c("Chaq_submit", "RNA3_submit",
                                            "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd")) +
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

heatmap_me3
ggsave("analyses/plots/me_sub_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_un3 <- ggplot(filter(heatmap_sub, location == "Unknown", 
                             segment %in% c("Chaq_submit", "RNA3_submit",
                                            "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#deebf7", "#9ecae1", "#3182bd")) +
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

heatmap_un3
ggsave("analyses/plots/un_sub_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)

heatmap_oh3 <- ggplot(filter(heatmap_sub, location == "Ohio", 
                             segment %in% c("Chaq_submit", "RNA3_submit",
                                            "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#9ecae1")) +
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

heatmap_oh3
ggsave("analyses/plots/oh_sub_heatmap_virus_pres.pdf", units = "in", width = 14, height = 8)
