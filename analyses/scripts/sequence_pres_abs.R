# This script was used to look at galbut virus and chaq virus segement presence 
# in a heat map style figure. The code used for figure 4.
library(tidyverse)
library(readxl)
library(forcats)

metadata <- read_xlsx("analyses/data/Metadata_Table_OH_removed.xlsx")
all_rpm <- read_csv("analyses/data/all_rpm_counts.csv")

metadata <- left_join(metadata, all_rpm, by = "sample_name")

heatmap_data <- metadata %>% 
  select(sample_name, location, RNA1_submit, RNA2_submit, RNA3_submit, Chaq_submit,
         RNA1_total, RNA2_total, RNA3_total, Chaq_total, submit_co) %>% 
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
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("white", "#9ecae1", "#3182bd", "#08519c")) +
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
         RNA1_total, RNA2_total, RNA3_total, Chaq_total, submit_co) %>% 
  pivot_longer(cols = c("RNA1_submit", "RNA2_submit", "RNA3_submit", "Chaq_submit",
                        "RNA1_total", "RNA2_total", "RNA3_total", "Chaq_total"),
               names_to = "segment", 
               values_to = "number") %>% 
  mutate(segment = factor(segment, levels = c("Chaq_submit", "Chaq_total", 
                                              "RNA3_submit", "RNA3_total", 
                                              "RNA2_submit", "RNA2_total",
                                              "RNA1_submit", "RNA1_total")),
         location = factor(location)) %>% 
  group_by(location)



heatmap_co <- heatmap_pres %>% 
  filter(submit_co == "y") %>% 
  filter(location == "Colorado")

heatmap_me <- heatmap_pres %>% 
  filter(submit_co == "y") %>% 
  filter(location == "Maine")

heatmap_oh <- heatmap_pres %>% 
  filter(submit_co == "y") %>% 
  filter(location == "Ohio")

heatmap_pe <- heatmap_pres %>% 
  filter(submit_co == "y") %>% 
  filter(location == "Pennsylvania")

heatmap_sin_co <- heatmap_pres %>% 
  filter(submit_co == "n")  %>% 
  filter(location == "Colorado")

heatmap_sin_me <- heatmap_pres %>% 
  filter(submit_co == "n") %>% 
  filter(location == "Maine")

heatmap_sin_oh <- heatmap_pres %>% 
  filter(submit_co == "n") %>% 
  filter(location == "Ohio")

heatmap_sin_pe <- heatmap_pres %>% 
  filter(submit_co == "n") %>% 
  filter(location == "Pennsylvania")

heatmap_all_co <- ggplot(filter(heatmap_co, segment %in% c("Chaq_submit", "RNA3_submit",
                                                               "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("white", "#9ecae1", "#3182bd", "#08519c")) +
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

heatmap_all_co
ggsave("analyses/plots/heatmap_co_co_sub.pdf", units = "in", width = 14, height = 8)

heatmap_me_co <- ggplot(filter(heatmap_me, segment %in% c("Chaq_submit", "RNA3_submit",
                                                           "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
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

heatmap_me_co
ggsave("analyses/plots/heatmap_co_me_sub.pdf", units = "in", width = 14, height = 8)

heatmap_pe_co <- ggplot(filter(heatmap_pe, segment %in% c("Chaq_submit", "RNA3_submit",
                                                          "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
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

heatmap_pe_co
ggsave("analyses/plots/heatmap_co_pe_sub.pdf", units = "in", width = 14, height = 8)

heatmap_sin_co <- ggplot(filter(heatmap_sin_co, segment %in% c("Chaq_submit", "RNA3_submit",
                                                               "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("white", "#9ecae1", "#3182bd", "#08519c")) +
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

heatmap_sin_co
ggsave("analyses/plots/heatmap_sin_co_sub.pdf", units = "in", width = 14, height = 8)

heatmap_me_sin <- ggplot(filter(heatmap_sin_me, segment %in% c("Chaq_submit", "RNA3_submit",
                                                               "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("white", "#9ecae1", "#3182bd", "#08519c")) +
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

heatmap_me_sin
ggsave("analyses/plots/heatmap_sin_me_sub.pdf", units = "in", width = 14, height = 8)

heatmap_pe_sin <- ggplot(filter(heatmap_sin_pe, segment %in% c("Chaq_submit", "RNA3_submit",
                                                               "RNA2_submit", "RNA1_submit"))) +
  geom_tile(aes(x = sample_name, y = segment, fill = factor(number), height = 0.95)) +
  scale_fill_manual(values = c("white", "#9ecae1")) +
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

heatmap_pe_sin
ggsave("analyses/plots/heatmap_sin_pe_sub.pdf", units = "in", width = 14, height = 8)

heatmap_oh_sin <- ggplot(filter(heatmap_sin_oh, segment %in% c("Chaq_submit", "RNA3_submit",
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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6.5)) +
  labs(x = "Sample ID", y = "Galbut virus Segement", fill = "Number of \nsequences")

heatmap_oh_sin
ggsave("analyses/plots/heatmap_sin_oh_sub.pdf", units = "in", width = 14, height = 8)
