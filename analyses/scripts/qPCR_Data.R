library(tidyverse)

qPCR <- read_csv("analyses/data/qPCR_TidyData_LK_20240403.csv") %>% 
  mutate(target = str_replace(target, "galbut_1600_1601", "Galbut A"),
         target = str_replace(target, "galbut_2165_2170", "Galbut A or B"),
         positive = str_replace(positive, "N", "No"),
         positive = str_replace(positive, "Y", "Yes"))

descriptive <- qPCR %>% 
  group_by(target, positive) %>% 
  summarise(n = n())

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





# wide format for genotype

qPCR_wide <- pivot_wider(qPCR, names_from = target, values_from = c(ct, tm1, positive))

qPCR_wide <- qPCR_wide %>% 
  rename(positive_galbut_A = `positive_Galbut A`,
         positive_galbut_A_B = `positive_Galbut A or B`,
         ct_galbut_A = `ct_Galbut A`,
         ct_galbut_A_B = `ct_Galbut A or B`,
         tm_galbut_A = `tm1_Galbut A`,
         tm_galbut_A_B = `tm1_Galbut A or B`)  %>% 
  filter(positive_RpL32 == "Yes") %>% 
  filter(positive_galbut_A == "Yes" | positive_galbut_A_B == "Yes") %>% 
  mutate(genotype = case_when(positive_galbut_A == "Yes" & positive_galbut_A_B == "Yes" ~ "A",
                              positive_galbut_A == "Yes" & positive_galbut_A_B == "No" ~ "A",
                              positive_galbut_A == "No" & positive_galbut_A_B == "Yes" ~ "B")) %>% 
  mutate(ct_use = ifelse(genotype == "A", ct_galbut_A, ct_galbut_A_B))

genotype_ct <- ggplot(filter(qPCR_wide, location %in% c("Adobe", "Briarwood", 
                                                 "CVID", "Elm", "James", 
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
ggsave("national_genotype.pdf", units = "in", width = 10, height = 8) 

# Normalized

qPCR_wide <- qPCR_wide %>% 
  mutate(delta_ct = ct_use - ct_RpL32)

genotype_ct_norm <- ggplot(filter(qPCR_wide, location %in% c("Adobe", "Briarwood", 
                                                        "CVID", "Elm", "James", 
                                                        "Linden", "Myrtle", 
                                                        "Wabash"))) +
  geom_jitter(aes(x = genotype, y = delta_ct, color = genotype, shape = genotype), 
              alpha = 0.6, size = 3, width = 0.25) +
  scale_colour_manual(values = c("brown3", "deeppink4")) +
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

genotype_ct_norm
ggsave("local_genotype_norm.pdf", units = "in", width = 10, height = 8)  

genotype_ct_nat_norm <- ggplot(filter(qPCR_wide, location %in% c("Unknown", "Maine", 
                                                            "Pennsylvania", 
                                                            "Ohio"))) +
  geom_jitter(aes(x = genotype, y = delta_ct, color = genotype, shape = genotype), 
              alpha = 0.7, size = 3, width = 0.25) +
  scale_colour_manual(values = c("brown3", "deeppink4")) +
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

genotype_ct_nat_norm
ggsave("national_genotype_norm.pdf", units = "in", width = 10, height = 8) 

# Genotype Over Time
over_time <- ggplot(filter(qPCR_wide, location %in% c("CVID", "Linden", "Wabash"))) +
  geom_jitter(aes(x = factor(date, levels = c("July", "August", "September", 
                                              "October")), y = delta_ct, 
                  color = genotype, shape = genotype), 
              alpha = 0.6, size = 2.75, width = 0.3) +
  scale_colour_manual(values = c("darkslategray", "darkturquoise")) +
  facet_wrap(~factor(location), nrow = 3) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Month Collected", y = "Delta Ct (Normalized to RpL32)", 
       color = "Galbut virus\n genotype", shape = "Galbut virus\n genotype")

over_time
ggsave("over_time_genotype.pdf", units = "in", width = 10, height = 8)  

over_time_ct <- ggplot(filter(qPCR_wide, location %in% c("CVID", "Linden", "Wabash"))) +
  geom_jitter(aes(x = factor(date, levels = c("July", "August", "September", 
                                              "October")), y = ct_use, 
                  color = genotype, shape = genotype), 
              alpha = 0.6, size = 2.75, width = 0.3) +
  scale_colour_manual(values = c("darkslategray", "darkturquoise")) +
  facet_wrap(~factor(location), nrow = 3) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Month Collected", y = "Cycle Threshold", 
       color = "Galbut virus\n genotype", shape = "Galbut virus\n genotype")

over_time_ct
ggsave("over_time_ct_genotype.pdf", units = "in", width = 10, height = 8) 


  