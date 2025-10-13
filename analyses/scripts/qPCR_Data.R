# This script was used to generate the qPCR relevant panels in figure 1. It reads in 
# qPCR data in a tidy format.
library(tidyverse)
library(FSA)
library(broom)
library(binom)
library(patchwork)

qPCR <- read_csv("analyses/data/qPCR_TidyData_LK_20240403.csv") %>% 
  mutate(target = str_replace(target, "galbut_1600_1601", "Galbut A"),
         target = str_replace(target, "galbut_2165_2170", "Galbut A or B"),
         positive = str_replace(positive, "N", "No"),
         positive = str_replace(positive, "Y", "Yes"),
         location = str_replace(location, "CVID", "Rampart"),
         location = str_replace(location, "Colorado", "JFK Parkway"))

descriptive <- qPCR %>% 
  group_by(target, positive) %>% 
  summarise(n = n())

# wide format

qPCR_wide <- pivot_wider(qPCR, names_from = target, values_from = c(ct, tm1, positive))

qPCR_wide <- qPCR_wide %>% 
  rename(positive_galbut_A = `positive_Galbut A`,
         positive_galbut_A_B = `positive_Galbut A or B`,
         ct_galbut_A = `ct_Galbut A`,
         ct_galbut_A_B = `ct_Galbut A or B`,
         tm_galbut_A = `tm1_Galbut A`,
         tm_galbut_A_B = `tm1_Galbut A or B`)  %>% 
  filter(positive_RpL32 == "Yes") %>% 
  mutate(ct_use = case_when(positive_galbut_A_B == "Yes" ~ ct_galbut_A_B,
                            positive_galbut_A_B == "No" ~ 0)) %>% 
  mutate(delta_ct = ct_use - ct_RpL32) %>% 
  mutate(relative = 2^-delta_ct) %>% 
  mutate(relative = case_when(positive_galbut_A_B == "Yes" ~ relative,
                              positive_galbut_A_B == "No" ~ 0)) %>% 
  mutate(dummy = "x")

descriptive <- qPCR_wide %>% 
  group_by(positive_galbut_A_B) %>% 
  summarise(n = n())

# All points
ct_all <- ggplot(qPCR_wide, aes(y = relative, x = "dummy")) +
  geom_violin() +
  geom_jitter(aes(color = factor(location, levels = c("Adobe", "Briarwood", 
                                                      "Elm", "James", 
                                                      "JFK Parkway", 
                                                      "Linden", "Myrtle", 
                                                      "Rampart", "Wabash", "Maine", 
                                                      "Ohio", "Pennsylvania"))), 
              alpha = 0.6, size = 2.5, width = 0.2, height = 0) +
  scale_colour_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                                 "darkorange4", "#00FFFF","#00CCCC","#006666")) +
  scale_y_continuous(trans="log10") +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        legend.position = "none",
        axis.text.x=element_blank(),
        panel.grid = element_blank()) +
  labs(y = "Normalized Galbut Virus Levels \nRelative to RpL32 mRNA", x = "")

ct_all
ggsave("analyses/plots/ct_all.pdf", units = "in", width = 4, height = 8)  

# Normalized RT-qPCR levels at all locations
ct_norm <- ggplot(qPCR_wide, aes(x = factor(location, levels = c("Adobe", "Briarwood", 
                                                                 "Elm", "James", 
                                                                 "JFK Parkway", 
                                                                 "Linden", "Myrtle", 
                                                                 "Rampart", "Wabash", "Maine", 
                                                                 "Ohio", "Pennsylvania")), y = relative)) +
  geom_violin() +
  geom_jitter(aes(color = factor(location, levels = c("Adobe", "Briarwood", 
                                                      "Elm", "James", 
                                                      "JFK Parkway", 
                                                      "Linden", "Myrtle", 
                                                      "Rampart", "Wabash", "Maine", 
                                                      "Ohio", "Pennsylvania"))), 
              alpha = 0.6, size = 2.5, width = 0.1) +
  scale_colour_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                                 "darkorange4", "#00FFFF","#00CCCC","#006666")) +
  scale_y_continuous(trans="log10") +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(x = "Sample Location", y = "Normalized Galbut Virus Levels \nRelative to RpL32 mRNA", 
       color = "Sample Location")

ct_norm
ggsave("ct_norm.pdf", units = "in", width = 10, height = 8)  


# Over time
over_time <- ggplot(filter(qPCR_wide, location %in% c("Rampart", "Wabash"))) +
  geom_jitter(aes(x = factor(date, levels = c("July", "August", "September", 
                                              "October")), y = relative, 
                  color = location), 
              alpha = 0.6, size = 2.75, width = 0.3) +
  scale_y_continuous(trans = "log10") +
  scale_colour_manual(values = c("darkorange2", "darkorange4")) +
  facet_wrap(~factor(location), nrow = 2) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Month Collected", y = "Normalized Galbut Virus Levels \nRelative to RpL32 mRNA", 
       color = "Sample Location")

over_time
ggsave("over_time.pdf", units = "in", width = 10, height = 8)  

descriptive <- qPCR_wide %>% 
  group_by(location, positive_galbut_A_B) %>% 
  summarise(n = n()) 

tot_loc <- qPCR_wide %>% 
  group_by(location) %>% 
  summarise(tot = n())

desc_loc <- left_join(descriptive, tot_loc, by = "location") %>% 
  mutate(prev = (n/tot)*100) %>% 
  filter(positive_galbut_A_B == "Yes")

# calculate binomial confidence intervals
prev_loc_conf_int <- binom.confint(desc_loc$n, 
                               desc_loc$tot, 
                               method="exact")

prev_loc_conf_int <- prev_loc_conf_int %>% 
  mutate(low_pre = (lower * 100),
         high_pre = (upper *100))

#All prevalance & binomial confidence int
all_tot <- sum(desc_loc$tot)
all_pos <- sum(desc_loc$n)

all_prev <- (all_pos/all_tot) * 100

loc <- 1

all_pre <- list(total = all_tot, positive = all_pos, prevalance = all_prev, 
                dummy = loc)
all_pre <- as.data.frame(all_pre)

all_prev_conf_int <- binom.confint(all_pre$positive, all_pre$total, 
                                   method = "exact")

all_prev_conf_int <- all_prev_conf_int %>% 
  mutate(low_pre = (lower *100),
         high_pre = (upper * 100))

prev_loc <- ggplot(desc_loc, aes(x = factor(location, levels = c("Adobe", "Briarwood", 
                                                                 "Elm", "James", 
                                                                 "JFK Parkway", 
                                                                 "Linden", "Myrtle", 
                                                                 "Rampart", "Wabash", "Maine", 
                                                                 "Ohio", "Pennsylvania")), 
                                 y = prev)) +
  geom_errorbar(aes(ymin = prev_loc_conf_int$low_pre, ymax = prev_loc_conf_int$high_pre), width = 0.1) +
  geom_point(aes(fill = factor(location, levels = c("Adobe", "Briarwood", 
                                                    "Elm", "James", 
                                                    "JFK Parkway", 
                                                    "Linden", "Myrtle", 
                                                    "Rampart", "Wabash", "Maine", 
                                                    "Ohio", "Pennsylvania"))),
             shape = 21, size = 6) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                                 "darkorange4", "#00FFFF","#00CCCC","#006666")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(x = "Sample Location", y = "Percent Positive", 
       fill = "Sample Location")

prev_loc
ggsave("analyses/plots/prev_loc.pdf", units = "in", width = 10, height = 8) 

prev_all <- ggplot(all_pre, aes(x = dummy, y = prevalance)) +
  geom_errorbar(aes(ymin = all_prev_conf_int$low_pre, ymax = all_prev_conf_int$high_pre), width = 0.1) +
  geom_point(shape = 21, size = 5) +
  ylim(0,100) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(y = "Percent Positive", x = "All Locations")

prev_all
ggsave("analyses/plots/prev_all.pdf", units = "in", width = 4, height = 8) 
  
desc_time <- qPCR_wide %>% 
  filter(location %in% c("Rampart", "Wabash")) %>% 
  group_by(location, positive_galbut_A_B, date) %>% 
  summarise(n = n()) 

     
tot_time <- qPCR_wide %>% 
  filter(location %in% c("Rampart", "Wabash")) %>% 
  group_by(location, date) %>% 
  summarise(tot = n())

desc_time <- left_join(desc_time, tot_time, by = c("location", "date")) %>% 
  mutate(prev = (n/tot)*100) %>% 
  filter(positive_galbut_A_B == "Yes")

# calculate binomial confidence intervals
prev_time_conf_int <- binom.confint(desc_time$n, 
                                   desc_time$tot, 
                                   method="exact")

prev_time_conf_int <- prev_time_conf_int %>% 
  mutate(low_pre = (lower * 100),
         high_pre = (upper *100))

prev_time <- ggplot(desc_time, aes(x = factor(date, levels = c("July", "August", "September", 
                                                               "October")), 
                                 y = prev)) +
  geom_errorbar(aes(ymin = prev_time_conf_int$low_pre, ymax = prev_time_conf_int$high_pre), width = 0.1) +
  geom_point(aes(fill = factor(location, levels = c("Rampart", "Wabash"))),
             shape = 21, size = 6) +
  scale_fill_manual(values = c("darkorange2", "darkorange4")) +
  ylim(0, 100) +
  facet_wrap(~location, nrow = 2) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(x = "Sample Location", y = "Percent Positive", 
       fill = "Sample Location")

prev_time
ggsave("analyses/plots/prev_time.pdf", units = "in", width = 10, height = 8) 

# Non parametric anova (Kruskal-Wallis test)- is there a difference between location
diff_loc <- qPCR_wide %>% 
  select(location, relative)

kruskal.test(relative ~ location, data = diff_loc)

diff_comp <- dunnTest(relative ~ location, data = diff_loc, method = "bonferroni", list = TRUE)
loc_comp_df <- diff_comp[['res']]

# Is there a difference over time?
diff_time_ramp <- qPCR_wide %>% 
  filter(location == "Rampart") %>% 
  select(relative, date, location)

kruskal.test(relative ~ date, data = diff_time_ramp)
diff_comp_ram <- dunnTest(relative ~ date, data = diff_time_ramp, method = "bonferroni", list = TRUE)
loc_comp_ram <- diff_comp_ram[['res']]

diff_time_wab <- qPCR_wide %>% 
  filter(location == "Wabash") %>% 
  select(relative, date, location)

kruskal.test(relative ~ date, data = diff_time_wab)
diff_comp_wab <- dunnTest(relative ~ date, data = diff_time_wab, method = "bonferroni", list = TRUE)
loc_comp_wab <- diff_comp_wab[['res']]



