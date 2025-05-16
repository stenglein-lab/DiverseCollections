library(tidyverse)
library(FSA)
library(broom)

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
  filter(positive_galbut_A != "Yes" | positive_galbut_A_B != "No") %>% 
  mutate(ct_use = case_when(positive_galbut_A_B == "Yes" ~ ct_galbut_A_B,
                            positive_galbut_A_B == "No" ~ 0)) %>% 
  mutate(delta_ct = ct_use - ct_RpL32) %>% 
  mutate(relative = 2^-delta_ct) %>% 
  mutate(relative = case_when(positive_galbut_A_B == "Yes" ~ relative,
                              positive_galbut_A_B == "No" ~ 0))


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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
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


# STATS
# logistic glm of positive and location

qPCR_wide <- qPCR_wide %>% 
  mutate(pres_abs = case_when(positive_galbut_A_B == "Yes" ~ 1,
                              positive_galbut_A_B == "No" ~ 0))

loc_stat <- glm(pres_abs ~ location, data = qPCR_wide, family = "binomial")
summary(loc_stat)

#logistic glm of positive and time

qPCR_time <- qPCR_wide %>% 
  filter(location == "Wabash")

time_stat_w <- glm(pres_abs ~ date, data = qPCR_time, family = "binomial")
summary(time_stat_w)

qPCR_time <- qPCR_wide %>% 
  filter(location == "Rampart")

time_stat_r <- glm(pres_abs ~ date, data = qPCR_time, family = "binomial")
summary(time_stat_r)

descriptive <- qPCR_wide %>% 
  group_by(location, positive_galbut_A_B) %>% 
  summarise(n = n()) 

tot_loc <- qPCR_wide %>% 
  group_by(location) %>% 
  summarise(tot = n())

desc_loc <- left_join(descriptive, tot_loc, by = "location") %>% 
  mutate(prev = (n/tot)*100) %>% 
  filter(positive_galbut_A_B == "Yes")

prev_loc <- ggplot(desc_loc, aes(x = factor(location, levels = c("Adobe", "Briarwood", 
                                                                 "Elm", "James", 
                                                                 "JFK Parkway", 
                                                                 "Linden", "Myrtle", 
                                                                 "Rampart", "Wabash", "Maine", 
                                                                 "Ohio", "Pennsylvania")), 
                                 y = prev)) +
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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(x = "Sample Location", y = "Percent Positive", 
       fill = "Sample Location")

prev_loc
ggsave("analyses/plots/prev_loc.pdf", units = "in", width = 10, height = 8) 
  
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

prev_time <- ggplot(desc_time, aes(x = factor(date, levels = c("July", "August", "September", 
                                                               "October")), 
                                 y = prev)) +
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

