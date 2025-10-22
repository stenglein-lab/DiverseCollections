library(tidyverse)
library(rstatix)
library(binom)
library(patchwork)
library(lme4)
library(multimode)
library(scales)

qPCR <- read_csv("analyses/data/qPCR_TidyData_LK_20240403.csv") %>% 
  mutate(target = str_replace(target, "galbut_1600_1601", "Galbut A"),
         target = str_replace(target, "galbut_2165_2170", "Galbut A or B"),
         positive = str_replace(positive, "N", "No"),
         positive = str_replace(positive, "Y", "Yes"),
         location = str_replace(location, "CVID", "Rampart"),
         location = str_replace(location, "Colorado", "JFK Parkway"))


# how many samples in total?
sample_names   <- qPCR %>% group_by(name) %>% summarize()
location_names <- tibble(location = str_match(sample_names$name, "([^_]+)_")[,2]) %>% group_by(location) %>% summarize()
n_samples      <- nrow(sample_names)

# wide format
qPCR_wide <- pivot_wider(qPCR, names_from = target, values_from = c(ct, tm1, positive))


# filter out Ohio samples without a valid RpL32 result: these are not melanogaster
# and we are excluding them from the paper

Ohio_samples <- filter(qPCR_wide, location == "Ohio")
Ohio_samples_to_exclude <- filter(Ohio_samples, positive_RpL32 == "No") %>% pull(name)
qPCR_wide <- filter(qPCR_wide, !name %in% Ohio_samples_to_exclude)

# total samples included in study
samples_in_study                  <- nrow(qPCR_wide)
samples_in_study_with_valid_rpl32 <- nrow(filter(qPCR_wide, positive_RpL32 == "Yes"))

# calculate delta Cts and relative RNA levels
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


# how many samples galbut positive
samples_galbut_positive <- nrow(filter(qPCR_wide, positive_galbut_A_B == "Yes"))
percent_galbut_positive <- samples_galbut_positive / samples_in_study_with_valid_rpl32

minimum_relative_galbut <- min(filter(qPCR_wide, relative > 0) %>% pull(relative))
maximum_relative_galbut <- max(filter(qPCR_wide, relative > 0) %>% pull(relative))

# calculate modes of distribution of galbut virus relative RNA levels
# see: https://stackoverflow.com/questions/27418461/calculate-the-modes-in-a-multimodal-distribution-in-r
mm <- log10(qPCR_wide %>% filter(is.finite(log10(relative))) %>% pull(relative))

# Select an appropriate bandwidth. Here we use a "rule of thumb" method
bandwidth <- bw.nrd0(mm) 

# Methods implemented in multimode are documented at ?bw.nrd
# Calculate number of modes
n <- nmodes(mm, bw = bandwidth)

# Calculate location of these modes
loc <- locmodes(mm, mod0 = n, display = TRUE)
# Remove `, display = TRUE` if you don't want to see the KDE
loc

# what is the relative RNA value at the 2 main nodes
low_positive_relative_exponent  <- loc$locations[1]
high_positive_relative_exponent <- loc$locations[3]

low_positive_relative  <- 10^abs(loc$locations[1])
high_positive_relative <- 10^loc$locations[3]

# output text in paper re: modes of distribution
print (paste0(
# "We identified the 2 modes of this distribution. ",
"The higher mode to correspond to ",
sprintf ("%0.1f ", high_positive_relative),
"galbut virus RNA copies per RpL32 mRNA, and the lower mode to ",
"1 galbut virus RNA per ",
sprintf ("%0.0f ", low_positive_relative),
"RpL32 mRNA copies."))

# Relative galbut virus RNA levels for all samples
ct_all <- ggplot(qPCR_wide, aes(y = relative, x = "All locations")) +
  geom_violin(linewidth=0.5, width=1) +
  geom_jitter(aes(fill = factor(location, levels = c("Adobe", "Briarwood", 
                                                      "Elm", "James", 
                                                      "JFK Parkway", 
                                                      "Linden", "Myrtle", 
                                                      "Rampart", "Wabash", "Maine", 
                                                      "Ohio", "Pennsylvania"))), 
              shape=21, alpha = 0.6, stroke=0.25, size = 1, width = 0.2, height = 0) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                                 "darkorange4", "#00FFFF","#00CCCC","#006666")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.title = element_text(size = 9),
        ) +
  labs(y = "Galbut virus RNA levels\nrelative to RpL32 mRNA",
       x = "")

ct_all

# Normalized RT-qPCR levels at all locations
ct_norm <- ggplot(qPCR_wide, aes(x = factor(location, levels = c("Adobe", "Briarwood", 
                                                                 "Elm", "James", 
                                                                 "JFK Parkway", 
                                                                 "Linden", "Myrtle", 
                                                                 "Rampart", "Wabash", "Maine", 
                                                                 "Ohio", "Pennsylvania")), y = relative)) +
  geom_violin(linewidth=0.25) +
  geom_jitter(aes(fill = factor(location, levels = c("Adobe", "Briarwood", 
                                                      "Elm", "James", 
                                                      "JFK Parkway", 
                                                      "Linden", "Myrtle", 
                                                      "Rampart", "Wabash", "Maine", 
                                                      "Ohio", "Pennsylvania"))), 
              shape=21, stroke=0.25, alpha = 0.6, size = 1, width = 0.2) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                                 "darkorange4", "#00FFFF","#00CCCC","#006666")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title = element_text(size=10),
        legend.position = "none") +
  labs(x = "", 
       y = "Galbut virus RNA levels\nrelative to RpL32 mRNA")
       
ct_norm


# Galbut virus RNA levels over time at locations sampled serially 

over_time <- ggplot(filter(qPCR_wide, location %in% c("Rampart", "Wabash"))) +
  geom_jitter(aes(x = factor(date, levels = c("July", "August", "September", 
                                              "October")), y = relative, 
                  fill = location), 
              shape=21, stroke=0.25, alpha = 0.6, size = 1, width = 0.3) +
  scale_fill_manual(values = c("darkorange2", "darkorange4")) +
  facet_wrap(~factor(location), nrow = 2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_rect(color="black", fill=NA, linewidth=0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        strip.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        ) +
  labs(x = "", y = "Galbut virus RNA levels\nrelative to RpL32 mRNA") 
       
over_time

# calculate prevalences

# qPCR_wide contains only valid samples (with a positive RpL32)
pos_neg_per_location <- qPCR_wide %>% 
  group_by(location, positive_galbut_A_B) %>% 
  summarise(n = n()) 

tot_loc <- qPCR_wide %>% 
  group_by(location) %>% 
  summarise(tot = n())

positives_per_location <- left_join(pos_neg_per_location, tot_loc, by = "location") %>% 
  mutate(prev = (n/tot)*100) %>% 
  filter(positive_galbut_A_B == "Yes")

# calculate binomial confidence intervals
prev_loc_conf_int <- binom.confint(positives_per_location$n, 
                               positives_per_location$tot, 
                               method="exact")

# prevalence with confidence intervals
prev_and_conf_int <- cbind(positives_per_location, select(prev_loc_conf_int, mean, lower, upper)) 

prev_loc_conf_int <- prev_loc_conf_int %>% 
  mutate(low_pre = (lower * 100),
         high_pre = (upper *100))

#All prevalance & binomial confidence int
all_tot <- sum(positives_per_location$tot)
all_pos <- sum(positives_per_location$n)

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

# prevalence plot for each location 
prev_loc <- ggplot(positives_per_location, aes(x = factor(location, levels = c("Adobe", "Briarwood", 
                                                                 "Elm", "James", 
                                                                 "JFK Parkway", 
                                                                 "Linden", "Myrtle", 
                                                                 "Rampart", "Wabash", "Maine", 
                                                                 "Ohio", "Pennsylvania")), 
                                 y = prev)) +
  geom_errorbar(aes(ymin = prev_loc_conf_int$low_pre, ymax = prev_loc_conf_int$high_pre), width = 0.25, linewidth=0.25) +
  geom_point(aes(fill = factor(location, levels = c("Adobe", "Briarwood", 
                                                    "Elm", "James", 
                                                    "JFK Parkway", 
                                                    "Linden", "Myrtle", 
                                                    "Rampart", "Wabash", "Maine", 
                                                    "Ohio", "Pennsylvania"))),
             shape = 21, size = 2, stroke=0.25) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                                 "darkorange4", "#00FFFF","#00CCCC","#006666")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title = element_text(size = 9),
        legend.position = "none") +
  labs(x = "", 
       y = "Galbut virus prevalence")
       

prev_loc

# prevalence for all samples combined
prev_all <- ggplot(all_pre, aes(x = "All locations", y = prevalance)) +
  geom_errorbar(aes(ymin = all_prev_conf_int$low_pre, ymax = all_prev_conf_int$high_pre), width = 0.5, linewidth=0.25) +
  geom_point(fill = "black", shape = 21, size = 2, stroke=0.25) +
  ylim(0,100) +
  theme_bw(base_size = 11) +
  theme(axis.title = element_text(size = 9),
        legend.position = "none") +
  labs(y = "Galbut virus prevalence", x = "")

prev_all
  
prevalence_by_time <- qPCR_wide %>% 
  filter(location %in% c("Rampart", "Wabash")) %>% 
  group_by(location, positive_galbut_A_B, date) %>% 
  summarise(n = n()) 
     
tot_time <- qPCR_wide %>% 
  filter(location %in% c("Rampart", "Wabash")) %>% 
  group_by(location, date) %>% 
  summarise(tot = n())

prevalence_by_time <- left_join(prevalence_by_time, tot_time, by = c("location", "date")) %>% 
  mutate(prev = (n/tot)*100) %>% 
  filter(positive_galbut_A_B == "Yes")

# calculate binomial confidence intervals
prev_time_conf_int <- binom.confint(prevalence_by_time$n, 
                                   prevalence_by_time$tot, 
                                   method="exact")

prev_time_conf_int <- prev_time_conf_int %>% 
  mutate(low_pre = (lower * 100),
         high_pre = (upper *100))

# calculate prevalence for samples that were sampled serially
prev_time <- ggplot(prevalence_by_time, 
                    aes(x = factor(date, 
                                   levels = c("July", "August", "September", "October")), 
                        y = prev)) +
  geom_errorbar(aes(ymin = prev_time_conf_int$low_pre, ymax = prev_time_conf_int$high_pre), width = 0.25, linewidth=0.25) +
  geom_point(aes(fill = factor(location, levels = c("Rampart", "Wabash"))),
             shape = 21, size = 2, stroke=0.25) +
  scale_fill_manual(values = c("darkorange2", "darkorange4")) +
  ylim(0, 100) +
  facet_wrap(~location, nrow = 2) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9),
        strip.background = element_rect(color="black", fill=NA, linewidth=0.5),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
  ) +
  labs(x = "", y = "Galbut virus prevalence") 
       

prev_time


# use patchwork to combine all panels
# this will label panels A->F
# we will relabel panels in Affinity Designer (B->...)
ct_norm + ct_all + over_time + prev_loc + prev_all + prev_time + 
  plot_layout(ncol = 3, widths = c(5,1,4)) +
  plot_annotation(tag_levels = "A")

# save PDF of these figure panels 
# next will open in Affinity Designer and combine with map images to make figure 1
ggsave("analyses/plots/Fig_1_B-E.pdf", units="in", width=8.5, height=7)


# ------------------------------------------------------------------------------
# Statistcal testing - does prevalence vary as a function of location and time?
# ------------------------------------------------------------------------------

# GLM to test for differences in prevalence by: 
# location
# location-date
# sex 

glm_input <- 
  qPCR_wide %>% 
  filter(positive_RpL32 == "Yes") %>%
  select(name, location, sex, date, positive_galbut_A_B, relative) %>%
  mutate(infected_boolean = if_else(positive_galbut_A_B == "Yes", T, F)) %>%
  select(-positive_galbut_A_B)

# make Wabash the first level so it will be reference location in glm
glm_input$location <- as.factor(glm_input$location)
glm_input$location <- fct_relevel(glm_input$location, "Wabash")

glm_input %>% group_by(location, date) %>% summarize(n=n())


# model prevalence by location
glm_prev_by_loc <-  glm(infected_boolean ~ 
                          location,  
                        data = glm_input, 
                        family = binomial) 

summary(glm_prev_by_loc)

# model prevalence by date-location
glm_input_loc_with_multiple_dates <- 
  glm_input %>% 
  filter(location %in% c("Wabash", "Rampart"))

glm_prev_by_loc_date <-  glm(infected_boolean ~ 
                               location + date + location:date,  
                             data = glm_input_loc_with_multiple_dates, 
                             family = binomial) 
summary(glm_prev_by_loc_date)

# run an all locations-vs-all locations pairwise Fisher's exact test
location_contingency_table <- table(glm_input$location, glm_input$infected_boolean)

pairwise_fet <- 
  pairwise_fisher_test(
    location_contingency_table, 
    simulate.p.value = T, 
    p.adjust.method = "holm")

pairwise_fet
write.table(pairwise_fet, 
            file="figures/tables/supplemental_table_2_pairwise_fet_by_location.tsv", 
            sep="\t",
            row.names=F,
            quote=F)

# do Fisher's exact tests for locations with multiple timepoints

glm_input_wabash  <- glm_input %>%  filter(location == "Wabash")
glm_input_linden  <- glm_input %>%  filter(location == "Linden")
glm_input_rampart <- glm_input %>%  filter(location == "Rampart")

wabash_contingency_table  <- table(glm_input_wabash$date,  glm_input_wabash$infected_boolean)
rampart_contingency_table <- table(glm_input_rampart$date, glm_input_rampart$infected_boolean)
linden_contingency_table  <- table(glm_input_linden$date,  glm_input_linden$infected_boolean)

wabash_pairwise_fet <- 
  pairwise_fisher_test(
    wabash_contingency_table, 
    simulate.p.value = T, 
    p.adjust.method = "holm")

rampart_pairwise_fet <- 
  pairwise_fisher_test(
    rampart_contingency_table, 
    simulate.p.value = T, 
    p.adjust.method = "holm")

linden_pairwise_fet <- 
  pairwise_fisher_test(
    linden_contingency_table, 
    simulate.p.value = T, 
    p.adjust.method = "holm")

# combine Wabash and Rampart results
wabash_pairwise_fet$location  <- "Wabash"
rampart_pairwise_fet$location <- "Rampart"
location_fet <- rbind(wabash_pairwise_fet, rampart_pairwise_fet) %>% relocate(location)

write.table(location_fet, 
            file="figures/tables/supplemental_table_3_pairwise_fet_by_time.tsv", 
            sep="\t",
            row.names=F,
            quote=F)

