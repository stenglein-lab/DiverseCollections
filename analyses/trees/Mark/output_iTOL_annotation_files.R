library(tidyverse)

# This script reads in a file of sample IDS
# that correspond to tip labels in a tree
# and does a couple things:
#
# 1) It parses sequence IDs and in some cases infers shortened names
# 2) It outputs iTOL annotation files that can be used to color tip labels
#    in iTOL trees and relabel tips in iTOL trees.
#
#  It selects colors based on accession or sequence IDs.
#
# For more information on iTOL annotation files, see:
# https://itol.embl.de/help.cgi#annoTemplate


# read in sequence IDS
sample_ids <- read.delim("rna1_sequence_ids.txt", header=F)
colnames(sample_ids) <- c("sequence_id")

# try to parse out sequence IDs for new DC sequences
metadata <- str_match(sample_ids$sequence_id, "D_(sim|mel)_(\\S+)_(\\d{4})_(.*)_(Galbut|Chaq)")

# try to parse out RNA segment from sequence ID
rna_segment <- str_match(sample_ids$sequence_id, "(RNA\\d{1})")

sample_ids$parsed_spp <- metadata[,2]
sample_ids$parsed_loc <- metadata[,3]
sample_ids$parsed_date <- metadata[,4]
sample_ids$parsed_id <- metadata[,5]
sample_ids$parsed_virus <- metadata[,6]
sample_ids$parsed_segment <- rna_segment[,2]
sample_ids <- 
  sample_ids %>% mutate (rna_segment = case_when(
    parsed_virus == "Chaq" ~ "Chaq",
    !is.na(parsed_segment) ~ parsed_segment,
    .default = NA_character_
  ))

# use parsed ID or just the sequence ID (=accession)
sample_ids <- sample_ids %>% mutate(sample_id = if_else(!is.na(parsed_id), parsed_id, sequence_id))

# write out an iTOL annotation file that will relabel tree tips
sample_map <- sample_ids %>% select(sequence_id, sample_id)

# write out header lines to make this an iTOL annotation file
# in the style of labels_template.txt
write_itol_labels_header <- function(filename){
  fileConn<-file(filename)
  writeLines(c("LABELS","SEPARATOR SPACE","DATA"), fileConn)
  close(fileConn)
}
filename <- "galbut_RNA1_label_format.txt"
write_itol_labels_header(filename)
write.table(sample_map, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

# decide on label and branch colors
# based on sequence ID
#

# what are sample types (categories)
sample_ids <- sample_ids %>% mutate(
  sample_type = case_when(
    str_detect(sample_id, "ME")    ~ "Maine",
    str_detect(sample_id, "Penn")  ~ "Pennsylvania",
    !is.na(parsed_id)              ~ "Colorado_new",
    str_detect(sample_id, "MT742") ~ "Colorado",
    str_detect(sample_id, "OR729") ~ "Colorado_new",
    str_detect(sample_id, "MH384") ~ "Australia_D_mel",
    str_detect(sample_id, "MW97")  ~ "Australia_D_sim",
    str_detect(sample_id, "KX88")  ~ "China_Diptera",
    str_detect(sample_id, "OR820") ~ "OC",
    str_detect(sample_id, "EVE")   ~ "EVE",
    .default = "Others"
  ))

# how should tree tip labels be formatted
# colors here are RGB
sample_ids <- sample_ids %>% mutate(label_format = case_when(
  sample_type == "Maine"               ~ "label #CC6600 normal 1",
  sample_type == "Pennsylvania"        ~ "label #CC6600 normal 1",
  str_detect(sample_type, "Colorado")  ~ "label #0000FF normal 1",
  sample_type == "OC"                  ~ "label #000000 normal 1",
  sample_type == "EVE"                 ~ "label #006600 normal 1",
  str_detect(sample_type, "Australia") ~ "label #A0A0A0 normal 1",
  str_detect(sample_type, "China")     ~ "label #A0A0A0 normal 1",
  sample_type == "Others"              ~ "label #A0A0A0 normal 1",
  .default    =                        "label #FFFF33 normal 1" 
))

# how should branches be formatted?
sample_ids <- sample_ids %>% mutate(branch_format = case_when(
  sample_type == "Maine"               ~ "branch #CC6600 normal 1",
  sample_type == "Pennsylvania"        ~ "branch #CC6600 normal 1",
  str_detect(sample_type, "Colorado")  ~ "branch #0000FF normal 1",
  sample_type == "OC"                  ~ "branch #000000 normal 1",
  sample_type == "EVE"                 ~ "branch #006600 normal 1",
  str_detect(sample_type, "Australia") ~ "branch #A0A0A0 normal 1",
  str_detect(sample_type, "China")     ~ "branch #A0A0A0 normal 1",
  sample_type == "Others"              ~ "branch #A0A0A0 normal 1",
  .default    =                        "branch #FFFF33 normal 1" 
))

# how many of each type?
sample_ids %>% group_by(sample_type) %>% summarize(n=n())

# write out header lines to make this an iTOL annotation file
# in the style of colors_styles_template.txt
write_itol_colors_style_header <- function(filename){
  fileConn<-file(filename)
  writeLines(c("TREE_COLORS","SEPARATOR SPACE","DATA"), fileConn)
  close(fileConn)
}

# write out label format info to space-delimited file
itol_labels <- sample_ids %>% select(sequence_id, label_format)
filename <- "galbut_RNA1_label_format.txt"
write_itol_colors_style_header(filename)
write.table(itol_labels, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

# write out branch format info to space-delimited file
itol_branch <- sample_ids %>% select(sequence_id, branch_format)
filename <- "galbut_RNA1_branch_format.txt"
write_itol_colors_style_header(filename)
write.table(itol_branch, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)


