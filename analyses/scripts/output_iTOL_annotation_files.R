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

RNA1 <- read.csv("analyses/trees/RNA1/RNA1_Seq_IDs.csv")
RNA2 <- read.csv("analyses/trees/RNA2/RNA2_Seq_IDs.csv")
RNA3 <- read.csv("analyses/trees/RNA3/RNA3_SeqIDs.csv")
Chaq <- read.csv("analyses/trees/Chaq/Chaq_Seq_IDs.csv")

RNA1$sample_id <- paste(RNA1$location, RNA1$date, RNA1$accession, sep = "_")
RNA1_ids <- RNA1 %>% 
  select(accession, sample_id)

RNA2$sample_id <- paste(RNA2$location, RNA2$date, RNA2$accession, sep = "_")
RNA2_ids <- RNA2 %>% 
  select(accession, sample_id)

RNA3$sample_id <- paste(RNA3$location, RNA3$date, RNA3$accession, sep = "_")
RNA3_ids <- RNA3 %>% 
  select(accession, sample_id)

Chaq$sample_id <- paste(Chaq$location, Chaq$date, Chaq$accession, sep = "_")
Chaq_ids <- Chaq %>% 
  select(accession, sample_id)

# write out header lines to make this an iTOL annotation file
# in the style of labels_template.txt
write_itol_labels_header <- function(filename){
  fileConn<-file(filename)
  writeLines(c("LABELS","SEPARATOR SPACE","DATA"), fileConn)
  close(fileConn)
}
filename <- "galbut_RNA1_label_format.txt"
write_itol_labels_header(filename)
write.table(RNA1_ids, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

filename <- "galbut_RNA2_label_format.txt"
write_itol_labels_header(filename)
write.table(RNA2_ids, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

filename <- "galbut_RNA3_label_format.txt"
write_itol_labels_header(filename)
write.table(RNA3_ids, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

filename <- "chaq_label_format.txt"
write_itol_labels_header(filename)
write.table(Chaq_ids, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

# decide on label and branch colors
# based on sequence ID

# how should tree tip labels be formatted
# colors here are RGB
RNA1 <- RNA1 %>% mutate(label_format = case_when(
  location == "Maine"               ~ "label #f44336 normal 1",
  location == "Pennsylvania"        ~ "label #5b5b5b normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #351c75 normal 1",
  location == "Minnesota"           ~ "label #134f5c normal 1",
  location == "California"          ~ "label #0b5394 normal 1",
  accession == "EVE"                ~ "label #006600 normal 1",
  location == "Australia"           ~ "label #b45f06 normal 1",
  location == "China"               ~ "label #bf9000 normal 1",
))

RNA2 <- RNA2 %>% mutate(label_format = case_when(
  location == "Maine"               ~ "label #f44336 normal 1",
  location == "Pennsylvania"        ~ "label #5b5b5b normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #351c75 normal 1",
  location == "Minnesota"           ~ "label #134f5c normal 1",
  location == "California"          ~ "label #0b5394 normal 1",
  location == "Australia"           ~ "label #b45f06 normal 1",
))

RNA3 <- RNA3 %>% mutate(label_format = case_when(
  location == "Maine"               ~ "label #f44336 normal 1",
  location == "Pennsylvania"        ~ "label #5b5b5b normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #351c75 normal 1",
  location == "Minnesota"           ~ "label #134f5c normal 1",
  location == "California"          ~ "label #0b5394 normal 1",
  location == "Australia"           ~ "label #b45f06 normal 1",
))

Chaq <- Chaq %>% mutate(label_format = case_when(
  location == "Maine"               ~ "label #f44336 normal 1",
  location == "Pennsylvania"        ~ "label #5b5b5b normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #351c75 normal 1",
  location == "Minnesota"           ~ "label #134f5c normal 1",
  location == "California"          ~ "label #0b5394 normal 1",
  location == "Australia"           ~ "label #b45f06 normal 1",
))

# how should branches be formatted?
RNA1 <- RNA1 %>% mutate(branch_format = case_when(
  date == "1963"                    ~ "label #000000 normal 1",
  location == "Maine"               ~ "label #CC6600 normal 1",
  location == "Pennsylvania"        ~ "label #CC6600 normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #000000 normal 1",
  location == "Minnesota"           ~ "label #000000 normal 1",
  location == "California"          ~ "label #000000 normal 1",
  accession == "EVE"                ~ "label #006600 normal 1",
  location == "Australia"           ~ "label #A0A0A0 normal 1",
  location == "China"               ~ "label #A0A0A0 normal 1",
))

RNA2 <- RNA2 %>% mutate(branch_format = case_when(
  date == "1963"                    ~ "label #000000 normal 1",
  location == "Maine"               ~ "label #CC6600 normal 1",
  location == "Pennsylvania"        ~ "label #CC6600 normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #000000 normal 1",
  location == "Minnesota"           ~ "label #000000 normal 1",
  location == "California"          ~ "label #000000 normal 1",
  location == "Australia"           ~ "label #A0A0A0 normal 1",
))

RNA3 <- RNA3 %>% mutate(branch_format = case_when(
  date == "1963"                    ~ "label #000000 normal 1",
  location == "Maine"               ~ "label #CC6600 normal 1",
  location == "Pennsylvania"        ~ "label #CC6600 normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #000000 normal 1",
  location == "Minnesota"           ~ "label #000000 normal 1",
  location == "California"          ~ "label #000000 normal 1",
  location == "Australia"           ~ "label #A0A0A0 normal 1",
))

Chaq <- Chaq %>% mutate(branch_format = case_when(
  location == "Maine"               ~ "label #CC6600 normal 1",
  location == "Pennsylvania"        ~ "label #CC6600 normal 1",
  location == "Colorado"            ~ "label #0000FF normal 1",
  location == "Ohio"                ~ "label #CC6600 normal 1",
  location == "Illinois"            ~ "label #000000 normal 1",
  location == "Minnesota"           ~ "label #000000 normal 1",
  location == "California"          ~ "label #000000 normal 1",
  location == "Australia"           ~ "label #A0A0A0 normal 1",
))

# write out header lines to make this an iTOL annotation file
# in the style of colors_styles_template.txt
write_itol_colors_style_header <- function(filename){
  fileConn<-file(filename)
  writeLines(c("TREE_COLORS","SEPARATOR SPACE","DATA"), fileConn)
  close(fileConn)
}

# write out label format info to space-delimited file
itol_colors_RNA1 <- RNA1 %>% select(accession, label_format)
filename <- "galbut_RNA1_color_format.txt"
write_itol_colors_style_header(filename)
write.table(itol_colors_RNA1, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

itol_colors_RNA2 <- RNA2 %>% select(accession, label_format)
filename <- "galbut_RNA2_color_format.txt"
write_itol_colors_style_header(filename)
write.table(itol_colors_RNA2, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

itol_colors_RNA3 <- RNA3 %>% select(accession, label_format)
filename <- "galbut_RNA3_color_format.txt"
write_itol_colors_style_header(filename)
write.table(itol_colors_RNA3, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

itol_colors_Chaq <- Chaq %>% select(accession, label_format)
filename <- "Chaq_color_format.txt"
write_itol_colors_style_header(filename)
write.table(itol_colors_Chaq, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)


# write out branch format info to space-delimited file
itol_branch_RNA1 <-RNA1 %>% select(accession, branch_format)
filename <- "galbut_RNA1_branch_format.txt"
write_itol_colors_style_header(filename)
write.table(itol_branch_RNA1, file=filename, append=T, sep=" ", col.names=F, row.names=F, quote=F)

write.csv(RNA1, "analyses/trees/RNA1/RNA1_annotation_metadata.csv", row.names = FALSE, quote = FALSE)
write.csv(RNA2, "analyses/trees/RNA2/RNA2_annotation_metadata.csv", row.names = FALSE, quote = FALSE)
write.csv(RNA3, "analyses/trees/RNA3/RNA3_annotation_metadata.csv", row.names = FALSE, quote = FALSE)
write.csv(Chaq, "analyses/trees/Chaq/Chaq_annotation_metadata.csv", row.names = FALSE, quote = FALSE)
