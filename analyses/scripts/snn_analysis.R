library(tidyverse)
library(ape)
library(patchwork)

# this script reads in sets of galbut virus sequences and 
# calculates Snn, the shared location nearest neighbor statistics 
# according to Hudson (2011) 
#
# Hudson RR. A new statistic for detecting genetic differentiation. 
# Genetics. 2000 Aug;155(4):2011-4. 
# doi: 10.1093/genetics/155.4.2011. PMID: 10924493; PMCID: PMC1461195.
#
# Mark Stenglein 8/28/2025

# -------------------------------
# Metadata (location of samples)
# -------------------------------

# first, read in map of accession --> sample ID
id_acc_map <- read.delim("analyses/data/sample_accession_map.txt", header=F)
colnames(id_acc_map) <- c("sample_id", "segment", "accession", "full_name")

# second, read in sample metadata including location
metadata <- read_excel("analyses/data/Metadata_Table.xlsx")

# third, location metadata for genbank sequences
genbank_metadata <- read_excel("analyses/data/genbank_metadata.xlsx")

# pull out just the columns we need
id_acc_map <- id_acc_map %>% select(sample_id, accession)
metadata <- metadata %>% select(sample_id = sample_name, location)
genbank_metadata <- genbank_metadata %>% select(accession, location)

# make consistent names w/r/t _ and -
metadata$sample_id <- str_replace_all(metadata$sample_id, "_", "-")

# join (merge) metadata
location_metadata <- left_join(id_acc_map, metadata) %>% select(accession, location)

# add in genbank metadata (rbind)
# to get a table of acccession -> location
location_metadata <- rbind(location_metadata, genbank_metadata)


num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "analyses/trees/alignments/RNA1_nucleotide_alignment.fasta", 
    prefix="RNA1") {
    
    # DEBUG
    fasta_msa = "analyses/trees/alignments/Chaq_nucleotide_alignment.fasta"
    prefix="Chaq"
      
    # read in sequences
    msa <- read.dna(fasta_msa, format="fasta")
    
    # relabel sequences with just the first part up to whitespace, i.e. the accession
    accessions <- str_match(labels(msa), "(\\S+) ")[,2]
    
    # update labels: shorten to just accessions
    rownames(msa) <- accessions
    
    # drop sequences for which we could not identify an accession
    msa <- msa[!(is.na(accessions)), ]
    accessions <- accessions[!(is.na(accessions))]
    
    # calculate distances
    # TN93 distance model
    dist_model = "TN93" 
    msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
    
    acc_loc <- tibble(accession = rownames(msa_dist))
    # this will make locations in order of matrix accession
    acc_loc <- left_join(acc_loc, location_metadata)  
    # create a vector of locations matching the order of the distance matrix
    locations <- acc_loc %>% pull(location)
    
    # calculate Snn from distance matrix and location metadata
    this_snn <- calculate_snn(msa_dist, locations)
    
    perm_snns <- c()
    
    # do permutation testing
    for (p in 1:num_permutations) {
      # first, scramble locations using sampling without replacement
      # using sampling without replacement will keep the # of samples from each
      # location the same between permutations, as specified in Hudson 2011
      scrambled_locations  <- sample(locations, replace=F)
      perm_snn             <- calculate_snn(msa_dist, scrambled_locations)
      perm_snns            <- c(perm_snns, perm_snn)
    }
    
    # calculate p-value from Snns from permutation tests
    this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
    
    if (this_snn_p_val == 0) {
      this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
      this_snn_p_val <- paste0("<", this_snn_p_val_sci)
    }
    
    snn_text <- paste0("p = ", 
                       this_snn_p_val)
    print(snn_text)
    
    return(list(snn           = this_snn,
                permuted_snns = perm_snns,
                pvalue        = this_snn_p_val,
                text          = snn_text))
}

plot_snns <- function( snns_df ) {
    
  snns_df$segment  <- fct_relevel(snns_df$segment, "Chaq", after = 3)
  ggplot(snns_df) +
    geom_point  (data = filter(snns_df, snn_type == "calculated"),
                 aes(x=segment, y = snn) , 
                 shape=21, size=2, 
                 fill="slateblue", color="black", stroke=0.25) + 
    geom_jitter(data = filter(snns_df, snn_type == "permuted"),
                aes(x=segment, y = snn),
                shape=21, size=2, 
                fill="grey80", color="black", stroke=0.25,
                alpha=0.25,
                width=0.25, height=0) +
    # geom_violin(data = filter(snns_df, snn_type == "permuted"),
                # aes(x=segment, y = snn),
                # fill=NA, color="black", linewidth=.25,
                # width=0.75) +
    xlab("") + 
    ylab("Nearest neighbor statistic, Snn") +
    ylim(c(0,1)) +
    theme_bw(base_size = 14)  +
    theme(
      plot.title = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none")
}

calculate_snn <- function(distance_matrix, locations_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = chaq_msa_dist
  # locations_vector = chaq_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same location as sample i
    num_shared_loc <- sum(locations_vector[nn_i] == locations_vector[i])
    
    # this is the number of nearest neighbors (regardless of location)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_loc <- (num_shared_loc / num_NN) / num_samples
    
    if (is.na(fraction_shared_loc)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_loc
  }
  
  # return Snn
  Snn 
}

# caclulate SNNs and do permutation testing for all segments
rna1_snn <- process_alignment("analyses/trees/alignments/RNA1_nucleotide_alignment.fasta", "RNA1")
rna2_snn <- process_alignment("analyses/trees/alignments/RNA2_nucleotide_alignment.fasta", "RNA2")
rna3_snn <- process_alignment("analyses/trees/alignments/RNA3_nucleotide_alignment.fasta", "RNA3")
chaq_snn <- process_alignment("analyses/trees/alignments/Chaq_nucleotide_alignment.fasta", "Chaq")

# make a combined tibble of SNN values (calculated and from permuted labels)
snns_df_calc_snn <- tibble(segment = c("RNA1", "RNA2", "RNA3", "Chaq"),
                           snn_type = "calculated",
                           snn = c(rna1_snn$snn, rna2_snn$snn, rna3_snn$snn, chaq_snn$snn))
snns_df_rna1_perm <- tibble(segment ="RNA1",
                            snn_type = "permuted",
                            snn = rna1_snn$permuted_snns)
snns_df_rna2_perm <- tibble(segment ="RNA2",
                            snn_type = "permuted",
                            snn = rna2_snn$permuted_snns)
snns_df_rna3_perm <- tibble(segment ="RNA3",
                            snn_type = "permuted",
                            snn = rna3_snn$permuted_snns)
snns_df_chaq_perm <- tibble(segment ="Chaq",
                            snn_type = "permuted",
                            snn = chaq_snn$permuted_snns)

snns_df <- rbind(snns_df_calc_snn, snns_df_rna1_perm, snns_df_rna2_perm, snns_df_rna3_perm, snns_df_chaq_perm)

plot_snns(snns_df)
ggsave("./analyses/plots/supplemental_figure_X_snn.pdf", width=7.5, height=5, units="in")

snn_supp_table <- 
  tibble( segment = c("galbut virus RNA1", "galbut virus RNA 2", "galbut virus RNA 3", "chaq virus"),
          snn  = c(rna1_snn$snn, rna2_snn$snn, rna3_snn$snn, chaq_snn$snn),
          pval = c(rna1_snn$pval, rna2_snn$pval, rna3_snn$pval, chaq_snn$pval))

write.table(snn_supp_table, 
            file="analyses/plots/supplemental_table_X_snn.txt",
            quote=F, sep="\t", row.names=F)

# text for paper
output_text <- 
  paste0(
    "The Snn for galbut virus RNA 1, RNA 2, RNA 3, and chaq virus were ",
    sprintf("%0.2f, ", rna1_snn$snn),
    sprintf("%0.2f, ", rna2_snn$snn),
    sprintf("%0.2f, ", rna3_snn$snn),
    sprintf("and %0.2f, respectively. ", chaq_snn$snn),
    "These high Snn values are consistent with a high degree of geographical population structure ",
    "and were all statistically signifcant with p-values ",
    rna1_snn$pvalue # one p-value OK here since all the same
  )

print(output_text)

# make a combined plot
# combined_snn_plot <-  rna1_snn$plot + rna2_snn$plot +  rna3_snn$plot +  chaq_snn$plot +
  # plot_layout(ncol = 2) +
  # plot_annotation(tag_levels = 'A')

# combined_snn_plot
# ggsave("./analyses/plots/supplemental_figure_X_snn.pdf", width=7.5, height=6, units="in")
