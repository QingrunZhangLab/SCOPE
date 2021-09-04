library(data.table)
library(stringr)

# Default values for parameters
DIFFCORR_PERCENTILE_THRESHOLD = 0.975

# This script, just like 6_pathway_overlaps.R requires all pathways to be
# in the "pathways/" directory.
all_pathways <- c()
for(curr_file in list.files("pathways/", full.names = TRUE)){
  curr_data <- fread(curr_file)
  all_pathways <- rbindlist(list(all_pathways, curr_data))
}

all_corrs <- c()
# It also requires the SecondaryGenes files to be in the "SecondaryGenes/"
# directory. Null distributions of differential co-expression should also 
# be in their own directory.
for(curr_file in list.files("SecondaryGenes/", full.names = TRUE)){
  curr_data <- fread(curr_file)
  
  # Relies on the cancer or disease to be the first part of the file name
  # For example: BRCA_SecondaryGenes.csv
  curr_cancer <- str_split(basename(curr_file), pattern = "_", simplify = TRUE)[1]
  curr_data$Cancer <- curr_cancer
  
  # Relies on the null distribution files also to use the same pattern of
  # file name. These can be easily modified as required.
  diffcorr_null_dist <- fread(paste0("NullDistributions/", curr_cancer, "_2222_100_1000_DiffCorrNull.csv"))
  
  diffcorr_cutoff <- c()
  
  print("Calculating median percentile from differential co-expression data...")
  for(i in 1:ncol(diffcorr_null_dist)){
    
    diffcorr_cutoff <- c(diffcorr_cutoff, quantile(diffcorr_null_dist[, ..i][[1]], DIFFCORR_PERCENTILE_THRESHOLD, na.rm = TRUE)[[1]])
    
  }
  
  diffcorr_cut <- median(diffcorr_cutoff)
  
  curr_data$`Differential Correlation Sig` <- diffcorr_cut
  curr_data$`Is Differential Correlation Sig` <- ifelse(curr_data$`Differential Correlation` > curr_data$`Differential Correlation Sig`, "Yes", "No")
  
  all_corrs <- rbindlist(list(all_corrs, curr_data))
}

final_data <- c()
for(i in 1:nrow(all_pathways)){
  genes <- unique(str_split(all_pathways$userId[i], ";")[[1]])
  curr_data <- data.table(Cancer = rep(all_pathways$Cancer[i], length(genes)),
                          `KEGG GeneSet` = rep(all_pathways$geneSet[i], length(genes)),
                          `Pathway Class` = rep(all_pathways$pathway_class[i], length(genes)),
                          `Pathway Name` = rep(all_pathways$description[i], length(genes)),
                          `Core Gene` = rep(all_pathways$`Core Gene`[i], length(genes)),
                          `Secondary Gene` = genes)
  
  final_data <- rbindlist(list(final_data, curr_data))
}

final_data$`Core Gene Name` <- all_corrs$`Core Gene Name`[match(final_data$`Core Gene`, all_corrs$`Core Gene`)]
final_data$`Secondary Gene Name` <- all_corrs$`Secondary Gene Name`[match(final_data$`Secondary Gene`, all_corrs$`Secondary Gene`)]

all_diffs <- fread("../AllDiffs.csv")
all_diffs$chromosome_name <- NULL

final_data <- merge(final_data, all_corrs, by = c("Cancer", "Core Gene", "Secondary Gene", "Core Gene Name", "Secondary Gene Name"), no.dups = TRUE)

final_data <- merge(final_data, all_diffs, by.x = c("Cancer", "Secondary Gene", "Secondary Gene Name"),
                    by.y = c("Cancer", "ensembl_gene_id", "external_gene_name"))

final_data <- merge(final_data, all_diffs, by.x = c("Cancer", "Core Gene", "Core Gene Name"),
                    by.y = c("Cancer", "ensembl_gene_id", "external_gene_name"), suffixes = c("", "_Core"))

fwrite(final_data, "SCOPE_Gene_Level.csv")
