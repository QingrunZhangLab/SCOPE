library(data.table)
library(biomaRt)

# This script will properly format and sort the outputs from the original lasso runs.
# Requires specifying the number of iterations used in the original runs.
ITERS = 1000

# Change biomart Mart according to requirement
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Assumes all '_Summary.csv' files are in the directory below to be formatted
for(x in list.files("LassoResultsSummaries/", full.names = TRUE)){
  curr_data <- fread(x)
  gene_coords <- as.data.table(getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                                     filters="ensembl_gene_id", values=curr_data$probe, mart=human))
  curr_data$gene_name <- gene_coords$external_gene_name[match(curr_data$probe, gene_coords$ensembl_gene_id)]
  curr_data$prop <- curr_data$model.count/ITERS
  curr_data$model.count <- NULL
  
  curr_data <- setcolorder(curr_data, c("probe", "gene_name", "prop", "min.val", "max.val"))
  colnames(curr_data) <- c("Probe", "Gene Name", "Proportion", "Minimum Coefficient", "Maximum Coefficient")
  curr_data <- curr_data[order(Proportion, decreasing = TRUE)]
  
  fwrite(curr_data, paste0("Formatted_", basename(x)))
}