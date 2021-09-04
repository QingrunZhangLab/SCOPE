library(data.table)
library(biomaRt)

args = commandArgs(trailingOnly=TRUE)

# Default values for parameters
CORE_CUTOFF = 950
CORR_PERCENTILE_THRESHOLD = 0.975
DIFFCORR_PERCENTILE_THRESHOLD = 0.975
PREFIX = "stabilized_lasso"

if(length(args)<4){
  stop("Expression file, summary file and null distribution files MUST be specified.")
}else if(length(args)==4){
  EXPR_FILE = args[1]
  SUMMARY_FILE = args[2]
  CORR_NULL_FILE = args[3]
  DIFFCORR_NULL_FILE = args[4]
}else{
  EXPR_FILE = args[1]
  SUMMARY_FILE = args[2]
  CORE_CUTOFF = as.numeric(args[3])
  CORR_NULL_FILE = args[4]
  CORR_PERCENTILE_THRESHOLD = as.numeric(args[5])
  DIFFCORR_NULL_FILE = args[6]
  DIFFCORR_PERCENTILE_THRESHOLD = as.numeric(args[7])
  PREFIX = args[8]
}

print(paste0("Loading gene summary file: ", SUMMARY_FILE))
gene_summary <- fread(SUMMARY_FILE)

print("Filtering for core genes...")
core_genes <- gene_summary$probe[gene_summary$model.count > CORE_CUTOFF & gene_summary$probe != "(Intercept)"]

print("Setting up biomaRt to get external gene names...")
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords <- as.data.table(getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                                   filters="ensembl_gene_id", values=core_genes, mart=human))

print("Reading null distribution of co-expressions...")
corr_null_dist <- fread(CORR_NULL_FILE)

positive_cutoff <- c()
negative_cutoff <- c()

print("Calculating median percentile from co-expression data...")
for(i in 1:ncol(corr_null_dist)){
  
  positive_cutoff <- c(positive_cutoff, quantile(corr_null_dist[as.logical(corr_null_dist[, ..i] > 0), ..i][[1]], CORR_PERCENTILE_THRESHOLD)[[1]])
  negative_cutoff <- c(negative_cutoff, -quantile(-corr_null_dist[as.logical(corr_null_dist[, ..i] < 0), ..i][[1]], CORR_PERCENTILE_THRESHOLD)[[1]])

}

pos_cut <- median(positive_cutoff)
neg_cut <- median(negative_cutoff)

print("Reading null distribution of differential co-expressions...")
diffcorr_null_dist <- fread(DIFFCORR_NULL_FILE)

diffcorr_cutoff <- c()

print("Calculating median percentile from differential co-expression data...")
for(i in 1:ncol(diffcorr_null_dist)){
  
  diffcorr_cutoff <- c(diffcorr_cutoff, quantile(diffcorr_null_dist[, ..i][[1]], DIFFCORR_PERCENTILE_THRESHOLD, na.rm = TRUE)[[1]])
  
}

diffcorr_cut <- median(diffcorr_cutoff)

print("Loading expression data to identify secondary genes...")
expr.data <- fread(EXPR_FILE)
print(paste0("Data loaded. Dimensions: ", paste0(dim(expr.data), collapse = ",")))

# Filtering low variance and constant expression values
print(paste0("Calculating variances..."))
vars <- sapply(expr.data, var)
constants <- names(vars)[vars!=0 | is.na(vars)]
print(paste0("Identified ", sum(vars==0, na.rm=TRUE), " transcripts with 0 variance. Removing..."))
expr.data <- expr.data[, ..constants]
print(paste0("Remaining data dimensions: ", paste0(dim(expr.data), collapse = ",")))

print(paste0("Calculating variances..."))
vars <- sapply(expr.data, var)
low_vars <- names(vars)[is.na(vars) | vars >= quantile(vars, 0.25, na.rm = TRUE)]
print(paste0("Identified ", (ncol(expr.data) - length(low_vars) - 2), " transcripts with lowest 25% variance. Removing..."))
expr.data <- expr.data[, ..low_vars]
print(paste0("Remaining data dimensions: ", paste0(dim(expr.data), collapse = ",")))

# Taking a look at available phenotypes
print("Current phenotypes: ")
print(table(expr.data$phen))

# Selecting only Primary Tumor and Solid Tissue Normal samples
print("Filtering to include only 'Primary Tumor' and 'Solid Tissue Normal' samples.")
expr.data <- expr.data[expr.data$phen %in% c("Primary Tumor", "Solid Tissue Normal")]
expr.data <- expr.data[, -c(1)]
tumor <- expr.data[expr.data$phen == "Primary Tumor",]
normal <- expr.data[expr.data$phen == "Solid Tissue Normal",]
expr.data <- expr.data[, !"phen"]

results <- c()

# For each core gene, co-expression values in the complete dataset,
# and each of the phenotypes are calculated.
for(i in 1:length(gene_coords$ensembl_gene_id)){
  
  core_gene <- gene_coords$ensembl_gene_id[i]
  g1 <- expr.data[[core_gene]]
  g1_tum <- tumor[[core_gene]]
  g1_norm <- normal[[core_gene]]
  
  for(j in 1:ncol(expr.data)){
    
    if(j %% 1000 == 0){
      print(paste0("Currently calculating for ", j, " of ", ncol(expr.data), " (Core Gene ", i, " of ", length(gene_coords$ensembl_gene_id), ")..."))  
    }
    
    curr_gene <- colnames(expr.data)[j]
    g2 <- expr.data[[curr_gene]]
    g2_tum <- tumor[[curr_gene]]
    g2_norm <- normal[[curr_gene]]
    
    cor_test <- cor.test(g1, g2, method = "pearson")
    tum_test <- cor.test(g1_tum, g2_tum, method = "pearson")
    norm_test <- cor.test(g1_norm, g2_norm, method = "pearson")
    
    curr_result <- data.table(core_gene = core_gene,
                              secondary_gene = curr_gene,
                              correlation = cor_test$estimate,
                              pval = cor_test$p.value,
                              tum.corr = tum_test$estimate,
                              tum.pval = tum_test$p.value,
                              norm.corr = norm_test$estimate,
                              norm.pval = norm_test$p.value)
    
    if(is.null(results)){
      results <- curr_result
    }else{
      results <- rbindlist(list(results, curr_result))
    }

  }
  
}

print("Calculating differential co-expression...")
results$corrdiff <- abs(results$tum.corr - results$norm.corr)
results <- results[results$core_gene != results$secondary_gene, ]

print("Setting up biomaRt to retrieve gene names...")
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords <- as.data.table(getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                                   filters="ensembl_gene_id", values=unique(c(results$core_gene, results$secondary_gene)), mart=human))
results$core_gene_name <- gene_coords$external_gene_name[match(results$core_gene, gene_coords$ensembl_gene_id)]
results$secondary_gene_name <- gene_coords$external_gene_name[match(results$secondary_gene, gene_coords$ensembl_gene_id)]
print(paste0("Unable to retrieve gene names for ",
             length(gene_coords$external_gene_name[gene_coords$external_gene_name==""]),
             " probes."))

results <- setcolorder(results, c("core_gene", "secondary_gene", "core_gene_name", "secondary_gene_name",
                                  "correlation", "pval", "tum.corr", "tum.pval", "norm.corr", "norm.pval",
                                  "corrdiff"))
colnames(results) <- c("Core Gene", "Secondary Gene", "Core Gene Name", "Secondary Gene Name",
                       "Correlation", "Correlation p-value", "Tumor Correlation", "Tumor Correlation p-value",
                       "Normal Correlation", "Normal Correlation p-value", "Differential Correlation")

# Identifying CGNs using secondary genes that are co-expressed or 
# differentially co-expressed.
filtered_results <- results[results$Correlation > pos_cut | results$Correlation < neg_cut | results$`Differential Correlation` > diffcorr_cut,]

print("Writing all results.")

fwrite(results, paste0(PREFIX, "_AllCoreCorrelations.csv"))
fwrite(filtered_results, paste0(PREFIX, "_SecondaryGenes.csv"))