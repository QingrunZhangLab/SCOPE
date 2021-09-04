library(data.table)
library(WebGestaltR)
library(KEGGREST)

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
  stop("Secondary gene file must be specified!")
}else{
  SECONDARY_GENE_FILE = args[1]
  CANCER = args[2]
}

full_data <- fread(SECONDARY_GENE_FILE)

# For each core gene, the CGN (core gene + secondary genes) are
# submitted for pathway enrichment using ORA.
for(curr_core_gene in unique(full_data$`Core Gene`)){
  
  curr_gene_list <- c(curr_core_gene, full_data$`Secondary Gene`[full_data$`Core Gene` == curr_core_gene])
  
  enrichment_results <- NULL
  
  tryCatch(enrichment_results <- as.data.table(WebGestaltR(interestGene = curr_gene_list,
                                                           organism = "hsapiens",
                                                           enrichDatabase="pathway_KEGG",
                                                           interestGeneType="ensembl_gene_id", referenceSet = "genome",
                                                           referenceGeneType = "ensembl_gene_id", isOutput = FALSE)),
           error = function(c){
             print("Error in enrichment. Skipping...")
             enrichment_results <<- NULL
           })
  
  if(is.null(enrichment_results) | nrow(enrichment_results) == 0){
    print("No pathways enriched. Skipping...")
    next
  } 
  
  enrichment_results$`Core Gene` <- curr_core_gene
  enrichment_results$`Core Gene Name` <- full_data$`Core Gene Name`[match(curr_core_gene, full_data$`Core Gene`)]
  enrichment_results$Cancer <- CANCER
  
  if(exists("final_results")){
    final_results <- rbindlist(list(final_results, enrichment_results))
  }else{
    final_results <- enrichment_results
  }
}

final_results$pathway_class <- NA

pathways <- keggLink("pathway", "hsa")

# Pathway classes of each pathway discovered using enrichment are obtained
for(i in 1:nrow(final_results)){
  query <- keggGet(paste0("path:", final_results$geneSet[i]))
  final_results$pathway_class[i] <- ifelse(exists("CLASS", query[[1]]), query[[1]]$CLASS, "N/A")
}

fwrite(final_results, paste0(CANCER, "_Pathways.csv"))
