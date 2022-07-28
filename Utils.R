results_to_datatable <- function(results, core.genes, core.pathways, core.pathway.classes, extra.causal.genes, dataset.info, extra.info){
  
  if(length(results$selected_genes)==1 & results$selected_genes[1]=="na"){
    selected_genes <- c()
  }else{
    selected_genes <- results$selected_genes
  }
  if((is.data.table(results$selected_pathways) && nrow(results$selected_pathways)==0) || length(results$selected_pathways)==1 && length(results$selected_pathways[1])==1 && results$selected_pathways[1]=="na"){
    selected_pathways <- c()
  }else{
    selected_pathways <- results$selected_pathways
  }
  
  return_table <- as.data.table(c(dataset.info, extra.info))
  
  core.genes <- str_split(core.genes, ";")[[1]]
  extra.causal.genes <- str_split(extra.causal.genes, ";")[[1]]
  core.pathways <- str_split(core.pathways, ";")[[1]]
  core.pathway.classes <- str_split(core.pathway.classes, "\\|")[[1]]
  
  #### GENES HIT
  cores_hit <- sum(selected_genes %in% core.genes)
  causal_hit <- sum(selected_genes %in% extra.causal.genes)
  
  return_table$cores_hit <- cores_hit
  return_table$causal_hit <- causal_hit
  
  #### CORE PATHWAYS
  pathways_hit <- unique(selected_pathways$geneSet)
  core_pathways_hit <- sum(pathways_hit %in% core.pathways)
  
  return_table$core_pathways_hit <- core_pathways_hit
  
  #### CORE PATHWAY CLASSES
  pathway_classes_hit <- unique(selected_pathways$pathway_class)
  core_pathway_classes <- unique(core.pathway.classes)
  core_pathway_classes_count <- length(core_pathway_classes)
  core_pathway_classes_hit <- sum(pathway_classes_hit %in% core_pathway_classes)
  
  return_table$core_pathway_classes_count <- core_pathway_classes_count
  return_table$core_pathway_classes_hit <- core_pathway_classes_hit
  
  #### EXTRAS
  extra_genes <- selected_genes[!(selected_genes %in% c(core.genes, extra.causal.genes))]
  extra_genes_count <- length(extra_genes)
  
  extra_pathways <- unique(selected_pathways$geneSet[!(selected_pathways$geneSet %in% core.pathways)])
  extra_pathway_count <- length(extra_pathways)
  
  extra_pathway_classes <- unique(selected_pathways$pathway_class[!(selected_pathways$pathway_class %in% core_pathway_classes)])
  extra_pathway_classes_count <- length(extra_pathway_classes)
  
  return_table$extra_genes_count <- extra_genes_count
  return_table$extra_pathway_count <- extra_pathway_count
  return_table$extra_pathway_classes_count <- extra_pathway_classes_count
  
  return(return_table)
}


check_files_stop <- function(file_names){
  files_present <- file.exists(file_names)
  if(any(!files_present)){
    stop(paste0("The following files do not exist:\n", paste0(file_names[!files_present], collapse = "\n"), "\nStopping run.\n"))
  } 
}


check_dirs_stop <- function(dir_names){
  dirs_present <- dir.exists(dir_names)
  if(any(!dirs_present)){
    stop(paste0("The following directories do not exist:\n", paste0(dir_names[!dirs_present], collapse = "\n"), "\nStopping run.\n"))
  } 
}


check_dirs_create <- function(dir_names){
  dirs_present <- dir.exists(dir_names)
  if(any(!dirs_present)){
    sapply(dir_names[!dirs_present], dir.create, recursive = TRUE)
    cat(paste0("The following directories did not exist and were created:\n", 
               paste0(dir_names[!dirs_present], collapse = "\n"), "\n"))
  } 
}


remove_zero_variance <- function(expr_data){
  vars <- sapply(expr_data, var)
  non_constants <- names(vars)[vars!=0 | is.na(vars)]
  return(expr_data[, ..non_constants])
}


remove_low_variance <- function(expr_data, quant){
  vars <- sapply(expr_data, var)
  high_vars <- names(vars)[is.na(vars) | vars >= quantile(vars, quant, na.rm = TRUE)]
  return(expr_data[, ..high_vars])
}


flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}