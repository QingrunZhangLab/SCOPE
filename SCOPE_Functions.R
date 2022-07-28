null_corr <- function(expr.data, iters = 100, probes.per.iter = 1000, quantile.cutoff = 0.975){
  positive_cutoffs <- c()
  negative_cutoffs <- c()
  
  for(i in 1:iters){
    rnd_smp <- sample(colnames(expr.data), probes.per.iter)
    
    cor_mat <- cor(expr.data[, ..rnd_smp])
    corrs <- flattenCorrMatrix(cor_mat)
    
    positive_cutoffs <- c(positive_cutoffs, quantile(corrs$cor[corrs$cor > 0], quantile.cutoff))
    negative_cutoffs <- c(negative_cutoffs, -1*quantile(-1*corrs$cor[corrs$cor < 0], quantile.cutoff))
  }
  
  return(list(positive_cut = median(positive_cutoffs),
              negative_cut = median(negative_cutoffs)))
}

null_diffcorr <- function(expr.data, phen, iters = 100, probes.per.iter = 1000, quantile.cutoff = 0.975){
  expr.data$phen <- phen
  
  diffcorr_cutoffs <- c()
  
  case <- expr.data[phen == 1,]
  control <- expr.data[phen == 0,]
  
  for(diffcorr_iter in 1:iters){
    
    rnd_smp <- sample(colnames(expr.data)[1:(ncol(expr.data)-1)], probes.per.iter)
    
    case_sel <- case[, ..rnd_smp]
    control_sel <- control[, ..rnd_smp]
    
    case_cor <- cor(case_sel)
    control_cor <- cor(control_sel)
    
    case_corflat <- flattenCorrMatrix(case_cor)
    control_corflat <- flattenCorrMatrix(control_cor)
    
    merged <- merge(case_corflat, control_corflat, by = c("row", "column"))
    merged$diffcor <- abs(merged$cor.x - merged$cor.y)
    
    diffcorr_cutoffs <- c(diffcorr_cutoffs, quantile(merged$diffcor, quantile.cutoff, na.rm = TRUE))
  }
  
  return(median(diffcorr_cutoffs))
}

slasso_run <- function(expr.data, phen, formula = NULL, 
                       seed_start = 2222,
                       iters = 1000, split_ratio = 0.7, 
                       cv.folds = 10, parallel = FALSE,
                       prop_cutoff = 0.8, 
                       c_pos_cutoff, c_neg_cutoff, 
                       d_cutoff, 
                       pathway.data,
                       log.file, 
                       info){
  
  expr.data$phen <- phen
  
  results <- c()
  # pred_metrics <- c()
  
  if(is.null(formula)){
    formula <- "phen ~ ."
  }
  
  set.seed(seed_start)
  for(k in 1:iters){
    
    # Splitting data proportional to phenotype composition
    inds  <- sample.split(expr.data$phen, SplitRatio = split_ratio)
    
    # print(table(expr.data$phen[inds]))
    
    cat(paste0("Currently fitting iteration: ", k, " (", info , ") \n"), file = log.file, append = TRUE)
    
    # tic("Iteration time")
    # Fitting a glmnet LASSO model for the training data
    cvfit <- glmnetUtils::cv.glmnet(formula = as.formula(formula), data = expr.data[inds,], alpha = 1,
                                    family = "binomial",  
                                    nfolds = cv.folds, intercept = FALSE,
                                    parallel = parallel)
    
    # Obtaining fitted coefficients
    tmp_coeffs <- coef(cvfit)
    # print("Fitted...")
    
    # Keeping track of results from each of the LASSO models fit in each iteration
    if(is.null(results)){
      results <- data.table(probe = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                            coefficient_1 = tmp_coeffs@x)
    }else{
      results <- merge(results, data.table(probe = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                                           coefficient = tmp_coeffs@x), 
                       all.x=TRUE, all.y=TRUE, 
                       by="probe", suffixes = c("",paste0("_", k)))
    }
    
  }
  
  # print("Combining LASSO results...")
  cof_sum <- data.table(probe = results$probe,
                        model.count = iters - rowSums(is.na(results[,2:ncol(results)])),
                        min.val = apply(results[,2:ncol(results)], 1, min, na.rm = TRUE),
                        max.val = apply(results[,2:ncol(results)], 1, max, na.rm = TRUE))
  
  identified_genes <- cof_sum$probe[(cof_sum$model.count/iters) > prop_cutoff]
  identified_genes <- unique(unlist(str_split(identified_genes, ":")))
  
  corr_results <- c()
  final_results <- c()
  
  case.data <- expr.data[phen == 1,]
  control.data <- expr.data[phen == 0,]
  
  cat(paste0("Currently calculating CGNs (", info , ") \n"), file = log.file, append = TRUE)
  if(length(identified_genes)>0){
    
    # For each core gene, co-expression values in the complete dataset,
    # and each of the phenotypes are calculated.
    for(i in 1:length(identified_genes)){
      
      core_gene <- identified_genes[i]
      g1 <- expr.data[[core_gene]]
      g1_case <- case.data[[core_gene]]
      g1_control <- control.data[[core_gene]]
      
      for(j in 1:(ncol(expr.data)-1)){
        
        curr_gene <- colnames(expr.data)[j]
        g2 <- expr.data[[curr_gene]]
        g2_case <- case.data[[curr_gene]]
        g2_control <- control.data[[curr_gene]]
        
        cor_test <- cor.test(g1, g2, method = "pearson")
        case_test <- cor.test(g1_case, g2_case, method = "pearson")
        control_test <- cor.test(g1_control, g2_control, method = "pearson")
        
        curr_result <- data.table(core_gene = core_gene,
                                  secondary_gene = curr_gene,
                                  correlation = cor_test$estimate,
                                  pval = cor_test$p.value,
                                  case.corr = case_test$estimate,
                                  case.pval = case_test$p.value,
                                  control.corr = control_test$estimate,
                                  control.pval = control_test$p.value)
        
        if(is.null(corr_results)){
          corr_results <- curr_result
        }else{
          corr_results <- rbindlist(list(corr_results, curr_result))
        }
      }
    }
    
    corr_results$corrdiff <- abs(corr_results$case.corr - corr_results$control.corr)
    corr_results <- corr_results[corr_results$core_gene != corr_results$secondary_gene, ]
    
    colnames(corr_results) <- c("Core Gene", "Secondary Gene", 
                                "Correlation", "Correlation p-value", "Case Correlation", "Case Correlation p-value",
                                "Control Correlation", "Control Correlation p-value", "Differential Correlation")
    
    # Identifying CGNs using secondary genes that are co-expressed or 
    # differentially co-expressed.
    filtered_results <- corr_results[corr_results$Correlation > c_pos_cutoff | 
                                       corr_results$Correlation < c_neg_cutoff | 
                                       corr_results$`Differential Correlation` > d_cutoff,]
    
    cat(paste0("Currently conducting SLASSO pathway enrichment (", info , ") \n"), file = log.file, append = TRUE)
    for(curr_core_gene in unique(filtered_results$`Core Gene`)){
      
      curr_gene_list <- c(curr_core_gene, filtered_results$`Secondary Gene`[filtered_results$`Core Gene` == curr_core_gene])
      curr_gene_list <- gsub("\\..*","", curr_gene_list)
      
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
      
      if(is.null(enrichment_results) || nrow(enrichment_results) == 0){
        print("No pathways enriched. Skipping...")
        next
      } 
      
      enrichment_results$`Core Gene` <- curr_core_gene
      
      final_results <- rbindlist(list(final_results, enrichment_results))
    }
    
    if(!is.null(final_results)){
      final_results$pathway_class <- pathway.data$PathwayClass[match(final_results$geneSet, pathway.data$PathwayID)]
    }
  }
  
  
  return_results <- list()
  if(length(identified_genes)==0){
    return_results[["selected_genes"]] <- c("na")
  }else{
    return_results[["selected_genes"]] <- identified_genes
  }
  if(!is.null(final_results)){
    return_results[["selected_pathways"]] <- final_results
  }else{
    return_results[["selected_pathways"]] <- c("na")
  }
  
  return(return_results)
}


