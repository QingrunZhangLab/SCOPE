if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, 
               stringr,
               glmnetUtils,
               WebGestaltR,
               caTools)

flattenCorrMatrix <- function(corMat) {
  ut <- upper.tri(corMat)
  data.table(
    row = rownames(corMat)[row(corMat)[ut]],
    column = rownames(corMat)[col(corMat)[ut]],
    cor  = (corMat)[ut]
  )
}

removeZeroVariance <- function(exprData){
  vars <- sapply(exprData, var)
  non_constants <- names(vars)[vars!=0 | is.na(vars)]
  return(exprData[, ..non_constants])
}

removeLowVariance <- function(exprData, quantileToRemove = 0.25){
  vars <- sapply(exprData, var)
  high_vars <- names(vars)[is.na(vars) | vars >= quantile(vars, quantileToRemove, na.rm = TRUE)]
  return(exprData[, ..high_vars])
}

nullCorr <- function(exprData, iterations = 100, probesPerIter = 1000, quantileCutoff = 0.975){
  positive_cutoffs <- c()
  negative_cutoffs <- c()
  
  for(i in 1:iterations){
    rnd_smp <- sample(colnames(exprData), probesPerIter)
    
    cor_mat <- cor(exprData[, ..rnd_smp])
    corrs <- flattenCorrMatrix(cor_mat)
    
    positive_cutoffs <- c(positive_cutoffs, quantile(corrs$cor[corrs$cor > 0], quantileCutoff))
    negative_cutoffs <- c(negative_cutoffs, -1*quantile(-1*corrs$cor[corrs$cor < 0], quantileCutoff))
  }
  
  return(list(positiveCut = median(positive_cutoffs),
              negativeCut = median(negative_cutoffs)))
}

nullDiffCorr <- function(exprData, phenotype, iterations = 100, probesPerIter = 1000, quantileCutoff = 0.975){
  exprData$phenotype <- phenotype
  
  diffCorrCutoffs <- c()
  
  case <- exprData[phenotype == 1,]
  control <- exprData[phenotype == 0,]
  
  for(diffcorr_iter in seq(1, iterations)){
    
    rndSample <- sample(colnames(exprData)[seq(1,(ncol(exprData)-1))], probesPerIter)
    
    caseSelected <- case[, ..rndSample]
    controlSelected <- control[, ..rndSample]
    
    caseCorr <- cor(caseSelected)
    controlCorr <- cor(controlSelected)
    
    caseCorrFlat <- flattenCorrMatrix(caseCorr)
    controlCorrFlat <- flattenCorrMatrix(controlCorr)
    
    merged <- merge(caseCorrFlat, controlCorrFlat, by = c("row", "column"))
    merged$diffCorr <- abs(merged$cor.x - merged$cor.y)
    
    diffCorrCutoffs <- c(diffCorrCutoffs, quantile(merged$diffCorr, quantileCutoff, na.rm = TRUE))
  }
  
  return(median(diffCorrCutoffs))
}

pathwayEnrichment <- function(interestGene,
                              organism = "hsapiens", enrichDatabase="pathway_KEGG",
                              interestGeneType="ensembl_gene_id", referenceGeneType = "ensembl_gene_id",
                              referenceSet = "genome", isOutput = FALSE, ...){
  
  enrichment_results <- tryCatch(as.data.table(WebGestaltR(interestGene = interestGene,
                                                           organism = organism,
                                                           enrichDatabase = enrichDatabase,
                                                           interestGeneType = interestGeneType, 
                                                           referenceSet = referenceSet,
                                                           referenceGeneType = referenceGeneType, 
                                                           isOutput = FALSE, ...)),
                                 error = function(c){
                                   print("Error in enrichment. Skipping...")
                                   NULL
                                 })
}

stabilizedLasso <- function(exprData, phenotype, formula = NULL, 
                            seed = 2222,
                            iterations = 1000, splitRatio = 0.7, 
                            cvFolds = 10, parallel = FALSE,
                            propCutoff = 0.8, 
                            nullCorr, 
                            nullDiffCorr,
                            enrichmentOrganism = "hsapiens", 
                            enrichmentEnrichDatabase="pathway_KEGG",
                            enrichmentInterestGeneType="ensembl_gene_id", 
                            enrichmentReferenceGeneType = "ensembl_gene_id",
                            enrichmentReferenceSet = "genome", 
                            enrichmentIsOutput = FALSE){
  
  exprData$phenotype <- phenotype
  
  lassoResults <- c()
  # pred_metrics <- c()
  
  if(is.null(formula)){
    formula <- "phenotype ~ ."
  }
  
  set.seed(seed)
  for(k in seq(1,iterations)){
    
    # Splitting data proportional to phenotypeotype composition
    inds <- sample.split(exprData$phenotype, SplitRatio = splitRatio)
    
    #TODO: Add progressbar
    print(paste0("Currently fitting iteration: ", k, " of ", iterations))
    
    # Fitting a glmnet LASSO model for the training data
    cvfit <- glmnetUtils::cv.glmnet(formula = as.formula(formula), data = exprData[inds,], alpha = 1,
                                    family = "binomial",  
                                    nfolds = cvFolds, intercept = FALSE,
                                    parallel = parallel)
    
    # Obtaining fitted coefficients
    tmpCoeffs <- coef(cvfit)
    
    # Keeping track of results from each of the LASSO models fit in each iteration
    if(is.null(lassoResults)){
      lassoResults <- data.table(probe = tmpCoeffs@Dimnames[[1]][tmpCoeffs@i + 1], 
                            coefficient_1 = tmpCoeffs@x)
    }else{
      lassoResults <- merge(lassoResults, data.table(probe = tmpCoeffs@Dimnames[[1]][tmpCoeffs@i + 1], 
                                           coefficient = tmpCoeffs@x), 
                       all.x=TRUE, all.y=TRUE, 
                       by="probe", suffixes = c("",paste0("_", k)))
    }
    
  }
  
  coefficientSummary <- data.table(probe = lassoResults$probe,
                                   model.count = iterations - rowSums(is.na(lassoResults[, 2:ncol(lassoResults)])),
                                   min.val = apply(lassoResults[, 2:ncol(lassoResults)], 1, min, na.rm = TRUE),
                                   max.val = apply(lassoResults[, 2:ncol(lassoResults)], 1, max, na.rm = TRUE))
  
  coefficientSummary <- coefficientSummary[order(coefficientSummary$model.count, decreasing = TRUE)]
  
  identifiedGenes <- coefficientSummary$probe[(coefficientSummary$model.count/iterations) > propCutoff]
  identifiedGenes <- unique(unlist(str_split(identifiedGenes, ":")))
  
  corrResults <- c()

  caseData <- exprData[phenotype == 1,]
  controlData <- exprData[phenotype == 0,]
  
  print(paste0("Currently calculating CGNs"))
  if(length(identifiedGenes)>0){
    
    # For each core gene, co-expression values in the complete dataset,
    # and each of the phenotypes are calculated.
    for(i in seq(1, length(identifiedGenes))){
      
      coreGene <- identifiedGenes[i]
      gene1 <- exprData[[coreGene]]
      gene1Case <- caseData[[coreGene]]
      gene1Control <- controlData[[coreGene]]
      
      for(j in seq(1, (ncol(exprData)-1))){
        
        currGene <- colnames(exprData)[j]
        gene2 <- exprData[[currGene]]
        gene2Case <- caseData[[currGene]]
        gene2Control <- controlData[[currGene]]
        
        corTest <- cor.test(gene1, gene2, method = "pearson")
        caseTest <- cor.test(gene1Case, gene2Case, method = "pearson")
        controlTest <- cor.test(gene1Control, gene2Control, method = "pearson")
        
        currResult <- data.table(core_gene = coreGene,
                                 secondary_gene = currGene,
                                 correlation = corTest$estimate,
                                 pval = corTest$p.value,
                                 case.corr = caseTest$estimate,
                                 case.pval = caseTest$p.value,
                                 control.corr = controlTest$estimate,
                                 control.pval = controlTest$p.value)
        
        corrResults <- rbindlist(list(corrResults, currResult))
      }
    }
    
    corrResults$corrdiff <- abs(corrResults$case.corr - corrResults$control.corr)
    corrResults <- corrResults[corrResults$core_gene != corrResults$secondary_gene, ]
    
    colnames(corrResults) <- c("Core Gene", "Secondary Gene", 
                                "Correlation", "Correlation p-value", "Case Correlation", "Case Correlation p-value",
                                "Control Correlation", "Control Correlation p-value", "Differential Correlation")
    
    # Identifying CGNs using secondary genes that are co-expressed or 
    # differentially co-expressed.
    filteredResults <- corrResults[corrResults$Correlation > nullCorr$positiveCut |
                                     corrResults$Correlation < nullCorr$negativeCut |
                                     corrResults$`Differential Correlation` > nullDiffCorr,]
    
    pathwayResults <- c()
    
    print(paste0("Currently conducting SLASSO pathway enrichment"))
    for(currCoreGene in unique(filteredResults$`Core Gene`)){
      
      currGeneList <- c(currCoreGene, filteredResults$`Secondary Gene`[filteredResults$`Core Gene` == currCoreGene])
      currGeneList <- gsub("\\..*","", currGeneList)
      
      enrichmentResults <- pathwayEnrichment(interestGene = currGeneList,
                                             organism = enrichmentOrganism, 
                                             enrichDatabase = enrichmentEnrichDatabase,
                                             interestGeneType = enrichmentInterestGeneType, 
                                             referenceGeneType = enrichmentReferenceGeneType,
                                             referenceSet = enrichmentReferenceSet, 
                                             isOutput = enrichmentIsOutput)
      
      if(is.null(enrichmentResults) || nrow(enrichmentResults) == 0){
        print("No pathways enriched. Skipping...")
        next
      } 
      
      enrichmentResults$`Core Gene` <- currCoreGene
      
      pathwayResults <- rbindlist(list(pathwayResults, enrichmentResults))
    }
  }
  
  return(list(coefficientSummary = coefficientSummary,
              coreGenes = identifiedGenes,
              pathwayResults = pathwayResults))
}


