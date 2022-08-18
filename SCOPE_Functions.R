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

nullCorr <- function(exprData, iterations = 100, probesPerIter = 1000, 
                     quantileCutoff = 0.975){
  positive_cutoffs <- c()
  negative_cutoffs <- c()
  
  print("Calculating null distribution for coexpression", quote = FALSE)
  pb <- txtProgressBar(min = 0, 
                       max = iterations, 
                       initial = 0,
                       style = 3)
  
  for(i in 1:iterations){
    setTxtProgressBar(pb, i)
    
    rnd_smp <- sample(colnames(exprData), probesPerIter)
    
    cor_mat <- cor(exprData[, ..rnd_smp])
    corrs <- flattenCorrMatrix(cor_mat)
    
    positive_cutoffs <- c(positive_cutoffs, quantile(corrs$cor[corrs$cor > 0], quantileCutoff))
    negative_cutoffs <- c(negative_cutoffs, -1*quantile(-1*corrs$cor[corrs$cor < 0], quantileCutoff))
  }
  
  close(pb)
  
  return(list(positiveCut = median(positive_cutoffs),
              negativeCut = median(negative_cutoffs)))
}

nullDiffCorr <- function(exprData, phenotype, iterations = 100, probesPerIter = 1000, 
                         quantileCutoff = 0.975){
  exprData$phenotype <- phenotype
  
  diffCorrCutoffs <- c()
  
  case <- exprData[phenotype == 1,]
  control <- exprData[phenotype == 0,]
  
  print("Calculating null distribution for differential coexpression", quote = FALSE)
  pb <- txtProgressBar(min = 0, 
                       max = iterations, 
                       initial = 0,
                       style = 3)
  
  for(diffcorr_iter in seq(1, iterations)){
    setTxtProgressBar(pb, diffcorr_iter)
    
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
  
  close(pb)
  
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

stabilizedLASSO <- function(exprData, phenotype, formula = NULL, 
                            seed = 2222,
                            iterations = 1000, splitRatio = 0.7,
                            cvFolds = 10, parallel = FALSE, family = "binomial",
                            propCutoff = 0.9,
                            ...){
  
  exprData$phenotype <- phenotype
  
  print("Running Stabilized LASSO", quote = FALSE)
  pb <- txtProgressBar(min = 0, 
                       max = iterations, 
                       initial = 0,
                       style = 3)
  
  lassoResults <- c()
  # pred_metrics <- c()
  
  if(is.null(formula)){
    formula <- "phenotype ~ ."
  }
  
  set.seed(seed)
  for(k in seq(1,iterations)){
    
    setTxtProgressBar(pb, k)
    
    # Splitting data proportional to phenotypeotype composition
    inds <- sample.split(exprData$phenotype, SplitRatio = splitRatio)
    
    #TODO: Add progressbar
    # print(paste0("Currently fitting iteration: ", k, " of ", iterations))
    
    # Fitting a glmnet LASSO model for the training data
    cvfit <- glmnetUtils::cv.glmnet(formula = as.formula(formula), data = exprData[inds,], alpha = 1,
                                    family = family,  
                                    nfolds = cvFolds, intercept = FALSE,
                                    parallel = parallel,
                                    ...)
    
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
  
  close(pb)
  
  coefficientSummary <- data.table(probe = lassoResults$probe,
                                   model.count = iterations - rowSums(is.na(lassoResults[, 2:ncol(lassoResults)])),
                                   min.val = apply(lassoResults[, 2:ncol(lassoResults)], 1, min, na.rm = TRUE),
                                   max.val = apply(lassoResults[, 2:ncol(lassoResults)], 1, max, na.rm = TRUE))
  
  coefficientSummary <- coefficientSummary[order(coefficientSummary$model.count, decreasing = TRUE)]
  
  identifiedGenes <- coefficientSummary$probe[coefficientSummary$model.count > (propCutoff*iterations)]
  
  return(list(coefficientSummary = coefficientSummary, 
              identifiedGenes = identifiedGenes))
  
}

identifyCoreGeneNetworks <- function(identifiedGenes,
                             exprData, 
                             phenotype,
                             nullCorr,
                             nullDiffCorr,
                             corrGenes = TRUE,
                             diffCorrGenes = TRUE){
  corrResults <- c()
  
  if(!corrGenes & !diffCorrGenes){
    stop("No CGNs can result from excluding BOTH correlated and differentially correlated genes.")
  }
  
  if(diffCorrGenes){
    caseData <- exprData[phenotype == 1,]
    controlData <- exprData[phenotype == 0,]
  }

  print("Identifying Core Gene Networks", quote = FALSE)
  pb <- txtProgressBar(min = 0, 
                       max = (length(identifiedGenes)*ncol(exprData)), 
                       initial = 0,
                       style = 3)
  
  prog <- 0
  for(i in seq(1, length(identifiedGenes))){
  
    coreGene <- identifiedGenes[i]
    
    if(corrGenes){
      gene1 <- exprData[[coreGene]]
    }

    if(diffCorrGenes){
      gene1Case <- caseData[[coreGene]]
      gene1Control <- controlData[[coreGene]]
    }
    
    for(j in seq(1, ncol(exprData))){
      prog <- prog + 1
      
      setTxtProgressBar(pb, prog)
      
      currGene <- colnames(exprData)[j]
      
      corTest <- list()
      corTest$estimate <- NA
      corTest$p.value <- NA
      
      caseTest <- list()
      caseTest$estimate <- NA
      caseTest$p.value <- NA
      controlTest <- list()
      controlTest$estimate <- NA
      controlTest$p.value <- NA
      
      if(corrGenes){
        gene2 <- exprData[[currGene]]
        corTest <- cor.test(gene1, gene2, method = "pearson")
      }
      
      if(diffCorrGenes){
        gene2Case <- caseData[[currGene]]
        gene2Control <- controlData[[currGene]]
        caseTest <- cor.test(gene1Case, gene2Case, method = "pearson")
        controlTest <- cor.test(gene1Control, gene2Control, method = "pearson")
      }

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
  
  close(pb)
  
  corrResults$corrdiff <- abs(corrResults$case.corr - corrResults$control.corr)
  corrResults <- corrResults[corrResults$core_gene != corrResults$secondary_gene, ]
  
  colnames(corrResults) <- c("Core Gene", "Secondary Gene", 
                             "Correlation", "Correlation p-value", "Case Correlation", "Case Correlation p-value",
                             "Control Correlation", "Control Correlation p-value", "Differential Correlation")
  
  filteredResults <- c()
  
  # Identifying CGNs using secondary genes that are co-expressed or 
  # differentially co-expressed.
  if(corrGenes){
    filteredResults <- rbindlist(list(filteredResults, 
                                      corrResults[corrResults$Correlation > nullCorr$positiveCut |
                                                                     corrResults$Correlation < nullCorr$negativeCut,]))
  }
  
  if(diffCorrGenes){
    filteredResults <- rbindlist(list(filteredResults, 
                                  corrResults[corrResults$`Differential Correlation` > nullDiffCorr,]))
  }
  
  CGNs <- data.table(`Core Gene` = unique(filteredResults$`Core Gene`))
  CGNs$`Secondary Genes` <- sapply(CGNs$`Core Gene`, function(x) {paste0(filteredResults$`Secondary Gene`[filteredResults$`Core Gene`==x], collapse = ";")})

  return(list(CGNs = CGNs, fullCGNs = filteredResults))
}

enrichCGNs <- function(coreGeneNetworks,
                       enrichmentOrganism = "hsapiens", 
                       enrichmentEnrichDatabase="pathway_KEGG",
                       enrichmentInterestGeneType="ensembl_gene_id", 
                       enrichmentReferenceGeneType = "ensembl_gene_id",
                       enrichmentReferenceSet = "genome", 
                       enrichmentIsOutput = FALSE,
                       ...){
  
  pathwayResults <- c()
  
  print("Conducting pathway enrichment using WebGestaltR", quote = FALSE)
  
  for(currCoreGene in unique(coreGeneNetworks$CGNs$`Core Gene`)){
    
    currGeneList <- c(currCoreGene, str_split(coreGeneNetworks$CGNs$`Secondary Genes`[coreGeneNetworks$CGNs$`Core Gene` == currCoreGene], ";", simplify = TRUE))
    
    enrichmentResults <- pathwayEnrichment(interestGene = currGeneList,
                                           organism = enrichmentOrganism, 
                                           enrichDatabase = enrichmentEnrichDatabase,
                                           interestGeneType = enrichmentInterestGeneType, 
                                           referenceGeneType = enrichmentReferenceGeneType,
                                           referenceSet = enrichmentReferenceSet, 
                                           isOutput = enrichmentIsOutput,
                                           ...)
    
    if(is.null(enrichmentResults) || nrow(enrichmentResults) == 0){
      print("No pathways enriched. Skipping...")
      next
    } 
    
    enrichmentResults <- data.table(`Core Gene` = currCoreGene,
                                    enrichmentResults)
    
    pathwayResults <- rbindlist(list(pathwayResults, enrichmentResults))
  }
  
  return(pathwayResults)
}

pathwayOverlaps <- function(enrichedCGNs){
  
  print("Calculating pathway overlaps", quote = FALSE)
  overlappedResults <- unique(enrichedCGNs[, c("geneSet", "description")])
  for(currCoreGene in unique(enrichedCGNs$`Core Gene`)){
    overlappedResults[[currCoreGene]] <- sapply(overlappedResults$geneSet, function(x) {
      if(nrow(enrichedCGNs[`Core Gene` == currCoreGene & geneSet == x]) == 0){
        return(0)
      }else{
        return(enrichedCGNs[`Core Gene` == currCoreGene & geneSet == x, overlap] / enrichedCGNs[`Core Gene` == currCoreGene & geneSet == x, size])
      }
      }
      )
  }
  
  overlappedResults$POS <- rowSums(as.matrix(overlappedResults[ , -c(1,2)]))
  overlappedResults <- overlappedResults[order(POS, decreasing = TRUE)]
  
  return(overlappedResults)
  
}

scope <- function(exprData, phenotype,
                  removeZeroVar = TRUE,
                  removeLowVar = TRUE,
                  lowVarPercentile = 0.25,
                  formula = NULL, 
                  seed = 2222,
                  iterations = 1000, splitRatio = 0.7, 
                  cvFolds = 10, parallel = FALSE, family = "binomial",
                  propCutoff = 0.8,
                  corrIterations = 100, corrProbesPerIter = 1000, 
                  corrQuantileCutoff = 0.975,
                  diffCorrIterations = 100, diffCorrProbesPerIter = 1000, 
                  diffCorrQuantileCutoff = 0.975,
                  corrGenesCGN = TRUE,
                  diffCorrGenesCGN = TRUE,
                  enrichmentOrganism = "hsapiens",
                  enrichmentEnrichDatabase="pathway_KEGG",
                  enrichmentInterestGeneType="ensembl_gene_id",
                  enrichmentReferenceGeneType = "ensembl_gene_id",
                  enrichmentReferenceSet = "genome",
                  enrichmentIsOutput = FALSE,
                  ...){
  
  if(removeZeroVar) exprData <- removeZeroVariance(exprData = exprData)
  if(removeLowVar) exprData <- removeLowVariance(exprData = exprData, quantileToRemove = lowVarPercentile)
  
  slassoResults <- stabilizedLASSO(exprData = exprData,
                                   phenotype = phenotype,
                                   formula = formula,
                                   seed = seed,
                                   iterations = iterations,
                                   splitRatio = splitRatio,
                                   cvFolds = cvFolds,
                                   parallel = parallel,
                                   family = family,
                                   propCutoff = propCutoff)
  
  nullCorrelations <- nullCorr(exprData = exprData, 
                               iterations = corrIterations, 
                               probesPerIter = corrProbesPerIter, 
                               quantileCutoff = corrQuantileCutoff)
  
  nullDiffCorrelations <- nullDiffCorr(exprData = exprData,
                                       phenotype = phenotype, 
                                       iterations = diffCorrIterations, 
                                       probesPerIter = diffCorrProbesPerIter,
                                       quantileCutoff = diffCorrQuantileCutoff)
  
  coreGeneNetworks <- identifyCoreGeneNetworks(identifiedGenes = slassoResults$identifiedGenes,
                                               exprData = exprData,
                                               phenotype = phenotype, 
                                               nullCorr = nullCorrelations, 
                                               nullDiffCorr = nullDiffCorrelations,
                                               corrGenes = corrGenesCGN,
                                               diffCorrGenes = diffCorrGenesCGN)
  
  enrichedCGNs <- enrichCGNs(coreGeneNetworks = coreGeneNetworks)
  
  pathwayOverlapResults <- pathwayOverlaps(enrichedCGNs = enrichedCGNs)
  
  return(list(stabilizedLassoResults = slassoResults,
              nullCorrelations = nullCorrelations,
              nullDiffCorrelations = nullDiffCorrelations,
              coreGeneNetworks = coreGeneNetworks,
              enrichedCGNs = enrichedCGNs,
              pathwayOverlapResults = pathwayOverlapResults))
}


