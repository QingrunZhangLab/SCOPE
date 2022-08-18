source("SCOPE_Functions.R")

library(doParallel)

doParallel::registerDoParallel(cores = 4)

expr.data <- fread("SampleData/SampleExpressionData.csv")
phen <- fread("SampleData/SamplePhenotypes.csv")[[1]]

set.seed(1234)

scopeOutputs <- scope(exprData = expr.data, 
                      phenotype = phen,
                      removeZeroVar = TRUE,
                      removeLowVar = TRUE,
                      lowVarPercentile = 0.25,
                      formula = NULL, 
                      seed = 2222,
                      iterations = 10, splitRatio = 0.7, 
                      cvFolds = 10, parallel = TRUE, family = "binomial",
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
                      enrichmentIsOutput = FALSE)