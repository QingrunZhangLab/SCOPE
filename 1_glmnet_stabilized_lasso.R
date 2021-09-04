library(data.table)
library(caret)
library(caTools)
library(glmnetUtils)
library(doParallel)
library(tictoc)

# Specify number of threads to use for each cross-validation fold in glmnet package.
# Please note that the higher, the number, the larger the memory usage and that there 
# is a time cost associated with copying the data needed for each thread. So even if
# sufficient RAM is available, a lower number of threads may be ideal.
registerDoParallel(2)

args = commandArgs(trailingOnly=TRUE)

# Default values for parameters
ITERS = 1000
FOLDS = 10
START_SEED = 2222
PREFIX = "stabilized_lasso"

if(length(args)==0) {
  stop("Please specify at minimum the expression data to be processed. Script can be
       run with defaults for other arguments.")
}else if(length(args)==1){
  print("Since no optional arguments were specified, default values will be used.")
}else{
  ITERS = as.numeric(args[2])
  FOLDS = as.numeric(args[3])
  START_SEED = as.numeric(args[4])
  PREFIX = args[5]
}

if(!file.exists(args[1])){
  stop(paste0("Input file ", args[1], " does not exist! Please check the file path and retry."))
}

print(paste0("Loading file: ", args[1]))
expr.data <- fread(args[1])
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

# Selecting only Primary Tumor and Solid Tissue Normal samples
print("Filtering to exclude unwanted phenotypes...")
expr.data <- expr.data[expr.data$phen %in% c("Primary Tumor", "Solid Tissue Normal")]
expr.data$phen  <- ifelse(expr.data$phen == "Solid Tissue Normal", 0, 1)
expr.data <- expr.data[, -c(1)]

results <- c()
pred_metrics <- c()
set.seed(START_SEED)

for(k in 1:ITERS){
  
  # Splitting data proportional to phenotype composition
  inds  <- sample.split(expr.data$phen, SplitRatio = 0.7)
  
  print(paste0("Currently fitting iteration: ", k))
  
  tic("Iteration time")
  # Fitting a glmnet LASSO model for the training data
  cvfit <- glmnetUtils::cv.glmnet(phen ~ ., data = expr.data[inds,], alpha = 1,
                                  family = "binomial",  
                                  nfolds = FOLDS, intercept = FALSE,
                                  parallel = TRUE)
  
  # Obtaining fitted coefficients
  tmp_coeffs <- coef(cvfit)
  print("Fitted...")
  
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
  
  print("Predicting...")
  tokeep <- which(sapply(expr.data, is.numeric))
  test_data <- as.matrix(expr.data[!inds, ..tokeep])
  cvpred <- predict(cvfit, test_data)
  
  toc(quiet = FALSE)
  
  print("Calculate prediction metrics...")
  preds <- ifelse(cvpred < 0.5, 0, 1)
  acc <- mean(preds == expr.data$phen[!inds])
  sens <- sensitivity(as.factor(preds), as.factor(expr.data$phen[!inds]), positive = "1", negative = "0")
  spef <- specificity(as.factor(preds), as.factor(expr.data$phen[!inds]), positive = "1", negative = "0")
  
  pred_metrics <- rbindlist(list(pred_metrics, 
                                 data.table(model_num = k,
                                            Accuracy = acc,
                                            Sensitivity = sens,
                                            Specificity = spef)))
}

print("Combining LASSO results...")
cof_sum <- data.table(probe = results$probe,
                      model.count = ITERS - rowSums(is.na(results[,2:ncol(results)])),
                      min.val = apply(results[,2:ncol(results)], 1, min, na.rm = TRUE),
                      max.val = apply(results[,2:ncol(results)], 1, max, na.rm = TRUE))

print("Writing all results.")
fwrite(cof_sum, paste0(PREFIX, "_Summary.csv"))
fwrite(results, paste0(PREFIX, "_AllCoefs.csv"))
fwrite(pred_metrics, paste0(PREFIX, "_PredMetrics.csv"))