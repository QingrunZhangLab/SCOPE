library(data.table)

# This script will generate the null distribution for co-expression given
# the expression matrix.

# Function to flatten the correlation matrix
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}

args = commandArgs(trailingOnly=TRUE)

# Default values for parameters
ITERS = 100
PROBES_PER_ITER = 1000
START_SEED = 2222
PREFIX = "stabilized_lasso"

if(length(args)==0) {
  stop("Please specify at minimum the expression data to be processed. Script can be
       run with defaults for other arguments.")
}else if(length(args)==1){
  print("Since no optional arguments were specified, default values will be used.")
}else{
  ITERS = as.numeric(args[2])
  PROBES_PER_ITER = as.numeric(args[3])
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

# Taking a look at available phenotypes
print("Current phenotypes: ")
print(table(expr.data$phen))

# Selecting only Primary Tumor and Solid Tissue Normal samples
print("Filtering to include only 'Primary Tumor' and 'Solid Tissue Normal' samples.")
expr.data <- expr.data[expr.data$phen %in% c("Primary Tumor", "Solid Tissue Normal")]

null_dist <- c()

set.seed(START_SEED)

# In each iteration 'PROBES_PER_ITER' number of columns are selected
# from the expression matrix and the correlations among these 
# calculated. Correlations are then stored to be later used for
# determining the null distribution.
for(i in 1:ITERS){
  print(paste0("Current trial ", i))
  
  rnd_smp <- sample(colnames(expr.data)[2:(ncol(expr.data)-1)], PROBES_PER_ITER)
  
  cor_mat <- cor(expr.data[, ..rnd_smp])
  corrs <- flattenCorrMatrix(cor_mat)
  
  if(is.null(null_dist)){
    null_dist <- as.data.table(corrs$cor)
  }else{
    null_dist <- cbind(null_dist, corrs$cor)
  }
}

print("Writing all results.")
fwrite(null_dist, paste0(PREFIX, "_CorrNull.csv"))