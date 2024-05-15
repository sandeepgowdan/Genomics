# Load required libraries
library(bigmemory)
library(biganalytics)
library(BGLR)
library(rrBLUP)

# Set the working directory
workdir <- "C://Users//windows//OneDrive//Desktop//GenoOPT"
setwd(workdir)

# Load the phenotypic data
Y <- read.csv("phenotypedata.csv", head = TRUE, row.names = 1)

# Load the genotypic data
X <- read.csv("genotypedata.csv", head = TRUE, row.names = 1, stringsAsFactors=FALSE)

# Calculate the Genomic Relationship Matrix (GRM)
molecular_data <- as.matrix(X[, -1])  # Convert to matrix if necessary
GRM <- A.mat(molecular_data)

# Select the phenotype
y <- Y[, 2]

# Create folds for cross-validation
set.seed(456)  # Set seed for reproducibility
folds <- sample(rep(1:5, length.out = length(y)))

# Initialize a vector to store Mean Squared Errors (MSE) for each fold
mse_list <- numeric(5)

# Perform 5-fold cross-validation
for (k in 1:5) {
  # Split the data into training and test sets
  train_indices <- which(folds != k)
  test_indices <- which(folds == k)
  
  y_train <- y[train_indices]
  y_test <- y[test_indices]
  
  X_train <- molecular_data[train_indices, ]
  X_test <- molecular_data[test_indices, ]
  
  # Specify linear predictor for training data
  EtaR_train <- list(markers = list(X = X_train, model = "BRR"))
  
  # Fit the model on training data
  fmR <- BGLR(y = y_train, ETA = EtaR_train, nIter = 10000, burnIn = 5000, thin = 10, verbose = FALSE)
  
  # Predict on test data
  y_pred <- predict(fmR, ETA = list(markers = list(X = X_test)))
  
  # Calculate Mean Squared Error (MSE) for the current fold
  mse_list[k] <- mean((y_test - y_pred)^2)
}

# Calculate the average MSE across all folds
avg_mse <- mean(mse_list)
cat("Average Mean Squared Error (MSE) across 5 folds:", avg_mse, "\n")

# Plot the estimated marker effects of BRR from the last fold model
betaHat <- fmR$ETA$markers$b
plot(betaHat, xlab = "Marker", ylab = "Estimated marker effect")
abline(h = 0, col = "red", lwd = 2)
