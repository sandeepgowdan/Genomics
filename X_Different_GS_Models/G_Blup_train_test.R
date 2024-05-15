### train test split on GBLUP model

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

# Scale genotype data
Z <- scale(molecular_data, center = TRUE, scale = TRUE)
G <- tcrossprod(Z) / ncol(Z)

# Specify linear predictor
EtaG <- list(markers = list(K = G, model = "RKHS"))

# Train-test split
set.seed(123)  # Set seed for reproducibility
train_indices <- sample(1:nrow(Y), size = 0.8 * nrow(Y))  # 80% for training
test_indices <- setdiff(1:nrow(Y), train_indices)  # Remaining 20% for testing

Y_train <- Y[train_indices, ]
Y_test <- Y[test_indices, ]

# Corresponding genotype data for training and test sets
X_train <- X[train_indices, ]
X_test <- X[test_indices, ]

# Fit the model on training data
y_train <- Y_train[, 2]  # Assuming the phenotype is in the second column
set.seed(789)
fmG <- BGLR(y = y_train, ETA = EtaG, nIter = 10000, burnIn = 5000, thin = 10, verbose = FALSE)

# BLUPs
uHat <- fmG$ETA$markers$u
uHat

# Variance parameters
varE <- fmG$varE  # Residual
varU <- fmG$ETA$markers$varU  # Genotypes

# Make predictions for the test set
# Calculate the predictions using the BLUPs
Z_test <- scale(X_test[, -1], center = TRUE, scale = TRUE)
predictions <- Z_test %*% uHat

# Evaluate the model (e.g., calculate Mean Squared Error)
y_test <- Y_test[, 2]  # Actual test phenotypes
MSE <- mean((y_test - predictions)^2)
cat("Mean Squared Error on test set:", MSE, "\n")

# Optional: Plot the predictions vs actual values
plot(y_test, predictions, xlab = "Actual Phenotype", ylab = "Predicted Phenotype", main = "GBLUP Predictions vs Actual")
abline(0, 1, col = "red", lwd = 2)

