##### Byaesian Ridge regression model  #########
##load the required libarary and packages
library(bigmemory)
library(biganalytics)
library(BGLR)
library(rrBLUP)

##set the working directory
getwd()
workdir<- "C://Users//windows//OneDrive//Desktop//GenoOPT"
setwd <- workdir

## load the pheotyoic data and save 
Y <- read.csv("phenotypedata.csv", head = TRUE,  row.names = 1)

## load the genotypic data and save
X <- read.csv("genotypedata.csv", head = TRUE,  row.names = 1, stringsAsFactors=FALSE)

## callcutae Genomic Relationship matrix or Kinship Matrix
molecular_data <- X[, -1]
GRM <- A.mat(molecular_data)

###Fitting Bayesian Ridge Regression (BRR)

#Specify linear predictor
EtaR<-list(markers=list(X=X,model="BRR"))

##selecting the phenotype
y <- Y[,2]

set.seed(456)

#Fit the model
fmR<- BGLR(y=y,ETA=EtaR,nIter=10000,burnIn=5000,
          thin=10,verbose=FALSE)
betaHat<-fmR$ETA$markers$b
head(betaHat)

# Initialize a new plot
plot.new()

# Plot the data estimated marker effetcs of BRR
plot(betaHat, xlab="Marker", ylab="Estimated marker effect")
abline(h=0, col="red", lwd=2)
