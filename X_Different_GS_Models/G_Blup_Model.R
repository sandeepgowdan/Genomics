
## Fitting Genomic Best Linear Unbiased Predictor (GBLUP)
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

## Fitting Genomic Best Linear Unbiased Predictor (GBLUP)
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

#Specify linear predictor
EtaG<-list(markers=list(K=G,model="RKHS"))

set.seed(789)
#Fit the model
fmG<-BGLR(y=y,ETA=EtaG,nIter=10000,
             burnIn=5000,thin=10,
             verbose=FALSE)
##BLUPs
uHat<-fmG$ETA$markers$u
uHat


#Variance parameters
fmG$varE #Residual
fmG$ETA$markers$varU #Genotypes

