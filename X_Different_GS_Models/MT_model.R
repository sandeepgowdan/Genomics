### Multi Trait Model (MTM)
##install.packages("MTM")
install.packages("remotes")
remotes::install_github("QuantGen/MTM")
library(MTM)
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

#Linear predictor
EtaM<-list(list(K = G, COV = list(type = "UN",df0 = 4,S0 =
                            diag(3))) )
#Residual
residual<-list(type = "UN",S0 = diag(3),df0 = 4)
fmM <- MTM(Y = Y, K=EtaM,resCov=residual,
              nIter=10000,burnIn=5000,thin=10)

#Predictions of phenotypical values
fmM$YHat

#Predictions of random effects
fmM$K[[1]]$U

#Residual covariance matrix
fmM$resCov$R

#Genetic covariance matrix
fmM$K[[1]]$G



