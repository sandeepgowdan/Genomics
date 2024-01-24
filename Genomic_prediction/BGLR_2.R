## Deriving Eigenvectors of the phenotype matrix

rm(list = ls())

library(BGLR)


# A phenotype matrix with 4 traits
Y2 <- read.csv("phenotype1.txt") 
Y2 <- as.numeric(Y2)  # Convert Y2 to numeric

# Check if the conversion was successful
if (!is.na(Y2)) {
  Y2 <- scale(Y2)  # Scale Y2
} else {
  # Handle the case when the conversion to numeric fails
  # Print an error message or take appropriate action
  print("Y2 contains non-numeric values.")
}

Y2 <- as.matrix(Y2)

Y2 <- scale(Y2)

# Singular value decompositions Y=UDV' Singular Value Decomposition of a Matrix

SVD = svd(Y2)
U=SVD$u
D=diag(SVD$d)
V=SVD$v


# Molecular markers
X <- read.csv("numeric2_geno_r1g1_male.csv", check.names = F, header = T)

X <- scale(X)

## regression coefficients

B=matrix(nrow=ncol(X),ncol=ncol(Y2))

ETA=list(list(X=X,model='BayesB'))

for(i in 1:ncol(Y2)){
  fm=BGLR(y=U[,i],ETA=ETA,verbose=F, nIter = 12000, burnIn = 2000) #use more iterations!
  B[,i]=fm$ETA[[1]]$b
}

# Rotating coefficients to put them in marker space
BETA=B%*%D%*%t(SVD$v)

# Prediction
YHat=X%*%BETA

# correlation in training
for(i in 1:ncol(Y2)){ print(cor(Y2[,i],YHat[,i])) }


#Training-testing comparing single and multi-trait

## Training/Testing partition

set.seed(195021)
N=nrow(X)
tst=sample(1:N,size=300)
Y.TRN=Y2[-tst,]
Y.TST=Y2[tst,]

X.TST=X[tst,]

ETA=list(list(X=X[-tst,],model='BayesB'))

## Prediction using single-trait models

YHAT.ST=matrix(nrow=nrow(Y.TST),ncol=5)

for(i in 1:5){
  fm=BGLR(y=Y.TRN[,i],ETA=ETA,verbose=F, nIter = 12000, burnIn = 2000) #use more iterations!
  YHAT.ST[,i]=X.TST%*%fm$ETA[[1]]$b
}	

## Prediction using eigenvectors aka MT models ??

SVD=svd(Y.TRN)
U=SVD$u
D=diag(SVD$d)
V=SVD$v

B=matrix(nrow=ncol(X), ncol=ncol(U))

for(i in 1:5){
  fm=BGLR(y=U[,i],ETA=ETA,verbose=F , nIter = 12000, burnIn = 2000) #use more iterations!
  B[,i]=fm$ETA[[1]]$b
}	

BETA=B%*%D%*%t(V)

YHAT.EV=X.TST%*%BETA

## Comparison

for(i in 1:5){
  message( round(cor(Y.TST[,i],YHAT.ST[,i]),3),"   ",round(cor(Y.TST[,i],YHAT.EV[,i]),3))
}


TABLE: Independent predictions with uni- and multitrait genomic selection models
