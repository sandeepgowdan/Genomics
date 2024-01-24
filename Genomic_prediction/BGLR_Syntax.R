library(BGLR)
data(wheat)

Y<- wheat.Y[,4] # Phenotypes
X<- wheat.X # Genotypes
dim(X)

length(Y)

X<-scale(wheat.X, center = TRUE, scale = TRUE) # Molecular data


ETA<-list(X=list(X=X,model="BRR")) # To run BayesB -> model="BayesB"

fm<-BGLR(y=Y,ETA=ETA,nIter=10000,burnIn=2000, thin = 5, verbose = FALSE)

accur1<-cor(fm$yHat, Y)

