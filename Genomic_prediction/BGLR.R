install.packages("bigmemory")
install.packages("biganlytics")
install.packages("BGLR")
library(bigmemory)
library(biganalytics)
library(BGLR)

getwd()
Y <- read.table("phenotype.txt", head = TRUE)
y <- Y[,2]
X <- read.table("genotype.txt", head = TRUE)
A <- read.table("Kinship.txt", head = TRUE)

dim(y)
length(y)
dim(X)
dim(A)

#Computing the genomic relationship matrix

X <- scale(X, center = TRUE, scale = TRUE)

# Calculate terossprod
G <- tcrossprod(X) / ncol(X)


#Computing the eigen-value decomposition of G
EVD <-eigen(G)

#Setting the linear predictor

ETA<-list(list(K=A, model='RKHS'),
          
          list(V=EVD$Svectors,d=EVD$Svalues, model='RKHS'))

#Fitting the model

ETA<- (list(list(K=G, model = 'RKHS')))
fm <- BGLR(y = y, ETA = ETA, nIter = 12000, burnIn = 2000, saveAt = "PGBLUP_")
save(fm,file='fmPG_BLUP.rda')



#Predictions

yHat<-fm$yHat
tmp<-range(c(y,yHat))

plot(yHat~y,xlab='Observed',ylab='Predicted',col=2, 
     xlim=tmp,ylim=tmp); abline(a=0,b=1,col=4,lwd=2)

#Exporting your Genomie prediction values 

write.table(yHat, "C:\\Users\\windows\\OneDrive\\Desktop\\file. st", sep="'/")

#Godness of fit and related statisties

fm$fit

fm$varE # compare to var(y)

#Variance components associated with the genomie and pedigree 

fm$ETA[[1]]$varU

fm$ETA[[2]]$varU

# Residual variance

varE<-scan('PGBLUP_varE.dat')
plot(varE,type='o', col= 2, cex= 0.5)






