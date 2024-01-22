## INSTALLING PACKAGES AND LIBRARY

install.packages("BGLR",repos="https://cran.r-project.org/")
install.packages("glmnet",repos="https://cran.r-project.org/")
install.packages("AGHmatrix", repos="https://cran.r-project.org/")
install.packages('learnr',repos="https://cran.r-project.org/")
install.packages('downloadthis',repos="https://cran.r-project.org/")
install.packages("remotes",repos="https://cran.r-project.org/")
remotes::install_github("rstudio/gradethis")
install.packages("tensorflow",repos="https://cran.r-project.org/")
install.packages("keras",repos="https://cran.r-project.org/")
install.packages("tfdatasets", repos="https://cran.r-project.org/")

# Install the package
install.packages("genetics")

# Load the package
library(genetics)
library(ggplot)


## DATA CURATION
# this is a dataset included with BGLR package
data(wheat)
# marker data, genotypes are coded as 0/1 as they are inbred lines
X = wheat.X
# phenotypes
Y = wheat.Y
# N obs
n = nrow(X)
# N markers
p = ncol(X)
# phenotype
y = Y[,2]

## DATA INSPECTION
# PC is computed on the correlation matrix between individual genotypes
pc = eigen(tcrossprod(scale(X))/ncol(X))

# variance explained by first two components
v1=round(pc$values[1]/sum(pc$values),3)*100
v2=round(pc$values[2]/sum(pc$values),3)*100

#Check allele frequency
# assumed genotypes coded as 0/1/2
if (max(X)==3) {X=X-1} # if genotypes 1/2/3 --> 0/1/2
if (max(X)==1) {X=X*2} # if genotypes 0/1 --> 0/2
p = apply(X, 2, mean) * 0.5
# this takes minor allele frequency (p<0.5)
p[p>0.5] = p[p>0.5]-0.5
hist(p, main='Allele frequencies')

# ggplot
p <- as_tibble(p)
ggplot(p, aes(x=value)) + 
  geom_histogram(bins=10, fill='grey', colour='black') + 
  theme_minimal()+ 
  labs(title='Minimum allele frequency')
plot(pc$vectors[,1], pc$vectors[,2], xlab = paste0('PC1 % = ', v1), ylab=paste0('PC2 % = ',v2))

## DATA PARTITION
# test set comprises 20% of data, randomly chosen

tst = sort(sample(1:n, size=n*0.2, replace=FALSE))

# train set
XTRN<-X[-tst,]
yTRN<-y[-tst]

# testing test
XTST<-X[tst,]
yTST<-y[tst]

##Given that structuring is present, you may want to check whether the two groups are represented in both training and test partitions.
# Train set PCA
# I may need to remove fixed markers, those for which SD is zero
X1 <- XTRN[, which(apply(XTRN, 2, sd) != 0)]
pcTRN = eigen(tcrossprod(scale(X1))/ncol(X))
v1=pcTRN$values[1]/sum(pcTRN$values)
v2=pcTRN$values[2]/sum(pcTRN$values)
plot(pcTRN$vectors[,1], pcTRN$vectors[,2], xlab=v1, ylab=v2, main='Train partition')

# Test PCA
X2 <- XTST[, which(apply(XTST, 2, sd) != 0)]
pcTST = eigen(tcrossprod(scale(X2))/ncol(X))
v1=pcTST$values[1]/sum(pcTST$values)
v2=pcTST$values[2]/sum(pcTST$values)
points(pcTST$vectors[,1], pcTST$vectors[,2], col=2)


##GBLUP
# phenotypes with tst individuals missing values
yNA = y
yNA[tst] = NA

# assumed heritability
h2 = 0.4

# this function computes GRM with Van Raden algorithm
vanraden = function(X){
  # allele frequencies (assumes genotypes coded as 0,1,2)
  p = apply(X, 2, mean) * 0.5
  G = tcrossprod(scale(X, scale = FALSE))
  G = G / (2 * sum(p*(1-p)))
  diag(G) = diag(G)*1.05
  return(G)
}


# genomic relationship matrix
GRM = vanraden(X)

# this function returns left and right hand side of mixed model equations
mme = function(y, V, h2, invert=FALSE){
  # y = phenotypes; V = covariance matrix; h2 = heritability
  if (invert) V = solve(V)
  n = length(y)
  x = rep(1,n)
  # replace missing with 0's
  x[is.na(y)] = 0
  y[is.na(y)] = 0
  Z = diag(x)
  LHS = matrix(nrow=n+1, ncol=n+1)
  V = V *(1-h2)/h2 + t(Z)%*%Z
  LHS[1,] = c(t(x)%*%(x), x)
  LHS[-1,] = cbind(x,V)
  RHS = c(t(x)%*%y, t(Z%*%y))
  return(list('LHS'=LHS, 'RHS'=RHS))
}


# mixed model equations, evaluated only in training set
MME = mme(y=yNA, V=GRM, h2=h2, invert=TRUE)

# the predicted breeding values are the last n solutions
uhat = solve(MME$LHS, MME$RHS)[-1]

# corr between observed and predicted phenotypes
yhatG = uhat[tst]
print(c('corr obs net', round(cor(yhatG,y[tst]),3)))

# plot
plot(y[tst],yhatG, main='Obs vs GBLUP', xlab='Obs', ylab = 'GBLUP')



# plot
ds <- data.frame( "true"=yTST, "predicted" = yhatG)
ggplot(ds, aes(x=true, y=predicted)) + 
  geom_point(col="red") + 
  stat_cor(method = "pearson",  aes(label = ..r.label..),label.x = -2.5, label.y = 1) +
  theme_minimal() +
  theme(axis.title.x = element_text("observed")) + 
  theme(axis.title.y = element_text("Gblup estimator"))

##BGLR
# removes phenotypes from TST partition, tst contains ids of testing set
yNA = y
yNA[tst] = NA

##Bayesian Ridge regression
### fit ridge model only to training set (yNA)
ETA<-list(list(X=X,model="BRR"))
fmRR<-BGLR(y=yNA,ETA=ETA,nIter=5000,burnIn=2000,verbose=FALSE)

# trace plot of the residual variance
varE<- scan("varE.dat")
plot(varE, ylab='Var E', xlab='Iteration', type="o",col=2,cex=.5)
abline(h=fmRR$varE,col=4, lwd=2)

## plots of estimates
plot(fmRR$ETA[[1]]$b,col=4,ylab='Estimate',main='BRR')

print(c('Corr Obs vs. BRR', round(cor(fmRR$yHat[tst], y[tst]),3)))

plot(y[tst], fmRR$yHat[tst],col=4,ylab='Observed y',xlab='BRR prediction')

##BAyes A
### BayesA(Scaled-t prior)
ETA<-list(list(X=X,model="BayesA"))
fmBA<-BGLR(y=yNA,ETA=ETA,nIter=5000,burnIn=2000,verbose=FALSE)

print(c('Corr Obs vs. Bayes A', round(cor(fmBA$yHat[tst], y[tst]),3)))

plot(y[tst], fmBA$yHat[tst],col=4, ylab='Observed y',xlab='Bayes A prediction')

### BayesB (point of mass at zero + scaled-t slab)
ETA<-list(list(X=X,model="BayesB"))
fmBB<-BGLR(y=yNA,ETA=ETA,nIter=5000,burnIn=2000,verbose=FALSE)

print(c('Corr Obs vs. Bayes B', round(cor(fmBB$yHat[tst], y[tst]),3)))

plot(y[tst], fmBB$yHat[tst],col=4, ylab='Observed y',xlab='Bayes B prediction')


### BayesC 
ETA=list(list(X=X,model='BayesC'))

fmBC = BGLR(y=yNA,ETA=ETA,nIter=5000,burnIn=2000,verbose=F)

### corr between predicted and observed
print(c('Corr Obs vs. Bayes C', round(cor(fmBC$yHat[tst], y[tst]),3)))

plot(y[tst], fmBC$yHat[tst],col=4,ylab='Observed y',xlab='Bayes C prediction')

# GRM
doGRM = function(X) {
  G = tcrossprod(scale(X))
  G = G/mean(diag(G))
  # recommended / required to avoid non positive definite matrices
  diag(G) = diag(G)*1.05
  return(G)
}

# Bayesian GBLUP, BRR with PC decomposition
G = doGRM(X)
EVD<-eigen(G)
PC = EVD$vectors%*%diag(sqrt(EVD$values))
ETA=list(list(X=PC,model='BRR'))

fmGB<-BGLR(y=yNA,ETA=ETA,nIter=5000,burnIn=2000,verbose=F)

### corr between predicted and observed
print(c('Corr Obs vs. Bayes GBLUP', round(cor(fmGB$yHat[tst], y[tst]),3)))

plot(y[tst], fmGB$yHat[tst],col=4,ylab='Observed y',xlab='BRR prediction')



##Single step GBLUP is one of the most important developments in genomic prediction, 
##since it allows considering both genotyped and ungenotyped individuals, 
##provided they are connected by a pedigree.

##SS is based on deriving a relationship matrix H that uses both marker and
##pedigree information. The most amazing property is that the inverse of H 
##is very easy to compute and requires to invert a matrix of size only the 
##number of genotyped individuals. This is quite convenient, since the number 
##of ungenotyped individuals istypically far larger than those with marker information

## Generating pedigree

# no of families
nf = 10

# cluster
h = hclust(dist(X),method="ward.D2")

# cuts the tree into nf clusters (families), hcut contains cluster id for each ind
hcut = cutree(h,nf)

# size is n plus two parents per family
N = nrow(X) + 2*nf
ped = matrix(0, nrow=N, ncol=3)
ped[,1] = seq(N)
ped[(2*nf+1):N,2] = hcut*2-1
ped[(2*nf+1):N,3] = hcut*2

# AGHmatrix package used to compute pedigree based relationship matrix
A = Amatrix(ped)


## Specifying genotyped individuala
# no. individuals in pedigree
N = nrow(A)
# no. individuals genotyped available
n = nrow(X)
# no. individuals not genotyped = founders + n/5
# %/% is the integer division
n0 = nf*2 + n %/% 5

# list of genotyped individuals (last n*4/5), ids refer to whole pedigree
idX = seq((n0 + 1), N)

# Xss: actual X to be used
# IMPORTANT: we need to subtract no. of founders since X contains genotypes of non-founders only (dimension n), idX refers to the global pedigree ids.
Xss = X[(idX-2*nf),]

# filter by non segregating SNPs
# not needed in vanraden option, but required for doGRM function
Xss = Xss[, which(apply(Xss, 2, sd) != 0)]


##H inverse
Gss = vanraden(Xss)
G_1 = solve(Gss)

# block in A for genotyped individuals, inverted
A22_1 = solve(A[idX,idX])

# inverse of A. Efficient rules to invert A are known since 1976, 
# here we simply use brute force, (sorry, Chuck).
A_1 = solve(A)

# H inverse
H_1 = A_1
H_1[idX,idX] = H_1[idX,idX] + G_1 - A22_1

# That's it!

##Direct H
H = matrix(nrow=N,ncol=N)
A22_1 = solve(A[idX,idX])
A12 = A[-idX,idX]

H11 = A[-idX,-idX] + A12 %*% A22_1 %*% (Gss-A[idX,idX]) %*% A22_1 %*% t(A12)
H12 = A12 %*% A22_1 %*% Gss

H[-idX,-idX] = H11
H[-idX,idX] = H12
H[idX,-idX] = t(H12)
H[idX,idX] = Gss

##Single step GBLUP
# The y vector now has dimension N = n + 2*nf, trn ids need to be updated
tst1 = tst + 2*nf
# we define a new vector of observations with NAs for founders
y1 = c(rep(NA,2*nf), y)

# phenotypes with tst individuals missing values
yNA = y1
yNA[tst1] = NA

# assumed heritability
h2 = 0.4

# mixed model equations, evaluated only in training set
MME = mme(y=yNA, V=H, h2=h2, invert=TRUE)

# or far more efficient
MME = mme(y=yNA, V=H_1, h2=h2, invert=FALSE)

# the predicted breeding values are the last n solutions
uhat = solve(MME$LHS, MME$RHS)[-1]

# corr between observed and predicted phenotypes
yhatH = uhat[tst1]

print(c('corr obs net', round(cor(yhatH,y1[tst1]),3)))

# plot
plot(y1[tst1],yhatH, main='Obs vs GBLUPss', xlab='Obs', ylab = 'GBLUPss')