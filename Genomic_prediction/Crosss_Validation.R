library(BGLR)
library(rrBLUP)
data(wheat)
dim(wheat.Y)
dim(wheat.X)
x = wheat.X
y = wheat.Y[,1]

# i) Creating the fold
n=559                               # original number of individuals
k=10                                # how many folds/partitions
set.seed(123)                       # same fold/partitions division every time
folds=sample(1:k,size=n,replace=T)  # creating random partitions 
folds                               # ind assigned to diff folds
table(folds)  

# Fitting Genomic selection using 10-fold CV scheme
r2 = vector()                        # empty vector to save the accuracy values

for(i in 1:max(folds)){              # loop across the folds
  print(i)
  # ii. Creating the test data set (individuals that we will predict)
  test = which(folds==i)
  y_train=y[-test]
  y_test = y[test]
  
  x_train= x[-test,]
  x_test= x[test,]
  
  # iii. Fitting the model in the Training data set
  fit = mixed.solve(y_train, x_train)
  
  # iv. Prediction in the test
  yhat = x_test %*% fit$u
  
  # vi. Pearson's correlaction
  r2[i] = cor(yhat, y_test,use="complete")
}

r2
mean(r2)
boxplot(r2)
