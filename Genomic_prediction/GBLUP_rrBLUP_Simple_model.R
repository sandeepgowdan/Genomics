library(BGLR)
library(rrBLUP)

# Download the data
data(wheat)
dim(wheat.A)
dim(wheat.X)
dim(wheat.Y)

# Fit RRBLUP
fit1 <- mixed.solve(y = wheat.Y[,2], Z = wheat.X) # RR-BLUP

gebv = wheat.X %*% fit1$u ## GEBV
gebv
cor(gebv,wheat.Y[,1])
plot(abs(fit1$u))

# Fit GBLUP
G =A.mat(wheat.X)
fit2 <- mixed.solve(y = wheat.Y[,1], K = G) # GBLUP
cor(fit2$u,wheat.Y[,1])
plot(gebv, fit2$u, main = "RR-BLUP vs. GBLUP")
fit2$Vu/(fit2$Ve+fit2$Vu) # h2

# Accuracy in TST data set  ##evaluate the model perfrormance (can also be done by cross validation)
y_tst = wheat.Y[1:100,1]
x_tst = wheat.X[1:100,]

y_trn = wheat.Y[-c(1:100),1]
x_trn = wheat.X[-c(1:100),]

fit4 <- mixed.solve(y = y_trn, Z = x_trn) # RR-BLUP
gebv_tst = x_tst %*% fit4$u ## GEBV
cor(gebv_tst,y_tst)
    