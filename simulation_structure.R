# Required Libraries
library(gasper)
library(igraph)
library(ggplot2)
library(optimx)
library(Matrix)
library(imager)
library(dplyr)
library(geosphere)
library(GNAR)
library(gridExtra)
library(arrow)
library(R.matlab)
library(rhdf5)
library(ncdf4)
library(scatterplot3d)
library(fields)
library(cccd)
library(locpolExpectile)
library(grid)
library(gridBase)
library(rgl)

source("utils.R")

###########################################################
### 01. 1D: line + line heterogeneous variance structure 
###########################################################
set.seed(11)
y.hetero.chisq.01 <- c(c(1:40)/2 + 3*(rchisq(40, df=5)-5), 30-c(41:80)/4 + (rchisq(40, df=5)-5))
plot(y.hetero.chisq.01)
y.hetero.chisq.01.true <- c(c(1:40)/2, 30-c(41:80)/4)
n.hetero.chisq.01 <- length(y.hetero.chisq.01)
W.hetero.chisq.01 <- matrix(0, nrow=n.hetero.chisq.01, ncol=n.hetero.chisq.01)
for(i in 1:(n.hetero.chisq.01-1)){
  W.hetero.chisq.01[i,(i+1)] <- 1
  W.hetero.chisq.01[(i+1),i] <- 1
}

L.hetero.chisq.01 <- gasper::laplacian_mat(W.hetero.chisq.01)

## mean-based
H_list_mean.hetero.chisq.01 <- list()
GCV_list.hetero.chisq.01 <- c()
imputed_CV_list.mean.hetero.chisq.01 <- c()
lambda_candidate.mean.hetero.chisq.01 <- exp(seq(-5, 4, by=0.3)) 

for(lambda in lambda_candidate.mean.hetero.chisq.01){
  H <- solve(diag(1,length(y.hetero.chisq.01))+n.hetero.chisq.01*lambda*L.hetero.chisq.01%*%L.hetero.chisq.01)
  H_list_mean.hetero.chisq.01[[length(H_list_mean.hetero.chisq.01)+1]] <- H
  GCV_list.hetero.chisq.01 <- c(GCV_list.hetero.chisq.01, GCV(lambda, y.hetero.chisq.01, H))
  imputed_CV <- imputedCV(lambda=lambda, y=y.hetero.chisq.01, tau=0.5, alpha=0.1, B=diag(n.hetero.chisq.01), 
                          P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), n.hetero.chisq.01), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.hetero.chisq.01, k=1, M=NULL, method="mean")
  imputed_CV_list.mean.hetero.chisq.01 <- c(imputed_CV_list.mean.hetero.chisq.01, imputed_CV)
}

plot(lambda_candidate.mean.hetero.chisq.01, GCV_list.hetero.chisq.01, type="l")
plot(lambda_candidate.mean.hetero.chisq.01, imputed_CV_list.mean.hetero.chisq.01, type="l")
opt.lambda.ind.mean.hetero.chisq.01 <- which.min(GCV_list.hetero.chisq.01)
opt.lambda.mean.hetero.chisq.01 <- lambda_candidate.mean.hetero.chisq.01[opt.lambda.ind.mean.hetero.chisq.01] 
opt.lambda.ind.imputedCV.mean.hetero.chisq.01 <- which.min(imputed_CV_list.mean.hetero.chisq.01)
opt.lambda.imputedCV.mean.hetero.chisq.01 <- lambda_candidate.mean.hetero.chisq.01[opt.lambda.ind.imputedCV.mean.hetero.chisq.01] 

## quantile-based

# tau=0.9
SIC_list90.hetero.chisq.01 <- c()
GACV_list90.hetero.chisq.01 <- c()
imputed_CV_list90.hetero.chisq.01 <- c()
imputed_CV_list.mean.hetero.chisq.01 <- c()
H_list_qt90.hetero.chisq.01 <- list()

lambda_candidate90.hetero.chisq.01 <- exp(seq(-5, 2, by=0.3)) # lambda 더 낮추면 SIC, GACV 계속 낮아짐

for(lambda in lambda_candidate90.hetero.chisq.01){
  cat("####### lambda = ", lambda, "####### \n")
  res90.hetero.chisq.01 <- PIRLS(lambda=lambda, y=y.hetero.chisq.01, tau=0.9, alpha=0.1, B=diag(n.hetero.chisq.01),
                                 P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), length(y.hetero.chisq.01)), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list90.hetero.chisq.01 <- c(SIC_list90.hetero.chisq.01, res90.hetero.chisq.01$SIC)
  GACV_list90.hetero.chisq.01 <- c(GACV_list90.hetero.chisq.01, res90.hetero.chisq.01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=y.hetero.chisq.01, tau=0.9, alpha=0.1, B=diag(n.hetero.chisq.01),
                          P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), n.hetero.chisq.01), max_iterations=10000, tol=1e-6,
                          check="Oh", L=L.hetero.chisq.01, k=1, M=NULL, method="quantile")
  imputed_CV_list90.hetero.chisq.01 <- c(imputed_CV_list90.hetero.chisq.01, imputed_CV)
  H_list_qt90.hetero.chisq.01[[length(H_list_qt90.hetero.chisq.01)+1]] <- res90.hetero.chisq.01$H
}

plot(lambda_candidate90.hetero.chisq.01, SIC_list90.hetero.chisq.01, type="l")
plot(lambda_candidate90.hetero.chisq.01, GACV_list90.hetero.chisq.01, type="l")
plot(lambda_candidate90.hetero.chisq.01, imputed_CV_list90.hetero.chisq.01, type="l")
opt.lambda.ind.SIC90.hetero.chisq.01 <- which.min(SIC_list90.hetero.chisq.01)
opt.lambda.SIC90.hetero.chisq.01 <- lambda_candidate90.hetero.chisq.01[opt.lambda.ind.SIC90.hetero.chisq.01] 
opt.lambda.ind.GACV90.hetero.chisq.01 <- which.min(GACV_list90.hetero.chisq.01)
opt.lambda.GACV90.hetero.chisq.01 <- lambda_candidate90.hetero.chisq.01[opt.lambda.ind.GACV90.hetero.chisq.01] 
opt.lambda.ind.imputedCV90.hetero.chisq.01 <- which.min(imputed_CV_list90.hetero.chisq.01)
opt.lambda.imputedCV90.hetero.chisq.01 <- lambda_candidate90.hetero.chisq.01[opt.lambda.ind.imputedCV90.hetero.chisq.01]

# tau=0.75
SIC_list75.hetero.chisq.01 <- c()
GACV_list75.hetero.chisq.01 <- c()
imputed_CV_list75.hetero.chisq.01 <- c()
imputed_CV_list.mean.hetero.chisq.01 <- c()
H_list_qt75.hetero.chisq.01 <- list()

lambda_candidate75.hetero.chisq.01 <- exp(seq(-5, 2, by=0.3)) # lambda 더 낮추면 SIC, GACV 계속 낮아짐

for(lambda in lambda_candidate75.hetero.chisq.01){
  cat("####### lambda = ", lambda, "####### \n")
  res75.hetero.chisq.01 <- PIRLS(lambda=lambda, y=y.hetero.chisq.01, tau=0.75, alpha=0.1, B=diag(n.hetero.chisq.01),
                                 P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), length(y.hetero.chisq.01)), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list75.hetero.chisq.01 <- c(SIC_list75.hetero.chisq.01, res75.hetero.chisq.01$SIC)
  GACV_list75.hetero.chisq.01 <- c(GACV_list75.hetero.chisq.01, res75.hetero.chisq.01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=y.hetero.chisq.01, tau=0.75, alpha=0.1, B=diag(n.hetero.chisq.01),
                          P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), n.hetero.chisq.01), max_iterations=10000, tol=1e-6,
                          check="Oh", L=L.hetero.chisq.01, k=1, M=NULL, method="quantile")
  imputed_CV_list75.hetero.chisq.01 <- c(imputed_CV_list75.hetero.chisq.01, imputed_CV)
  H_list_qt75.hetero.chisq.01[[length(H_list_qt75.hetero.chisq.01)+1]] <- res75.hetero.chisq.01$H
}

plot(lambda_candidate75.hetero.chisq.01, SIC_list75.hetero.chisq.01, type="l")
plot(lambda_candidate75.hetero.chisq.01, GACV_list75.hetero.chisq.01, type="l")
plot(lambda_candidate75.hetero.chisq.01, imputed_CV_list75.hetero.chisq.01, type="l")
opt.lambda.ind.SIC75.hetero.chisq.01 <- which.min(SIC_list75.hetero.chisq.01)
opt.lambda.SIC75.hetero.chisq.01 <- lambda_candidate75.hetero.chisq.01[opt.lambda.ind.SIC75.hetero.chisq.01] 
opt.lambda.ind.GACV75.hetero.chisq.01 <- which.min(GACV_list75.hetero.chisq.01)
opt.lambda.GACV75.hetero.chisq.01 <- lambda_candidate75.hetero.chisq.01[opt.lambda.ind.GACV75.hetero.chisq.01] 
opt.lambda.ind.imputedCV75.hetero.chisq.01 <- which.min(imputed_CV_list75.hetero.chisq.01)
opt.lambda.imputedCV75.hetero.chisq.01 <- lambda_candidate75.hetero.chisq.01[opt.lambda.ind.imputedCV75.hetero.chisq.01]

# tau=0.5
SIC_list50.hetero.chisq.01 <- c()
GACV_list50.hetero.chisq.01 <- c()
imputed_CV_list50.hetero.chisq.01 <- c()
imputed_CV_list.mean.hetero.chisq.01 <- c()
H_list_qt50.hetero.chisq.01 <- list()

lambda_candidate50.hetero.chisq.01 <- exp(seq(-5, 2, by=0.3)) # lambda 더 낮추면 SIC, GACV 계속 낮아짐

for(lambda in lambda_candidate50.hetero.chisq.01){
  cat("####### lambda = ", lambda, "####### \n")
  res50.hetero.chisq.01 <- PIRLS(lambda=lambda, y=y.hetero.chisq.01, tau=0.5, alpha=0.1, B=diag(n.hetero.chisq.01),
                                 P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), length(y.hetero.chisq.01)), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list50.hetero.chisq.01 <- c(SIC_list50.hetero.chisq.01, res50.hetero.chisq.01$SIC)
  GACV_list50.hetero.chisq.01 <- c(GACV_list50.hetero.chisq.01, res50.hetero.chisq.01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=y.hetero.chisq.01, tau=0.5, alpha=0.1, B=diag(n.hetero.chisq.01),
                          P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), n.hetero.chisq.01), max_iterations=10000, tol=1e-6,
                          check="Oh", L=L.hetero.chisq.01, k=1, M=NULL, method="quantile")
  imputed_CV_list50.hetero.chisq.01 <- c(imputed_CV_list50.hetero.chisq.01, imputed_CV)
  H_list_qt50.hetero.chisq.01[[length(H_list_qt50.hetero.chisq.01)+1]] <- res50.hetero.chisq.01$H
}

plot(lambda_candidate50.hetero.chisq.01, SIC_list50.hetero.chisq.01, type="l")
plot(lambda_candidate50.hetero.chisq.01, GACV_list50.hetero.chisq.01, type="l")
plot(lambda_candidate50.hetero.chisq.01, imputed_CV_list50.hetero.chisq.01, type="l")
opt.lambda.ind.SIC50.hetero.chisq.01 <- which.min(SIC_list50.hetero.chisq.01)
opt.lambda.SIC50.hetero.chisq.01 <- lambda_candidate50.hetero.chisq.01[opt.lambda.ind.SIC50.hetero.chisq.01] 
opt.lambda.ind.GACV50.hetero.chisq.01 <- which.min(GACV_list50.hetero.chisq.01)
opt.lambda.GACV50.hetero.chisq.01 <- lambda_candidate50.hetero.chisq.01[opt.lambda.ind.GACV50.hetero.chisq.01] 
opt.lambda.ind.imputedCV50.hetero.chisq.01 <- which.min(imputed_CV_list50.hetero.chisq.01)
opt.lambda.imputedCV50.hetero.chisq.01 <- lambda_candidate50.hetero.chisq.01[opt.lambda.ind.imputedCV50.hetero.chisq.01]

# tau=0.25
SIC_list25.hetero.chisq.01 <- c()
GACV_list25.hetero.chisq.01 <- c()
imputed_CV_list25.hetero.chisq.01 <- c()
imputed_CV_list.mean.hetero.chisq.01 <- c()
H_list_qt25.hetero.chisq.01 <- list()

lambda_candidate25.hetero.chisq.01 <- exp(seq(-5, 2, by=0.3)) # lambda 더 낮추면 SIC, GACV 계속 낮아짐

for(lambda in lambda_candidate25.hetero.chisq.01){
  cat("####### lambda = ", lambda, "####### \n")
  res25.hetero.chisq.01 <- PIRLS(lambda=lambda, y=y.hetero.chisq.01, tau=0.25, alpha=0.1, B=diag(n.hetero.chisq.01),
                                 P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), length(y.hetero.chisq.01)), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list25.hetero.chisq.01 <- c(SIC_list25.hetero.chisq.01, res25.hetero.chisq.01$SIC)
  GACV_list25.hetero.chisq.01 <- c(GACV_list25.hetero.chisq.01, res25.hetero.chisq.01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=y.hetero.chisq.01, tau=0.25, alpha=0.1, B=diag(n.hetero.chisq.01),
                          P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), n.hetero.chisq.01), max_iterations=10000, tol=1e-6,
                          check="Oh", L=L.hetero.chisq.01, k=1, M=NULL, method="quantile")
  imputed_CV_list25.hetero.chisq.01 <- c(imputed_CV_list25.hetero.chisq.01, imputed_CV)
  H_list_qt25.hetero.chisq.01[[length(H_list_qt25.hetero.chisq.01)+1]] <- res25.hetero.chisq.01$H
}

plot(lambda_candidate25.hetero.chisq.01, SIC_list25.hetero.chisq.01, type="l")
plot(lambda_candidate25.hetero.chisq.01, GACV_list25.hetero.chisq.01, type="l")
plot(lambda_candidate25.hetero.chisq.01, imputed_CV_list25.hetero.chisq.01, type="l")
opt.lambda.ind.SIC25.hetero.chisq.01 <- which.min(SIC_list25.hetero.chisq.01)
opt.lambda.SIC25.hetero.chisq.01 <- lambda_candidate25.hetero.chisq.01[opt.lambda.ind.SIC25.hetero.chisq.01] 
opt.lambda.ind.GACV25.hetero.chisq.01 <- which.min(GACV_list25.hetero.chisq.01)
opt.lambda.GACV25.hetero.chisq.01 <- lambda_candidate25.hetero.chisq.01[opt.lambda.ind.GACV25.hetero.chisq.01] 
opt.lambda.ind.imputedCV25.hetero.chisq.01 <- which.min(imputed_CV_list25.hetero.chisq.01)
opt.lambda.imputedCV25.hetero.chisq.01 <- lambda_candidate25.hetero.chisq.01[opt.lambda.ind.imputedCV25.hetero.chisq.01]

# tau=0.1
SIC_list10.hetero.chisq.01 <- c()
GACV_list10.hetero.chisq.01 <- c()
imputed_CV_list10.hetero.chisq.01 <- c()
imputed_CV_list.mean.hetero.chisq.01 <- c()
H_list_qt10.hetero.chisq.01 <- list()

lambda_candidate10.hetero.chisq.01 <- exp(seq(-5, 2, by=0.3)) # lambda 더 낮추면 SIC, GACV 계속 낮아짐

for(lambda in lambda_candidate10.hetero.chisq.01){
  cat("####### lambda = ", lambda, "####### \n")
  res10.hetero.chisq.01 <- PIRLS(lambda=lambda, y=y.hetero.chisq.01, tau=0.1, alpha=0.1, B=diag(n.hetero.chisq.01),
                                 P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), length(y.hetero.chisq.01)), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list10.hetero.chisq.01 <- c(SIC_list10.hetero.chisq.01, res10.hetero.chisq.01$SIC)
  GACV_list10.hetero.chisq.01 <- c(GACV_list10.hetero.chisq.01, res10.hetero.chisq.01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=y.hetero.chisq.01, tau=0.1, alpha=0.1, B=diag(n.hetero.chisq.01),
                          P=L.hetero.chisq.01%*%L.hetero.chisq.01, gamma_init=rep(mean(y.hetero.chisq.01), n.hetero.chisq.01), max_iterations=10000, tol=1e-6,
                          check="Oh", L=L.hetero.chisq.01, k=1, M=NULL, method="quantile")
  imputed_CV_list10.hetero.chisq.01 <- c(imputed_CV_list10.hetero.chisq.01, imputed_CV)
  H_list_qt10.hetero.chisq.01[[length(H_list_qt10.hetero.chisq.01)+1]] <- res10.hetero.chisq.01$H
}

plot(lambda_candidate10.hetero.chisq.01, SIC_list10.hetero.chisq.01, type="l")
plot(lambda_candidate10.hetero.chisq.01, GACV_list10.hetero.chisq.01, type="l")
plot(lambda_candidate10.hetero.chisq.01, imputed_CV_list10.hetero.chisq.01, type="l")
opt.lambda.ind.SIC10.hetero.chisq.01 <- which.min(SIC_list10.hetero.chisq.01)
opt.lambda.SIC10.hetero.chisq.01 <- lambda_candidate10.hetero.chisq.01[opt.lambda.ind.SIC10.hetero.chisq.01] 
opt.lambda.ind.GACV10.hetero.chisq.01 <- which.min(GACV_list10.hetero.chisq.01)
opt.lambda.GACV10.hetero.chisq.01 <- lambda_candidate10.hetero.chisq.01[opt.lambda.ind.GACV10.hetero.chisq.01] 
opt.lambda.ind.imputedCV10.hetero.chisq.01 <- which.min(imputed_CV_list10.hetero.chisq.01)
opt.lambda.imputedCV10.hetero.chisq.01 <- lambda_candidate10.hetero.chisq.01[opt.lambda.ind.imputedCV10.hetero.chisq.01]


par(mfrow=c(1,1))
par(mar=c(5,5.5,4,2)+0.1, oma=c(0,0,0,0))
plot(y.hetero.chisq.01, ylab="value", xlab="node index",cex.lab=3.2, cex.axis=2.5, cex=2, lwd=2)
lines(y.hetero.chisq.01.true, lwd=2)
lines(as.numeric(H_list_mean.hetero.chisq.01[[opt.lambda.ind.imputedCV.mean.hetero.chisq.01]]%*%y.hetero.chisq.01), col="blue", lwd=2) 

lines(as.numeric(H_list_qt50.hetero.chisq.01[[opt.lambda.ind.imputedCV50.hetero.chisq.01]]%*%y.hetero.chisq.01), col="red", lwd=2)
lines(as.numeric(H_list_qt90.hetero.chisq.01[[opt.lambda.ind.imputedCV90.hetero.chisq.01]]%*%y.hetero.chisq.01), col="red", lty=2, lwd=2)
lines(as.numeric(H_list_qt75.hetero.chisq.01[[opt.lambda.ind.imputedCV75.hetero.chisq.01]]%*%y.hetero.chisq.01), col="red", lty=3, lwd=2)
lines(as.numeric(H_list_qt25.hetero.chisq.01[[opt.lambda.ind.imputedCV25.hetero.chisq.01]]%*%y.hetero.chisq.01), col="red", lty=4, lwd=2)
lines(as.numeric(H_list_qt10.hetero.chisq.01[[opt.lambda.ind.imputedCV10.hetero.chisq.01]]%*%y.hetero.chisq.01), col="red", lty=5, lwd=2)

legend(x = "bottomright", lty = c(1,1,2,3,1,4,5), text.font = 4, cex=1.8, 
       col= c("black","blue","red","red","red","red","red"), text.col = "black", 
       legend=c("underlying","mean","tau = 0.9","tau = 0.75","tau = 0.5","tau = 0.25","tau = 0.1"), lwd=2)


#######################################################
### 02. 2D: heterogeneous variance structure (lattice)
#######################################################
g1 = make_lattice(c(15,15)) #make 30 x 30 grid graph
x = rep(0:14, 15) #define coordinates of vertices
y = rep(0:14, rep(15,15))
graph_attr(g1, 'xy') = data.frame(x,y) #define 'xy' attributes of graph g which represents coordinates
g1.ad_mat = get.adjacency(g1, sparse=T) #get unweighted adjacency matrix of graph g
graph_attr(g1, 'sA') = g1.ad_mat

lattice01 <- list()
lattice01$xy <- g1$xy
n.lattice01 <- nrow(lattice01$xy)
rownames(lattice01$xy) <- 1:n.lattice01
N <- 100

# weight matrix
lattice01$A <- g1$sA
edge.wt <- igraph::as_data_frame(g1, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              g1$sA[edge.wt[i,1], edge.wt[i,2]]))
}

# sparse weight matrix
lattice01$sA <- sp.wmat
L.lattice01 <- gasper::laplacian_mat(lattice01$A)

# signal assignment
lattice01$f.true <- (lattice01$xy$x+lattice01$xy$y)/3 # true sine signal

q1 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = lattice01$f.true, 
                         min=0, max=10, value="value") + ggtitle("true signal") # plot graph signal


set.seed(100)
lowertri <- rownames(lattice01$xy[lattice01$xy$x + lattice01$xy$y <= 13,])
lattice01$f <- lattice01$f.true
lattice01$f[as.numeric(lowertri)] <- lattice01$f[as.numeric(lowertri)] + 0.5*(rchisq(length(lowertri), df=4)-4)
lattice01$f[-as.numeric(lowertri)] <- lattice01$f[-as.numeric(lowertri)] + 0.1*(rchisq(225-length(lowertri), df=4)-4)


q2 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = lattice01$f, 
                         min=-1.5, max=9.5, value="value") + ggtitle("original signal") # plot graph signal


## mean-based
H_list_mean.lattice01 <- list()
GCV_list.lattice01 <- c()
imputed_CV_list.mean.lattice01 <- c()
lambda_candidate.mean.lattice01 <- exp(seq(-10, -2, by=0.3))
for(lambda in lambda_candidate.mean.lattice01){
  H <- solve(diag(1,n.lattice01)+n.lattice01*lambda*L.lattice01%*%L.lattice01)
  H_list_mean.lattice01[[length(H_list_mean.lattice01)+1]] <- H
  GCV_list.lattice01 <- c(GCV_list.lattice01, GCV(lambda, lattice01$f, H))
  imputed_CV <- imputedCV(lambda=lambda, y=lattice01$f, tau=0.5, alpha=0.1, B=diag(n.lattice01), 
                          P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice01, k=1, M=NULL, method="mean")
  imputed_CV_list.mean.lattice01 <- c(imputed_CV_list.mean.lattice01, imputed_CV)
}
plot(lambda_candidate.mean.lattice01, GCV_list.lattice01, type="l")
plot(lambda_candidate.mean.lattice01, imputed_CV_list.mean.lattice01, type="l")
opt.lambda.ind.mean.lattice01 <- which.min(GCV_list.lattice01)
opt.lambda.mean.lattice01 <- lambda_candidate.mean.lattice01[opt.lambda.ind.mean.lattice01] 
opt.lambda.ind.imputedCV.mean.lattice01 <- which.min(imputed_CV_list.mean.lattice01)
opt.lambda.imputedCV.mean.lattice01 <- lambda_candidate.mean.lattice01[opt.lambda.ind.imputedCV.mean.lattice01] 


# tau=0.9
SIC_list90.lattice01 <- c()
GACV_list90.lattice01 <- c()
imputed_CV_list90.lattice01 <- c()
H_list_qt90.lattice01 <- list()

lambda_candidate90.lattice01 <- exp(seq(-7, -3, by=0.4))
for(lambda in lambda_candidate90.lattice01){
  cat("####### lambda = ", lambda, "####### \n")
  res90.lattice01 <- PIRLS(lambda=lambda, y=lattice01$f, tau=0.9, alpha=0.1, B=diag(n.lattice01), 
                           P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list90.lattice01 <- c(SIC_list90.lattice01, res90.lattice01$SIC)
  GACV_list90.lattice01 <- c(GACV_list90.lattice01, res90.lattice01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice01$f, tau=0.9, alpha=0.1, B=diag(n.lattice01), 
                          P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice01, k=1, M=NULL, method="quantile")
  imputed_CV_list90.lattice01 <- c(imputed_CV_list90.lattice01, imputed_CV)
  H_list_qt90.lattice01[[length(H_list_qt90.lattice01)+1]] <- res90.lattice01$H
}
plot(lambda_candidate90.lattice01, SIC_list90.lattice01, type="l")
plot(lambda_candidate90.lattice01, GACV_list90.lattice01, type="l")
plot(lambda_candidate90.lattice01, imputed_CV_list90.lattice01, type="l")
opt.lambda.ind.SIC90.lattice01 <- which.min(SIC_list90.lattice01)
opt.lambda.SIC90.lattice01 <- lambda_candidate90.lattice01[opt.lambda.ind.SIC90.lattice01] 
opt.lambda.ind.GACV90.lattice01 <- which.min(GACV_list90.lattice01)
opt.lambda.GACV90.lattice01 <- lambda_candidate90.lattice01[opt.lambda.ind.GACV90.lattice01] 
opt.lambda.ind.imputedCV90.lattice01 <- which.min(imputed_CV_list90.lattice01)
opt.lambda.imputedCV90.lattice01 <- lambda_candidate90.lattice01[opt.lambda.ind.imputedCV90.lattice01] 


# tau=0.75
SIC_list75.lattice01 <- c()
GACV_list75.lattice01 <- c()
imputed_CV_list75.lattice01 <- c()
H_list_qt75.lattice01 <- list()

lambda_candidate75.lattice01 <- exp(seq(-7, -3, by=0.4))
for(lambda in lambda_candidate75.lattice01){
  cat("####### lambda = ", lambda, "####### \n")
  res75.lattice01 <- PIRLS(lambda=lambda, y=lattice01$f, tau=0.75, alpha=0.1, B=diag(n.lattice01), 
                           P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list75.lattice01 <- c(SIC_list75.lattice01, res75.lattice01$SIC)
  GACV_list75.lattice01 <- c(GACV_list75.lattice01, res75.lattice01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice01$f, tau=0.75, alpha=0.1, B=diag(n.lattice01), 
                          P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice01, k=1, M=NULL, method="quantile")
  imputed_CV_list75.lattice01 <- c(imputed_CV_list75.lattice01, imputed_CV)
  H_list_qt75.lattice01[[length(H_list_qt75.lattice01)+1]] <- res75.lattice01$H
}
plot(lambda_candidate75.lattice01, SIC_list75.lattice01, type="l")
plot(lambda_candidate75.lattice01, GACV_list75.lattice01, type="l")
plot(lambda_candidate75.lattice01, imputed_CV_list75.lattice01, type="l")
opt.lambda.ind.SIC75.lattice01 <- which.min(SIC_list75.lattice01)
opt.lambda.SIC75.lattice01 <- lambda_candidate75.lattice01[opt.lambda.ind.SIC75.lattice01] 
opt.lambda.ind.GACV75.lattice01 <- which.min(GACV_list75.lattice01)
opt.lambda.GACV75.lattice01 <- lambda_candidate75.lattice01[opt.lambda.ind.GACV75.lattice01] 
opt.lambda.ind.imputedCV75.lattice01 <- which.min(imputed_CV_list75.lattice01)
opt.lambda.imputedCV75.lattice01 <- lambda_candidate75.lattice01[opt.lambda.ind.imputedCV75.lattice01] 


# tau=0.5
SIC_list50.lattice01 <- c()
GACV_list50.lattice01 <- c()
imputed_CV_list50.lattice01 <- c()
H_list_qt50.lattice01 <- list()

lambda_candidate50.lattice01 <- exp(seq(-7, -3, by=0.4))
for(lambda in lambda_candidate50.lattice01){
  cat("####### lambda = ", lambda, "####### \n")
  res50.lattice01 <- PIRLS(lambda=lambda, y=lattice01$f, tau=0.5, alpha=0.1, B=diag(n.lattice01), 
                           P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list50.lattice01 <- c(SIC_list50.lattice01, res50.lattice01$SIC)
  GACV_list50.lattice01 <- c(GACV_list50.lattice01, res50.lattice01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice01$f, tau=0.5, alpha=0.1, B=diag(n.lattice01), 
                          P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice01, k=1, M=NULL, method="quantile")
  imputed_CV_list50.lattice01 <- c(imputed_CV_list50.lattice01, imputed_CV)
  H_list_qt50.lattice01[[length(H_list_qt50.lattice01)+1]] <- res50.lattice01$H
}
plot(lambda_candidate50.lattice01, SIC_list50.lattice01, type="l")
plot(lambda_candidate50.lattice01, GACV_list50.lattice01, type="l")
plot(lambda_candidate50.lattice01, imputed_CV_list50.lattice01, type="l")
opt.lambda.ind.SIC50.lattice01 <- which.min(SIC_list50.lattice01)
opt.lambda.SIC50.lattice01 <- lambda_candidate50.lattice01[opt.lambda.ind.SIC50.lattice01] 
opt.lambda.ind.GACV50.lattice01 <- which.min(GACV_list50.lattice01)
opt.lambda.GACV50.lattice01 <- lambda_candidate50.lattice01[opt.lambda.ind.GACV50.lattice01] 
opt.lambda.ind.imputedCV50.lattice01 <- which.min(imputed_CV_list50.lattice01)
opt.lambda.imputedCV50.lattice01 <- lambda_candidate50.lattice01[opt.lambda.ind.imputedCV50.lattice01] 


# tau=0.25
SIC_list25.lattice01 <- c()
GACV_list25.lattice01 <- c()
imputed_CV_list25.lattice01 <- c()
H_list_qt25.lattice01 <- list()

lambda_candidate25.lattice01 <- exp(seq(-7, -3, by=0.4))
for(lambda in lambda_candidate25.lattice01){
  cat("####### lambda = ", lambda, "####### \n")
  res25.lattice01 <- PIRLS(lambda=lambda, y=lattice01$f, tau=0.25, alpha=0.1, B=diag(n.lattice01), 
                           P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list25.lattice01 <- c(SIC_list25.lattice01, res25.lattice01$SIC)
  GACV_list25.lattice01 <- c(GACV_list25.lattice01, res25.lattice01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice01$f, tau=0.25, alpha=0.1, B=diag(n.lattice01), 
                          P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice01, k=1, M=NULL, method="quantile")
  imputed_CV_list25.lattice01 <- c(imputed_CV_list25.lattice01, imputed_CV)
  H_list_qt25.lattice01[[length(H_list_qt25.lattice01)+1]] <- res25.lattice01$H
}
plot(lambda_candidate25.lattice01, SIC_list25.lattice01, type="l")
plot(lambda_candidate25.lattice01, GACV_list25.lattice01, type="l")
plot(lambda_candidate25.lattice01, imputed_CV_list25.lattice01, type="l")
opt.lambda.ind.SIC25.lattice01 <- which.min(SIC_list25.lattice01)
opt.lambda.SIC25.lattice01 <- lambda_candidate25.lattice01[opt.lambda.ind.SIC25.lattice01] 
opt.lambda.ind.GACV25.lattice01 <- which.min(GACV_list25.lattice01)
opt.lambda.GACV25.lattice01 <- lambda_candidate25.lattice01[opt.lambda.ind.GACV25.lattice01] 
opt.lambda.ind.imputedCV25.lattice01 <- which.min(imputed_CV_list25.lattice01)
opt.lambda.imputedCV25.lattice01 <- lambda_candidate25.lattice01[opt.lambda.ind.imputedCV25.lattice01] 


# tau=0.1
SIC_list10.lattice01 <- c()
GACV_list10.lattice01 <- c()
imputed_CV_list10.lattice01 <- c()
H_list_qt10.lattice01 <- list()

lambda_candidate10.lattice01 <- exp(seq(-7, -3, by=0.4))
for(lambda in lambda_candidate10.lattice01){
  cat("####### lambda = ", lambda, "####### \n")
  res10.lattice01 <- PIRLS(lambda=lambda, y=lattice01$f, tau=0.1, alpha=0.1, B=diag(n.lattice01), 
                           P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list10.lattice01 <- c(SIC_list10.lattice01, res10.lattice01$SIC)
  GACV_list10.lattice01 <- c(GACV_list10.lattice01, res10.lattice01$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice01$f, tau=0.1, alpha=0.1, B=diag(n.lattice01), 
                          P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice01, k=1, M=NULL, method="quantile")
  imputed_CV_list10.lattice01 <- c(imputed_CV_list10.lattice01, imputed_CV)
  H_list_qt10.lattice01[[length(H_list_qt10.lattice01)+1]] <- res10.lattice01$H
}
plot(lambda_candidate10.lattice01, SIC_list10.lattice01, type="l")
plot(lambda_candidate10.lattice01, GACV_list10.lattice01, type="l")
plot(lambda_candidate10.lattice01, imputed_CV_list10.lattice01, type="l")
opt.lambda.ind.SIC10.lattice01 <- which.min(SIC_list10.lattice01)
opt.lambda.SIC10.lattice01 <- lambda_candidate10.lattice01[opt.lambda.ind.SIC10.lattice01] 
opt.lambda.ind.GACV10.lattice01 <- which.min(GACV_list10.lattice01)
opt.lambda.GACV10.lattice01 <- lambda_candidate10.lattice01[opt.lambda.ind.GACV10.lattice01] 
opt.lambda.ind.imputedCV10.lattice01 <- which.min(imputed_CV_list10.lattice01)
opt.lambda.imputedCV10.lattice01 <- lambda_candidate10.lattice01[opt.lambda.ind.imputedCV10.lattice01] 


# fitting results
q1 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = lattice01$f.true, 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("underlying signal") # plot graph signal
q2 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = lattice01$f, 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("noisy signal") # plot graph signal
q3 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f), 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("mean") # plot graph signal
q4 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt90.lattice01[[opt.lambda.ind.imputedCV90.lattice01]]%*%lattice01$f), 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.9") # plot graph signal
q5 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt75.lattice01[[opt.lambda.ind.imputedCV75.lattice01]]%*%lattice01$f), 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.75") # plot graph signal
q6 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f), 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.5") # plot graph signal
q7 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt25.lattice01[[opt.lambda.ind.imputedCV25.lattice01]]%*%lattice01$f), 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.25") # plot graph signal
q8 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt10.lattice01[[opt.lambda.ind.imputedCV10.lattice01]]%*%lattice01$f), 
                         min=-1.5, max=9.5, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.1") # plot graph signal

grid.arrange(q1,q3,q4,q7,q2,q6,q5,q8, ncol=4)

# # scaled persp plots
# par(mfrow=c(2,3))
# persp(1:15, 1:15, matrix(lattice01$f.true, nrow=15), theta = -40, phi = 20, expand = 0.5, col = "gray",
#       xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
#       cex.lab=3)
# title(main = "underlying signal", cex.main = 4, line=1)
# persp(1:15, 1:15, matrix(lattice01$f, nrow=15), theta = -40, phi = 20, expand = 0.5, col = "gray",
#       xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
#       cex.lab=3)
# title(main = "noisy signal", cex.main = 4, line=1)
# persp(1:15, 1:15, matrix(as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f), nrow=15),
#       theta = -40, phi = 20, expand = 0.5, col = "lightblue",
#       xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
#       cex.lab=3)
# title(main = "mean (GICV)", cex.main = 4, line=1)
# persp(1:15, 1:15, matrix(as.numeric(H_list_qt90.lattice01[[opt.lambda.ind.imputedCV90.lattice01]]%*%lattice01$f), nrow=15),
#       theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
#       xlab="column", ylab="row", zlab="value", zlim=range(t(matrix(lattice01$f, nrow=15))[15:1,]),
#       cex.lab=3)
# title(main = "tau = 0.9 (GIQCV)", cex.main = 4, line=1)
# persp(1:15, 1:15, matrix(as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f), nrow=15),
#       theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
#       xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
#       cex.lab=3)
# title(main = "tau = 0.5 (GIQCV)", cex.main = 4, line=1)
# persp(1:15, 1:15, matrix(as.numeric(H_list_qt10.lattice01[[opt.lambda.ind.imputedCV10.lattice01]]%*%lattice01$f), nrow=15),
#       theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
#       xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
#       cex.lab=3)
# title(main = "tau = 0.1 (GIQCV)", cex.main = 4, line=1)

# scaled persp plots
par(mfrow=c(2,3), mar=c(5,4,6,2)+0.1)
persp(1:15, 1:15, matrix(lattice01$f.true, nrow=15), theta = -40, phi = 20, expand = 0.5, col = "gray",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "underlying signal", cex.main = 4, line=2.5)
persp(1:15, 1:15, matrix(lattice01$f, nrow=15), theta = -40, phi = 20, expand = 0.5, col = "gray",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "noisy signal", cex.main = 4, line=2.5)
persp(1:15, 1:15, matrix(as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightblue",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "mean", cex.main = 4, line=2.5)
persp(1:15, 1:15, matrix(as.numeric(H_list_qt90.lattice01[[opt.lambda.ind.imputedCV90.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
      xlab="column", ylab="row", zlab="value", zlim=range(t(matrix(lattice01$f, nrow=15))[15:1,]),
      cex.lab=3)
title(main = "tau = 0.9", cex.main = 4, line=2.5)
persp(1:15, 1:15, matrix(as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "tau = 0.5", cex.main = 4, line=2.5)
persp(1:15, 1:15, matrix(as.numeric(H_list_qt10.lattice01[[opt.lambda.ind.imputedCV10.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "tau = 0.1", cex.main = 4, line=2.5)

