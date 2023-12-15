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


source("utils.R")

###################################
### 01. 1D: line w/ outliers 
###################################
y01.true <- c(1:100)/2 # true signal
plot(y01.true)
lines(y01.true)
N <- 100
n01 <- length(y01.true)
Y01 <- t(replicate(N, y01.true))
outlier.ind.mat <- matrix(0, nrow=N, ncol=round(n01/10))
set.seed(10)
for(i in 1:N){
  outlier.ind <- sort(sample(n01, round(n01/10)))
  outlier.ind.mat[i,] <- outlier.ind
  Y01[i,outlier.ind] <- Y01[i,outlier.ind] + rnorm(round(n01/10),0,20) # generate 10% outliers
  Y01[i,-outlier.ind] <- Y01[i,-outlier.ind] + rnorm(n01-round(n01/10),0,1) 
}

W01 <- matrix(0, nrow=n01, ncol=n01)
for(i in 1:(n01-1)){
  W01[i,(i+1)] <- 1
  W01[(i+1),i] <- 1
}

L01 <- gasper::laplacian_mat(W01) # compute graph Laplacian

MSE.mean.01 <- NULL
MSE.imputedCV.mean.01 <- NULL
MSE.SIC50.01 <- NULL
MSE.GACV50.01 <- NULL
MSE.imputedCV50.01 <- NULL 
opt.lambda.SIC50.01.list <- NULL
opt.lambda.GACV50.01.list <- NULL
opt.lambda.imputedCV50.01.list <- NULL
opt.lambda.mean.01.list <- NULL
opt.lambda.imputedCV.mean.01.list <- NULL
for(i in 1:N){
  y01 <- Y01[i,]
  # 50% quantile-based
  SIC_list50.01 <- c()
  GACV_list50.01 <- c()
  imputed_CV_list50.01 <- c()
  imputed_CV_list.mean.01 <- c()
  H_list_qt50.01 <- list()
  
  lambda_candidate50.01 <- exp(seq(-5, floor(log(n01)), by=0.3))
  for(lambda in lambda_candidate50.01){
    cat("####### lambda = ", lambda, "####### \n")
    res50.01 <- PIRLS(lambda=lambda, y=y01, tau=0.5, alpha=0.1, B=diag(n01),
                      P=L01%*%L01, gamma_init=rep(mean(y01), length(y01)), max_iterations=10000, tol=1e-6, check="Oh")
    SIC_list50.01 <- c(SIC_list50.01, res50.01$SIC)
    GACV_list50.01 <- c(GACV_list50.01, res50.01$GACV)
    imputed_CV <- imputedCV(lambda=lambda, y=y01, tau=0.5, alpha=0.1, B=diag(n01), 
                              P=L01%*%L01, gamma_init=rep(mean(y01), n01), max_iterations=10000, tol=1e-6, 
                              check="Oh", L=L01, k=1, M=NULL, method="quantile")
    imputed_CV_list50.01 <- c(imputed_CV_list50.01, imputed_CV)
    H_list_qt50.01[[length(H_list_qt50.01)+1]] <- res50.01$H
  }
  
  # plot(lambda_candidate50.01, SIC_list50.01, type="l")
  # plot(lambda_candidate50.01, GACV_list50.01, type="l")
  # plot(lambda_candidate50.01, imputed_CV_list50.01, type="l")
  opt.lambda.ind.SIC50.01 <- which.min(SIC_list50.01)
  opt.lambda.SIC50.01 <- lambda_candidate50.01[opt.lambda.ind.SIC50.01] 
  opt.lambda.ind.GACV50.01 <- which.min(GACV_list50.01)
  opt.lambda.GACV50.01 <- lambda_candidate50.01[opt.lambda.ind.GACV50.01] 
  opt.lambda.ind.imputedCV50.01 <- which.min(imputed_CV_list50.01)
  opt.lambda.imputedCV50.01 <- lambda_candidate50.01[opt.lambda.ind.imputedCV50.01] 
  
  # mean-based
  H_list_mean.01 <- list()
  GCV_list.01 <- c()
  lambda_candidate.mean.01 <- exp(seq(-5, floor(log(n01)), by=0.3)) 
  
  for(lambda in lambda_candidate.mean.01){
    H <- solve(diag(1,length(y01))+n01*lambda*L01%*%L01)
    H_list_mean.01[[length(H_list_mean.01)+1]] <- H
    GCV_list.01 <- c(GCV_list.01, GCV(lambda, y01, H))
    imputed_CV <- imputedCV(lambda=lambda, y=y01, tau=0.5, alpha=0.1, B=diag(n01), 
                              P=L01%*%L01, gamma_init=rep(mean(y01), n01), max_iterations=10000, tol=1e-6, 
                              check="Oh", L=L01, k=1, M=NULL, method="mean")
    imputed_CV_list.mean.01 <- c(imputed_CV_list.mean.01, imputed_CV)
  }
  
  # plot(lambda_candidate.mean.01, GCV_list.01, type="l")
  # plot(lambda_candidate.mean.01, imputed_CV_list.mean.01, type="l")
  opt.lambda.ind.mean.01 <- which.min(GCV_list.01)
  opt.lambda.mean.01 <- lambda_candidate.mean.01[opt.lambda.ind.mean.01] 
  opt.lambda.ind.imputedCV.mean.01 <- which.min(imputed_CV_list.mean.01)
  opt.lambda.imputedCV.mean.01 <- lambda_candidate.mean.01[opt.lambda.ind.imputedCV.mean.01]
  
  # calculate MSE
  MSE.mean.01 <- c(MSE.mean.01, mse(y01.true, H_list_mean.01[[opt.lambda.ind.mean.01]]%*%y01))
  MSE.imputedCV.mean.01 <- c(MSE.imputedCV.mean.01, mse(y01.true, H_list_mean.01[[opt.lambda.ind.imputedCV.mean.01]]%*%y01))
  MSE.SIC50.01 <- c(MSE.SIC50.01, mse(y01.true, H_list_qt50.01[[opt.lambda.ind.SIC50.01]]%*%y01))
  MSE.GACV50.01 <- c(MSE.GACV50.01, mse(y01.true, H_list_qt50.01[[opt.lambda.ind.GACV50.01]]%*%y01))
  MSE.imputedCV50.01 <- c(MSE.imputedCV50.01, mse(y01.true, H_list_qt50.01[[opt.lambda.ind.imputedCV50.01]]%*%y01))
  
  # save optimal lambda
  opt.lambda.SIC50.01.list <- c(opt.lambda.SIC50.01.list, opt.lambda.SIC50.01)
  opt.lambda.GACV50.01.list <- c(opt.lambda.GACV50.01.list, opt.lambda.GACV50.01)
  opt.lambda.imputedCV50.01.list <- c(opt.lambda.imputedCV50.01.list, opt.lambda.imputedCV50.01)
  opt.lambda.mean.01.list <- c(opt.lambda.mean.01.list, opt.lambda.mean.01)
  opt.lambda.imputedCV.mean.01.list <- c(opt.lambda.imputedCV.mean.01.list, opt.lambda.imputedCV.mean.01)
}


# visualize fit
par(mfrow=c(1,1))
par(mar=c(5,5.5,4,2)+0.1, oma=c(0,0,0,0))
plot(y01, ylab="value", xlab="node index", cex.lab=3.5, cex.axis=2.5, cex=2)
lines(y01.true, lwd=2)
lines(as.numeric(H_list_mean.01[[opt.lambda.ind.imputedCV.mean.01]]%*%y01), col="blue", lwd=2) 
lines(as.numeric(H_list_qt50.01[[opt.lambda.ind.imputedCV50.01]]%*%y01), col="red", lwd=2) 

legend(x = "topleft", lty = c(1,1,1), text.font = 4, cex=2, 
       col= c("black","blue","red"), text.col = "black", 
       legend=c("true","mean","tau = 0.5"), lwd=2)

# MSE plot
par(mar=c(5,5.5,4,2)+0.1)
plot(MSE.imputedCV.mean.01_1_100, type="l", lwd=2, ylab="MSE", xlab="replication index", ylim=c(0,7), cex.lab=3.2, cex.axis=2.5, cex=2)
points(MSE.imputedCV.mean.01_1_100, pch=19)
lines(MSE.imputedCV50.01_1_100, lwd=2, col="red")
points(MSE.imputedCV50.01_1_100, lwd=2, pch=19, col="red")

legend(x = "topleft", lty = c(1,1,1), text.font = 4, cex=2, 
       col= c("black","red"), text.col = "black", 
       legend=c("mean","tau = 0.5"), lwd=2)


par(mar=c(5,5.5,4,2)+0.1)
plot(MSE.imputedCV.mean.01, type="l", ylim=c(0,7),
     ylab="MSE", xlab="replication index", col="black", lwd=2, cex.lab=3.2, cex.axis=2.5, cex=2)
lines(MSE.mean.01, col="blue", lwd=2)
lines(MSE.imputedCV50.01, col="red", lwd=2)
lines(MSE.SIC50.01, col="magenta", lwd=2)
lines(MSE.GACV50.01, col="cyan", lwd=2)

legend(x = "topright", lty = c(1,1,1,1,1), text.font = 4, cex=2, 
       col= c("black","blue","red","magenta", "cyan"), text.col = "black", 
       legend=c("GICV","GCV","GIQCV", "SIC","GACV"), lwd=2)

#######################################
### 02. 2D: sine w/ outliers (lattice)  
#######################################
g1 = make_lattice(c(15,15)) #make 30 x 30 grid graph
x = rep(0:14, 15) #define coordinates of vertices
y = rep(0:14, rep(15,15))
graph_attr(g1, 'xy') = data.frame(x,y) #define 'xy' attributes of graph g which represents coordinates
g1.ad_mat = as.matrix(g1, matrix.type=c('adjacency')) #get unweighted adjacency matrix of graph g
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
lattice01$f.true <- sin(2*pi*lattice01$xy$x/7) # true sine signal
Y.lattice01 <- t(replicate(N, lattice01$f.true))
outlier.ind.mat <- matrix(0, nrow=N, ncol=round(n.lattice01/10))
set.seed(10)
for(i in 1:N){
  outlier.ind <- sort(sample(n.lattice01, round(n.lattice01/10)))
  outlier.ind.mat[i,] <- outlier.ind
  Y.lattice01[i,outlier.ind] <- Y.lattice01[i,outlier.ind] + rnorm(round(n.lattice01/10),0,4) # generate 10% outliers
  Y.lattice01[i,-outlier.ind] <- Y.lattice01[i,-outlier.ind] + rnorm(n.lattice01-round(n.lattice01/10),0,0.2) 
}
plot_graph_custom3(lattice01, e.size=1.3, v.size=6, vertex_color = Y.lattice01[100,], 
                   min=min(Y.lattice01[100,]), max=max(Y.lattice01[100,]), value="value") + ggtitle("original") # plot graph signal
plot(Y.lattice01[100,1:15], type="l")

# visualize
plot_graph(lattice01) # plot graph structure
plot_graph_custom3(lattice01, e.size=1.3, v.size=6, vertex_color = lattice01$f.true, 
                   min=min(lattice01$f.true), max=max(lattice01$f.true), value="value") + ggtitle("original") # plot graph signal

MSE.mean.lattice01 <- NULL
MSE.imputedCV.mean.lattice01 <- NULL
MSE.tps.lattice01 <- NULL
MSE.SIC50.lattice01 <- NULL
MSE.GACV50.lattice01 <- NULL
MSE.imputedCV50.lattice01 <- NULL
opt.lambda.SIC50.lattice01.list <- NULL
opt.lambda.GACV50.lattice01.list <- NULL
opt.lambda.imputedCV50.lattice01.list <- NULL
opt.lambda.mean.lattice01.list <- NULL
opt.lambda.imputedCV.mean.lattice01.list <- NULL
opt.lambda.tps.lattice01.list <- NULL
for(i in 1:N){
  lattice01$f <- Y.lattice01[i,]
  # 50% quantile-based
  SIC_list50.lattice01 <- c()
  GACV_list50.lattice01 <- c()
  imputed_CV_list50.lattice01 <- c()
  H_list_qt50.lattice01 <- list()
  
  lambda_candidate50.lattice01 <- exp(seq(-9, -5, by=0.3))
  for(lambda in lambda_candidate50.lattice01){
    cat("####### lambda = ", lambda,"iteration", i, "####### \n")
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
  # plot(lambda_candidate50.lattice01, SIC_list50.lattice01, type="l")
  # plot(lambda_candidate50.lattice01, GACV_list50.lattice01, type="l")
  # plot(lambda_candidate50.lattice01, imputed_CV_list50.lattice01, type="l")
  opt.lambda.ind.SIC50.lattice01 <- which.min(SIC_list50.lattice01)
  opt.lambda.SIC50.lattice01 <- lambda_candidate50.lattice01[opt.lambda.ind.SIC50.lattice01] 
  opt.lambda.ind.GACV50.lattice01 <- which.min(GACV_list50.lattice01)
  opt.lambda.GACV50.lattice01 <- lambda_candidate50.lattice01[opt.lambda.ind.GACV50.lattice01] 
  opt.lambda.ind.imputedCV50.lattice01 <- which.min(imputed_CV_list50.lattice01)
  opt.lambda.imputedCV50.lattice01 <- lambda_candidate50.lattice01[opt.lambda.ind.imputedCV50.lattice01] 
  
  # mean-based
  H_list_mean.lattice01 <- list()
  GCV_list.lattice01 <- c()
  imputed_CV_list.mean.lattice01 <- c()
  lambda_candidate.mean.lattice01 <- exp(seq(-8, -4, by=0.3))
  for(lambda in lambda_candidate.mean.lattice01){
    H <- solve(diag(1,n.lattice01)+n.lattice01*lambda*L.lattice01%*%L.lattice01)
    H_list_mean.lattice01[[length(H_list_mean.lattice01)+1]] <- H
    GCV_list.lattice01 <- c(GCV_list.lattice01, GCV(lambda, lattice01$f, H))
    imputed_CV <- imputedCV(lambda=lambda, y=lattice01$f, tau=0.5, alpha=0.1, B=diag(n.lattice01), 
                            P=L.lattice01%*%L.lattice01, gamma_init=rep(mean(lattice01$f), n.lattice01), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.lattice01, k=1, M=NULL, method="mean")
    imputed_CV_list.mean.lattice01 <- c(imputed_CV_list.mean.lattice01, imputed_CV)
  }
  # plot(lambda_candidate.mean.lattice01, GCV_list.lattice01, type="l")
  # plot(lambda_candidate.mean.lattice01, imputed_CV_list.mean.lattice01, type="l")
  opt.lambda.ind.mean.lattice01 <- which.min(GCV_list.lattice01)
  opt.lambda.mean.lattice01 <- lambda_candidate.mean.lattice01[[opt.lambda.ind.mean.lattice01]] 
  opt.lambda.ind.imputedCV.mean.lattice01 <- which.min(imputed_CV_list.mean.lattice01)
  opt.lambda.imputedCV.mean.lattice01 <- lambda_candidate.mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]] 
  
  # thin plate spline
  GCV.tps.lattice01 <- sapply(lambda_candidate.mean.lattice01, function(l){
    tmp <- Tps(lattice01$xy, lattice01$f, lambda=l)
    return(unlist(summary(tmp))$sum.gcv.lambda.GCV)})
  # plot(lambda_candidate.mean.lattice01, GCV.tps.lattice01, type="l")
  opt.lambda.ind.tps.lattice01 <- which.min(GCV.tps.lattice01)
  opt.lambda.tps.lattice01 <- lambda_candidate.mean.lattice01[opt.lambda.ind.tps.lattice01]
  fit.tps.lattice01 <- Tps(lattice01$xy, lattice01$f, lambda = opt.lambda.tps.lattice01)
  
  # calculate MSE
  MSE.mean.lattice01 <- c(MSE.mean.lattice01, mse(lattice01$f.true, H_list_mean.lattice01[[opt.lambda.ind.mean.lattice01]]%*%lattice01$f))
  MSE.imputedCV.mean.lattice01 <- c(MSE.imputedCV.mean.lattice01, mse(lattice01$f.true, H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f))
  MSE.tps.lattice01 <- c(MSE.tps.lattice01, mse(lattice01$f.true, fit.tps.lattice01$fitted.values))
  MSE.SIC50.lattice01 <- c(MSE.SIC50.lattice01, mse(lattice01$f.true, H_list_qt50.lattice01[[opt.lambda.ind.SIC50.lattice01]]%*%lattice01$f))
  MSE.GACV50.lattice01 <- c(MSE.GACV50.lattice01, mse(lattice01$f.true, H_list_qt50.lattice01[[opt.lambda.ind.GACV50.lattice01]]%*%lattice01$f))
  MSE.imputedCV50.lattice01 <- c(MSE.imputedCV50.lattice01, mse(lattice01$f.true, H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f))
  
  # save optimal lambda
  opt.lambda.SIC50.lattice01.list <- c(opt.lambda.SIC50.lattice01.list, opt.lambda.SIC50.lattice01)
  opt.lambda.GACV50.lattice01.list <- c(opt.lambda.GACV50.lattice01.list, opt.lambda.GACV50.lattice01)
  opt.lambda.imputedCV50.lattice01.list <- c(opt.lambda.imputedCV50.lattice01.list, opt.lambda.imputedCV50.lattice01)
  opt.lambda.mean.lattice01.list <- c(opt.lambda.mean.lattice01.list, opt.lambda.mean.lattice01)
  opt.lambda.imputedCV.mean.lattice01.list <- c(opt.lambda.imputedCV.mean.lattice01.list, opt.lambda.imputedCV.mean.lattice01)
  opt.lambda.tps.lattice01.list <- c(opt.lambda.tps.lattice01.list, opt.lambda.tps.lattice01)
}

# scores for selecting lambda
par(mfrow=c(2,3))
plot(lambda_candidate50.lattice01, SIC_list50.lattice01, type="l", main = "SIC (tau=0.5)", cex.main=2)
points(opt.lambda.SIC50.lattice01, SIC_list50.lattice01[opt.lambda.ind.SIC50.lattice01], col="red", pch=19, cex=2)
plot(lambda_candidate50.lattice01, GACV_list50.lattice01, type="l", main = "GACV (tau=0.5)", cex.main=2)
points(opt.lambda.GACV50.lattice01, GACV_list50.lattice01[opt.lambda.ind.GACV50.lattice01], col="red", pch=19, cex=2)
plot(lambda_candidate50.lattice01, imputed_CV_list50.lattice01, type="l", main = "Imputed CV (tau=0.5)", cex.main=2)
points(opt.lambda.imputedCV50.lattice01, imputed_CV_list50.lattice01[opt.lambda.ind.imputedCV50.lattice01], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.lattice01, GCV_list.lattice01, type="l", main = "GCV (mean)", cex.main=2)
points(opt.lambda.mean.lattice01, GCV_list.lattice01[opt.lambda.ind.mean.lattice01], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.lattice01, imputed_CV_list.mean.lattice01, type="l", main = "Imputed CV (mean)", cex.main=2)
points(opt.lambda.imputedCV.mean.lattice01, imputed_CV_list.mean.lattice01[opt.lambda.ind.imputedCV.mean.lattice01], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.lattice01, GCV.tps.lattice01, type="l", main = "TPS", cex.main=2)
points(opt.lambda.tps.lattice01, GCV.tps.lattice01[opt.lambda.ind.tps.lattice01], col="red", pch=19, cex=2)


# denoising results
q0 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = lattice01$f.true, 
                         min=-2.2, max=2.2, value="value", ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("true signal")
q1 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = lattice01$f, 
                         min=min(lattice01$f), max=max(lattice01$f), value="value", ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("noisy signal")
q2 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f),
                         min=-2.2, max=2.2, 
                         value="value", ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("mean")
q3 <- plot_graph_custom3(lattice01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f),
                         min=-2.2, max=2.2, 
                         value="value", ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.5")

grid.arrange(q0,q1,q2,q3, nrow=1)

# specific row
par(mfrow=c(1,1), mar=c(5,5.5,4,2)+0.1, oma=c(0,0,0,0))
plot(lattice01$f[46:60],type="l", col="black", lty=2, ylab="value",
     xlab="node column index", cex.lab=3.2, cex.axis=2.5, cex=2, lwd=2)
lines(lattice01$f.true[46:60], lwd=2)
lines(as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f)[46:60], col="red", lwd=2)
lines(as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f)[46:60], col="blue", lwd=2)
legend(x = "topleft", text.font = 4, cex=2, 
       col= c("black","black", "blue","red"), lty=c(1,2,1,1),text.col = "black", 
       legend=c("true","noisy","mean","tau = 0.5"), lwd=2)

# scaled persp plots
par(mfrow=c(1,3))
persp(1:15, 1:15, matrix(lattice01$f, nrow=15), theta = -40, phi = 20, expand = 0.5, col = "gray",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "noisy signal", cex.main = 4, line=-16)
persp(1:15, 1:15, matrix(as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightblue",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "mean (GICV)", cex.main = 4, line=-16)
persp(1:15, 1:15, matrix(as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "tau=0.5 (GIQCV)", cex.main = 4, line=-16)

par(mfrow=c(1,3))
persp(1:15, 1:15, matrix(lattice01$f, nrow=15), theta = -40, phi = 20, expand = 0.5, col = "gray",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "noisy signal", cex.main = 4, line=-16)
persp(1:15, 1:15, matrix(as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightblue",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "mean", cex.main = 4, line=-16)
persp(1:15, 1:15, matrix(as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f), nrow=15),
      theta = -40, phi = 20, expand = 0.5, col = "lightcoral",
      xlab="column", ylab="row", zlab="value", zlim=range(matrix(lattice01$f, nrow=15)),
      cex.lab=3)
title(main = "tau=0.5", cex.main = 4, line=-16)

# non-scaled persp plots
par(mfrow=c(1,3))
persp(1:15, 1:15, t(matrix(lattice01$f, nrow=15))[15:1,], theta = 50, phi = 20, expand = 0.5, col = "gray",
      xlab="row", ylab="column", zlab="value", cex.lab=2.5)
title(main = "noisy signal", cex.main = 3, line=-10)
persp(1:15, 1:15, t(matrix(as.numeric(H_list_mean.lattice01[[opt.lambda.ind.imputedCV.mean.lattice01]]%*%lattice01$f), nrow=15))[15:1,],
      theta = 50, phi = 20, expand = 0.5, col = "lightblue",
      xlab="row", ylab="column", zlab="value", cex.lab=2.5)
title(main = "mean denoised (GICV)", cex.main = 3, line=-10)
persp(1:15, 1:15, t(matrix(as.numeric(H_list_qt50.lattice01[[opt.lambda.ind.imputedCV50.lattice01]]%*%lattice01$f), nrow=15))[15:1,],
      theta = 50, phi = 20, expand = 0.5, col = "lightcoral",
      xlab="row", ylab="column", zlab="value", cex.lab=2.5)
title(main = "50% quantile (GIQCV)", cex.main = 3, line=-10)

# MSE plot
par(mar=c(5,5.5,4,2)+0.1)
plot(MSE.imputedCV.mean.lattice01, type="l", ylim=c(0,0.8),
     ylab="MSE", xlab="replication index", col="black", lwd=2, cex.lab=3.2, cex.axis=2.5, cex=2)
lines(MSE.mean.lattice01, col="blue", lwd=2)
lines(MSE.imputedCV50.lattice01, col="red", lwd=2)
lines(MSE.SIC50.lattice01, col="magenta", lwd=2)
lines(MSE.GACV50.lattice01, col="cyan", lwd=2)

legend(x = "topright", lty = c(1,1,1,1,1), text.font = 4, cex=2, 
       col= c("black","blue","red","magenta", "cyan"), text.col = "black", 
       legend=c("GICV","GCV","GIQCV", "SIC","GACV"), lwd=2)

par(mar=c(5,5.5,4,2)+0.1)
plot(MSE.imputedCV.mean.lattice01, type="l", ylim=c(0,0.65),
     ylab="MSE", xlab="replication index", col="black", lwd=2, cex.lab=3.2, cex.axis=2.5, cex=2)
points(MSE.imputedCV.mean.lattice01, pch=19)
lines(MSE.imputedCV50.lattice01, col="red", lwd=2)
points(MSE.imputedCV50.lattice01, col="red", pch=19)
legend(x = "topleft", lty = c(1,1), text.font = 4, cex=2, 
       col= c("black","red"), text.col = "black", 
       legend=c("mean", "tau = 0.5"), lwd=2)


#######################################
### 03. irregular graph w/ outliers   
#######################################
x = rep(0:14, 15) #define coordinates of vertices
y = rep(0:14, rep(15,15))

irregular01 <- list()
set.seed(1)
irregular01$xy <- data.frame(x,y) + rnorm(225,0,0.7)# define 'xy' attributes of graph g which represents coordinates
n.irregular01 <- nrow(irregular01$xy)
rownames(irregular01$xy) <- 1:n.irregular01

distmat.irregular01 <- as.matrix(dist(irregular01$xy))
A.irregular01 <- c()
for(i in 1:(nrow(distmat.irregular01)-1)){
  for(j in (i+1):ncol(distmat.irregular01)){
    val <- distmat.irregular01[i,j]
    A.irregular01 <- rbind(A.irregular01, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.irregular01, k=5), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.irregular01[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=n.irregular01, ncol=n.irregular01)

colnames(wmat) <- 1:n.irregular01
rownames(wmat) <- 1:n.irregular01

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}


# sparse weight matrix
# weight matrix
irregular01$A <- wmat

# sparse weight matrix
irregular01$sA <- sp.wmat

irregular01$dist <- distmat.irregular01
irregular01$sdist <- A.irregular01

L.irregular01 <- gasper::laplacian_mat(irregular01$A)
eigenres.irregular01 <- eigen(L.irregular01)


plot_graph(irregular01)
plot_graph_custom3(irregular01, e.size=1.3, v.size=6, vertex_color = eigenres.irregular01$vectors[,50],
                   min=-0.5, max=0.5, value="value")


N <- 100

# signal assignment
irregular01$f.true <- exp(sqrt((irregular01$xy$x-7)^2 + (irregular01$xy$y-7)^2)/5) # true sine signal
plot_graph_custom3(irregular01, e.size=1.3, v.size=6, vertex_color = irregular01$f.true,
                   min=min(irregular01$f.true), max=max(irregular01$f.true), value="value")

Y.irregular01 <- t(replicate(N, irregular01$f.true))
outlier.ind.mat <- matrix(0, nrow=N, ncol=round(n.irregular01/10))
set.seed(100)
for(i in 1:N){
  outlier.ind <- sort(sample(n.irregular01, round(n.irregular01/10)))
  outlier.ind.mat[i,] <- outlier.ind
  Y.irregular01[i,outlier.ind] <- Y.irregular01[i,outlier.ind] + rnorm(round(n.irregular01/10),0,5) # generate 10% outliers
  Y.irregular01[i,-outlier.ind] <- Y.irregular01[i,-outlier.ind] + rnorm(n.irregular01-round(n.irregular01/10),0,0.2) 
}
plot_graph_custom3(irregular01, e.size=1.3, v.size=6, vertex_color = Y.irregular01[100,], 
                   min=min(Y.irregular01[100,]), max=max(Y.irregular01[100,]), value="value") + ggtitle("original") # plot graph signal

irregular01$f <- Y.irregular01[100,]
# visualize
plot_graph(irregular01) # plot graph structure
plot_graph_custom3(irregular01, e.size=1.3, v.size=6, vertex_color = irregular01$f.true, 
                   min=min(irregular01$f.true), max=max(irregular01$f.true), value="value") + ggtitle("original") # plot graph signal

MSE.mean.irregular01 <- NULL
MSE.imputedCV.mean.irregular01 <- NULL
MSE.tps.irregular01 <- NULL
MSE.SIC50.irregular01 <- NULL
MSE.GACV50.irregular01 <- NULL
MSE.imputedCV50.irregular01 <- NULL 
opt.lambda.SIC50.irregular01.list <- NULL
opt.lambda.GACV50.irregular01.list <- NULL
opt.lambda.imputedCV50.irregular01.list <- NULL
opt.lambda.mean.irregular01.list <- NULL
opt.lambda.imputedCV.mean.irregular01.list <- NULL
opt.lambda.tps.irregular01.list <- NULL
for(i in 1:N){
  irregular01$f <- Y.irregular01[i,]
  # 50% quantile-based
  SIC_list50.irregular01 <- c()
  GACV_list50.irregular01 <- c()
  imputed_CV_list50.irregular01 <- c()
  H_list_qt50.irregular01 <- list()
  
  lambda_candidate50.irregular01 <- exp(seq(-6, floor(log(n.irregular01))-1, by=0.3))
  for(lambda in lambda_candidate50.irregular01){
    res50.irregular01 <- PIRLS(lambda=lambda, y=irregular01$f, tau=0.5, alpha=0.1, B=diag(n.irregular01), 
                               P=L.irregular01%*%L.irregular01, gamma_init=rep(mean(irregular01$f), n.irregular01), max_iterations=10000, tol=1e-6, check="Oh")
    SIC_list50.irregular01 <- c(SIC_list50.irregular01, res50.irregular01$SIC)
    GACV_list50.irregular01 <- c(GACV_list50.irregular01, res50.irregular01$GACV)
    imputed_CV <- imputedCV(lambda=lambda, y=irregular01$f, tau=0.5, alpha=0.1, B=diag(n.irregular01), 
                            P=L.irregular01%*%L.irregular01, gamma_init=rep(mean(irregular01$f), n.irregular01), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.irregular01, k=1, M=NULL, method="quantile")
    imputed_CV_list50.irregular01 <- c(imputed_CV_list50.irregular01, imputed_CV)
    H_list_qt50.irregular01[[length(H_list_qt50.irregular01)+1]] <- res50.irregular01$H
  }
  # plot(lambda_candidate50.irregular01, SIC_list50.irregular01, type="l")
  # plot(lambda_candidate50.irregular01, GACV_list50.irregular01, type="l")
  # plot(lambda_candidate50.irregular01, imputed_CV_list50.irregular01, type="l")
  opt.lambda.ind.SIC50.irregular01 <- which.min(SIC_list50.irregular01)
  opt.lambda.SIC50.irregular01 <- lambda_candidate50.irregular01[opt.lambda.ind.SIC50.irregular01] 
  opt.lambda.ind.GACV50.irregular01 <- which.min(GACV_list50.irregular01)
  opt.lambda.GACV50.irregular01 <- lambda_candidate50.irregular01[opt.lambda.ind.GACV50.irregular01] 
  opt.lambda.ind.imputedCV50.irregular01 <- which.min(imputed_CV_list50.irregular01)
  opt.lambda.imputedCV50.irregular01 <- lambda_candidate50.irregular01[opt.lambda.ind.imputedCV50.irregular01] 
  
  # mean-based
  H_list_mean.irregular01 <- list()
  GCV_list.irregular01 <- c()
  imputed_CV_list.mean.irregular01 <- c()
  lambda_candidate.mean.irregular01 <- exp(seq(-6, floor(log(n.irregular01))-1, by=0.3))
  for(lambda in lambda_candidate.mean.irregular01){
    H <- solve(diag(1,n.irregular01)+n.irregular01*lambda*L.irregular01%*%L.irregular01)
    H_list_mean.irregular01[[length(H_list_mean.irregular01)+1]] <- H
    GCV_list.irregular01 <- c(GCV_list.irregular01, GCV(lambda, irregular01$f, H))
    imputed_CV <- imputedCV(lambda=lambda, y=irregular01$f, tau=0.5, alpha=0.1, B=diag(n.irregular01), 
                            P=L.irregular01%*%L.irregular01, gamma_init=rep(mean(irregular01$f), n.irregular01), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.irregular01, k=1, M=NULL, method="mean")
    imputed_CV_list.mean.irregular01 <- c(imputed_CV_list.mean.irregular01, imputed_CV)
  }
  # plot(lambda_candidate.mean.irregular01, GCV_list.irregular01, type="l")
  # plot(lambda_candidate.mean.irregular01, imputed_CV_list.mean.irregular01, type="l")
  opt.lambda.ind.mean.irregular01 <- which.min(GCV_list.irregular01)
  opt.lambda.mean.irregular01 <- lambda_candidate.mean.irregular01[opt.lambda.ind.mean.irregular01] 
  opt.lambda.ind.imputedCV.mean.irregular01 <- which.min(imputed_CV_list.mean.irregular01)
  opt.lambda.imputedCV.mean.irregular01 <- lambda_candidate.mean.irregular01[opt.lambda.ind.imputedCV.mean.irregular01] 
  
  # thin plate spline
  GCV.tps.irregular01 <- sapply(lambda_candidate.mean.irregular01, function(l){
    tmp <- Tps(irregular01$xy, irregular01$f, lambda=l)
    return(unlist(summary(tmp))$sum.gcv.lambda.GCV)})
  # plot(lambda_candidate.mean.irregular01, GCV.tps.irregular01, type="l")
  opt.lambda.ind.tps.irregular01 <- which.min(GCV.tps.irregular01)
  opt.lambda.tps.irregular01 <- lambda_candidate.mean.irregular01[opt.lambda.ind.tps.irregular01]
  fit.tps.irregular01 <- Tps(irregular01$xy, irregular01$f, lambda = opt.lambda.tps.irregular01)
  
  # calculate MSE
  MSE.mean.irregular01 <- c(MSE.mean.irregular01, mse(irregular01$f.true, H_list_mean.irregular01[[opt.lambda.ind.mean.irregular01]]%*%irregular01$f))
  MSE.imputedCV.mean.irregular01 <- c(MSE.imputedCV.mean.irregular01, mse(irregular01$f.true, H_list_mean.irregular01[[opt.lambda.ind.imputedCV.mean.irregular01]]%*%irregular01$f))
  MSE.tps.irregular01 <- c(MSE.tps.irregular01, mse(irregular01$f.true, fit.tps.irregular01$fitted.values))
  MSE.SIC50.irregular01 <- c(MSE.SIC50.irregular01, mse(irregular01$f.true, H_list_qt50.irregular01[[opt.lambda.ind.SIC50.irregular01]]%*%irregular01$f))
  MSE.GACV50.irregular01 <- c(MSE.GACV50.irregular01, mse(irregular01$f.true, H_list_qt50.irregular01[[opt.lambda.ind.GACV50.irregular01]]%*%irregular01$f))
  MSE.imputedCV50.irregular01 <- c(MSE.imputedCV50.irregular01, mse(irregular01$f.true, H_list_qt50.irregular01[[opt.lambda.ind.imputedCV50.irregular01]]%*%irregular01$f))
  
  # save optimal lambda
  opt.lambda.SIC50.irregular01.list <- c(opt.lambda.SIC50.irregular01.list, opt.lambda.SIC50.irregular01)
  opt.lambda.GACV50.irregular01.list <- c(opt.lambda.GACV50.irregular01.list, opt.lambda.GACV50.irregular01)
  opt.lambda.imputedCV50.irregular01.list <- c(opt.lambda.imputedCV50.irregular01.list, opt.lambda.imputedCV50.irregular01)
  opt.lambda.mean.irregular01.list <- c(opt.lambda.mean.irregular01.list, opt.lambda.mean.irregular01)
  opt.lambda.imputedCV.mean.irregular01.list <- c(opt.lambda.imputedCV.mean.irregular01.list, opt.lambda.imputedCV.mean.irregular01)
  opt.lambda.tps.irregular01.list <- c(opt.lambda.tps.irregular01.list, opt.lambda.tps.irregular01)
}

# scores for selecting lambda
par(mfrow=c(2,3))
plot(lambda_candidate50.irregular01, SIC_list50.irregular01, type="l", main = "SIC (tau=0.5)", cex.main=2)
points(opt.lambda.SIC50.irregular01, SIC_list50.irregular01[opt.lambda.ind.SIC50.irregular01], col="red", pch=19, cex=2)
plot(lambda_candidate50.irregular01, GACV_list50.irregular01, type="l", main = "GACV (tau=0.5)", cex.main=2)
points(opt.lambda.GACV50.irregular01, GACV_list50.irregular01[opt.lambda.ind.GACV50.irregular01], col="red", pch=19, cex=2)
plot(lambda_candidate50.irregular01, imputed_CV_list50.irregular01, type="l", main = "Imputed CV (tau=0.5)", cex.main=2)
points(opt.lambda.imputedCV50.irregular01, imputed_CV_list50.irregular01[opt.lambda.ind.imputedCV50.irregular01], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.irregular01, GCV_list.irregular01, type="l", main = "GCV (mean)", cex.main=2)
points(opt.lambda.mean.irregular01, GCV_list.irregular01[opt.lambda.ind.mean.irregular01], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.irregular01, imputed_CV_list.mean.irregular01, type="l", main = "Imputed CV (mean)", cex.main=2)
points(opt.lambda.imputedCV.mean.irregular01, imputed_CV_list.mean.irregular01[opt.lambda.ind.imputedCV.mean.irregular01], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.irregular01, GCV.tps.irregular01, type="l", main = "TPS", cex.main=2)
points(opt.lambda.tps.irregular01, GCV.tps.irregular01[opt.lambda.ind.tps.irregular01], col="red", pch=19, cex=2)

# denoising results
q0 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=11, vertex_color = irregular01$f.true,
                         min=1.1, max=9.1, value="value", label.title.size=19, label.text.size = 17) + ggtitle("true signal")

q1 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=11, vertex_color = irregular01$f,
                         min=min(irregular01$f), max=max(irregular01$f), value="value", label.title.size=19, label.text.size = 17) + ggtitle("noisy signal")

q2 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_mean.irregular01[[opt.lambda.ind.imputedCV.mean.irregular01]]%*%irregular01$f),
                         min=1.1, max=9.1, value="value", label.title.size=19, label.text.size = 17) + ggtitle("mean")

q3 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt50.irregular01[[opt.lambda.ind.imputedCV50.irregular01]]%*%irregular01$f),
                         min=1.1, max=9.1, value="value", label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.5")

grid.arrange(q0,q1,q2,q3, nrow=1)


# MSE plot
par(mar=c(5,5.5,4,2)+0.1)
plot(MSE.imputedCV.mean.irregular01, type="l", ylim=c(0,3),
     ylab="MSE", xlab="replication index", col="black", lwd=2, cex.lab=3.2, cex.axis=2.5, cex=2)
lines(MSE.mean.irregular01, col="blue", lwd=2)
lines(MSE.imputedCV50.irregular01, col="red", lwd=2)
lines(MSE.SIC50.irregular01, col="magenta", lwd=2)
lines(MSE.GACV50.irregular01, col="cyan", lwd=2)

legend(x = "topright", lty = c(1,1,1,1,1,1), text.font = 4, cex=2, 
       col= c("black","blue", "red","magenta", "cyan"), text.col = "black", 
       legend=c("GICV","GCV","GIQCV", "SIC","GACV"), lwd=2)

par(mar=c(5,5.5,4,2)+0.1)
plot(MSE.imputedCV.mean.irregular01, type="l", ylim=c(0,3),
     ylab="MSE", xlab="replication index", col="black", lwd=2, cex.lab=3.2, cex.axis=2.5, cex=2)
points(MSE.imputedCV.mean.irregular01, pch=19)
lines(MSE.imputedCV50.irregular01, col="red", lwd=2)
points(MSE.imputedCV50.irregular01, col="red", pch=19)
legend(x = "topleft", lty = c(1,1), text.font = 4, cex=2, 
       col= c("black","red"), text.col = "black", 
       legend=c("mean", "tau = 0.5"), lwd=2)

