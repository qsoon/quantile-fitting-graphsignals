# Required Libraries
library(gasper)
library(igraph)
library(ggplot2)
library(optimx)
library(Matrix)
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
library(quantreg)

source("utils.R")

###################################
### 01. Manhattan taxi data
###################################
data("NYCdata") # A: NYC adjacency matrix, f: Total Amount observed from miles travelled to the given drop off point

manhattanzone <- read.csv("Manhattan_zones.csv") # Manhattan zones info in NYC
manhattan.zone.ind <- sort(unique(manhattanzone$LocationID))
W.manhattan <- NYCdata$A[manhattan.zone.ind, manhattan.zone.ind]
rowSums(W.manhattan) # remove zero places (disconnected nodes)
manhattan.zone.ind <- manhattan.zone.ind[manhattan.zone.ind!=103]
W.manhattan <- NYCdata$A[manhattan.zone.ind, manhattan.zone.ind] # weight matrix of Manhattan taxi data
manhattanzone <- manhattanzone[manhattanzone$LocationID %in% manhattan.zone.ind,]
manhattanzone <- manhattanzone[order(manhattanzone$LocationID),]

M.taxi <- list()
M.taxi$A <- W.manhattan
rownames(M.taxi$A) <- manhattan.zone.ind
colnames(M.taxi$A) <- manhattan.zone.ind
M.taxi$f <- NYCdata$f[manhattan.zone.ind]

avg_coords <- function(poly_string) {
  # Extract all numeric values (including negative signs)
  numbers <- as.numeric(unlist(regmatches(poly_string, gregexpr("(-?\\d+(\\.\\d+)?)", poly_string))))
  
  # Extract x and y coordinates
  x_coords <- numbers[seq(1, length(numbers), 2)]
  y_coords <- numbers[seq(2, length(numbers), 2)]
  
  # Calculate averages
  list(avg_x = mean(x_coords), avg_y = mean(y_coords))
}

results <- lapply(manhattanzone$the_geom, avg_coords)
avg_x <- sapply(results, `[[`, "avg_x")
avg_y <- sapply(results, `[[`, "avg_y")
M.taxi$xy <- cbind(avg_x, avg_y)
rownames(M.taxi$xy) <- manhattan.zone.ind
colnames(M.taxi$xy) <- c("x","y")

sp.wmat <- c()
for(i in 1:(length(manhattan.zone.ind)-1)){
  for(j in (i+1):length(manhattan.zone.ind)){
    if(M.taxi$A[i,j]!=0){
      sp.wmat <- rbind(sp.wmat, c(i, j, M.taxi$A[i,j]))
    }
  }
}
M.taxi$sA <- sp.wmat

# reduce edges
M.taxi.reduced <- M.taxi
M.taxi.reduced$A <- M.taxi.reduced$A*(M.taxi.reduced$A >= 0.94) 
tmp.ind <- which(rowSums(M.taxi.reduced$A)==0)
M.taxi.reduced$A <- M.taxi.reduced$A[-tmp.ind,-tmp.ind]
sp.wmat <- c()
for(i in 1:(length(manhattan.zone.ind[-tmp.ind])-1)){
  for(j in (i+1):length(manhattan.zone.ind[-tmp.ind])){
    if(M.taxi.reduced$A[i,j]!=0){
      sp.wmat <- rbind(sp.wmat, c(i, j, M.taxi.reduced$A[i,j]))
    }
  }
}
M.taxi.reduced$xy <- M.taxi.reduced$xy[-tmp.ind,] 
M.taxi.reduced$sA <- sp.wmat
M.taxi.reduced$f <- M.taxi.reduced$f[-tmp.ind]
L.manhattan.reduced <- gasper::laplacian_mat(M.taxi.reduced$A)
n.manhattan.reduced <- nrow(L.manhattan.reduced)

# visualize
plot_graph(M.taxi.reduced)
plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=6, vertex_color = M.taxi.reduced$f, 
                   min=min(M.taxi.reduced$f), max=max(M.taxi.reduced$f), value="Total Amount")

## mean-based
H_list_mean.manhattan.reduced <- list()
GCV_list.manhattan.reduced <- c()
imputed_CV_list.mean.manhattan.reduced <- c()
lambda_candidate.mean.manhattan.reduced <- exp(seq(-30, -5, by=0.3)) 

for(lambda in lambda_candidate.mean.manhattan.reduced){
  H <- solve(diag(1,n.manhattan.reduced)+n.manhattan.reduced*lambda*L.manhattan.reduced%*%L.manhattan.reduced)
  H_list_mean.manhattan.reduced[[length(H_list_mean.manhattan.reduced)+1]] <- H
  GCV_list.manhattan.reduced <- c(GCV_list.manhattan.reduced, GCV(lambda, M.taxi.reduced$f, H))
  imputed_CV <- imputedCV(lambda=lambda, y=M.taxi.reduced$f, tau=0.5, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.manhattan.reduced, k=1, M=NULL, method="mean")
  
  imputed_CV_list.mean.manhattan.reduced <- c(imputed_CV_list.mean.manhattan.reduced, imputed_CV)
}
plot(lambda_candidate.mean.manhattan.reduced, GCV_list.manhattan.reduced, type="l")
plot(lambda_candidate.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced, type="l")
opt.lambda.ind.mean.manhattan.reduced <- which.min(GCV_list.manhattan.reduced)
opt.lambda.mean.manhattan.reduced <- lambda_candidate.mean.manhattan.reduced[opt.lambda.ind.mean.manhattan.reduced] 
opt.lambda.ind.imputedCV.mean.manhattan.reduced <- which.min(imputed_CV_list.mean.manhattan.reduced)
opt.lambda.imputedCV.mean.manhattan.reduced <- lambda_candidate.mean.manhattan.reduced[opt.lambda.ind.imputedCV.mean.manhattan.reduced] 

# thin plate spline
par(mfrow=c(1,1))
GCV.tps.manhattan.reduced <- sapply(lambda_candidate.mean.manhattan.reduced, function(l){
  tmp <- Tps(M.taxi.reduced$xy, M.taxi.reduced$f, lambda=l)
  return(unlist(summary(tmp))$sum.gcv.lambda.GCV)})

plot(lambda_candidate.mean.manhattan.reduced, GCV.tps.manhattan.reduced, type="l")
opt.lambda.ind.tps.manhattan.reduced <- which.min(GCV.tps.manhattan.reduced)
opt.lambda.tps.manhattan.reduced <- lambda_candidate.mean.manhattan.reduced[opt.lambda.ind.tps.manhattan.reduced]

## quantile-based

# tau=0.9
SIC_list90.manhattan.reduced <- c()
GACV_list90.manhattan.reduced <- c()
imputed_CV_list90.manhattan.reduced <- c()
H_list_qt90.manhattan.reduced <- list()

lambda_candidate90.manhattan.reduced <- exp(seq(-27, -5, by=0.3))

for(lambda in lambda_candidate90.manhattan.reduced){
  res90.manhattan.reduced <- PIRLS(lambda=lambda, y=M.taxi.reduced$f, tau=0.9, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list90.manhattan.reduced <- c(SIC_list90.manhattan.reduced, res90.manhattan.reduced$SIC)
  GACV_list90.manhattan.reduced <- c(GACV_list90.manhattan.reduced, res90.manhattan.reduced$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=M.taxi.reduced$f, tau=0.9, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.manhattan.reduced, k=1, M=NULL, method="quantile")
  
  imputed_CV_list90.manhattan.reduced <- c(imputed_CV_list90.manhattan.reduced, imputed_CV)
  H_list_qt90.manhattan.reduced[[length(H_list_qt90.manhattan.reduced)+1]] <- res90.manhattan.reduced$H
}
plot(lambda_candidate90.manhattan.reduced, SIC_list90.manhattan.reduced, type="l")
plot(lambda_candidate90.manhattan.reduced, GACV_list90.manhattan.reduced, type="l")
plot(lambda_candidate90.manhattan.reduced, imputed_CV_list90.manhattan.reduced, type="l")
opt.lambda.ind.SIC90.manhattan.reduced <- which.min(SIC_list90.manhattan.reduced)
opt.lambda.SIC90.manhattan.reduced <- lambda_candidate90.manhattan.reduced[opt.lambda.ind.SIC90.manhattan.reduced] 
opt.lambda.ind.GACV90.manhattan.reduced <- which.min(GACV_list90.manhattan.reduced)
opt.lambda.GACV90.manhattan.reduced <- lambda_candidate90.manhattan.reduced[opt.lambda.ind.GACV90.manhattan.reduced] 
opt.lambda.ind.imputedCV90.manhattan.reduced <- which.min(imputed_CV_list90.manhattan.reduced)
opt.lambda.imputedCV90.manhattan.reduced <- lambda_candidate90.manhattan.reduced[opt.lambda.ind.imputedCV90.manhattan.reduced] 

# tau=0.75
SIC_list75.manhattan.reduced <- c()
GACV_list75.manhattan.reduced <- c()
imputed_CV_list75.manhattan.reduced <- c()
H_list_qt75.manhattan.reduced <- list()

lambda_candidate75.manhattan.reduced <- exp(seq(-28, -6, by=0.3))

for(lambda in lambda_candidate75.manhattan.reduced){
  res75.manhattan.reduced <- PIRLS(lambda=lambda, y=M.taxi.reduced$f, tau=0.75, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list75.manhattan.reduced <- c(SIC_list75.manhattan.reduced, res75.manhattan.reduced$SIC)
  GACV_list75.manhattan.reduced <- c(GACV_list75.manhattan.reduced, res75.manhattan.reduced$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=M.taxi.reduced$f, tau=0.75, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.manhattan.reduced, k=1, M=NULL, method="quantile")
  
  imputed_CV_list75.manhattan.reduced <- c(imputed_CV_list75.manhattan.reduced, imputed_CV)
  H_list_qt75.manhattan.reduced[[length(H_list_qt75.manhattan.reduced)+1]] <- res75.manhattan.reduced$H
}
plot(lambda_candidate75.manhattan.reduced, SIC_list75.manhattan.reduced, type="l")
plot(lambda_candidate75.manhattan.reduced, GACV_list75.manhattan.reduced, type="l")
plot(lambda_candidate75.manhattan.reduced, imputed_CV_list75.manhattan.reduced, type="l")
opt.lambda.ind.SIC75.manhattan.reduced <- which.min(SIC_list75.manhattan.reduced)
opt.lambda.SIC75.manhattan.reduced <- lambda_candidate75.manhattan.reduced[opt.lambda.ind.SIC75.manhattan.reduced]
opt.lambda.ind.GACV75.manhattan.reduced <- which.min(GACV_list75.manhattan.reduced)
opt.lambda.GACV75.manhattan.reduced <- lambda_candidate75.manhattan.reduced[opt.lambda.ind.GACV75.manhattan.reduced] 
opt.lambda.ind.imputedCV75.manhattan.reduced <- which.min(imputed_CV_list75.manhattan.reduced)
opt.lambda.imputedCV75.manhattan.reduced <- lambda_candidate75.manhattan.reduced[opt.lambda.ind.imputedCV75.manhattan.reduced] 

# tau=0.5
SIC_list50.manhattan.reduced <- c()
GACV_list50.manhattan.reduced <- c()
imputed_CV_list50.manhattan.reduced <- c()
H_list_qt50.manhattan.reduced <- list()

lambda_candidate50.manhattan.reduced <- exp(seq(-30, -8, by=0.3))

for(lambda in lambda_candidate50.manhattan.reduced){
  res50.manhattan.reduced <- PIRLS(lambda=lambda, y=M.taxi.reduced$f, tau=0.5, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list50.manhattan.reduced <- c(SIC_list50.manhattan.reduced, res50.manhattan.reduced$SIC)
  GACV_list50.manhattan.reduced <- c(GACV_list50.manhattan.reduced, res50.manhattan.reduced$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=M.taxi.reduced$f, tau=0.5, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.manhattan.reduced, k=1, M=NULL, method="quantile")
  
  imputed_CV_list50.manhattan.reduced <- c(imputed_CV_list50.manhattan.reduced, imputed_CV)
  H_list_qt50.manhattan.reduced[[length(H_list_qt50.manhattan.reduced)+1]] <- res50.manhattan.reduced$H
}
plot(lambda_candidate50.manhattan.reduced, SIC_list50.manhattan.reduced, type="l")
plot(lambda_candidate50.manhattan.reduced, GACV_list50.manhattan.reduced, type="l")
plot(lambda_candidate50.manhattan.reduced, imputed_CV_list50.manhattan.reduced, type="l")
opt.lambda.ind.SIC50.manhattan.reduced <- which.min(SIC_list50.manhattan.reduced)
opt.lambda.SIC50.manhattan.reduced <- lambda_candidate50.manhattan.reduced[opt.lambda.ind.SIC50.manhattan.reduced] 
opt.lambda.ind.GACV50.manhattan.reduced <- which.min(GACV_list50.manhattan.reduced)
opt.lambda.GACV50.manhattan.reduced <- lambda_candidate50.manhattan.reduced[opt.lambda.ind.GACV50.manhattan.reduced]
opt.lambda.ind.imputedCV50.manhattan.reduced <- which.min(imputed_CV_list50.manhattan.reduced)
opt.lambda.imputedCV50.manhattan.reduced <- lambda_candidate50.manhattan.reduced[opt.lambda.ind.imputedCV50.manhattan.reduced] 

# tau=0.25
SIC_list25.manhattan.reduced <- c()
GACV_list25.manhattan.reduced <- c()
imputed_CV_list25.manhattan.reduced <- c()
H_list_qt25.manhattan.reduced <- list()

lambda_candidate25.manhattan.reduced <- exp(seq(-27, -4, by=0.3))

for(lambda in lambda_candidate25.manhattan.reduced){
  res25.manhattan.reduced <- PIRLS(lambda=lambda, y=M.taxi.reduced$f, tau=0.25, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list25.manhattan.reduced <- c(SIC_list25.manhattan.reduced, res25.manhattan.reduced$SIC)
  GACV_list25.manhattan.reduced <- c(GACV_list25.manhattan.reduced, res25.manhattan.reduced$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=M.taxi.reduced$f, tau=0.25, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.manhattan.reduced, k=1, M=NULL, method="quantile")
  
  imputed_CV_list25.manhattan.reduced <- c(imputed_CV_list25.manhattan.reduced, imputed_CV)
  H_list_qt25.manhattan.reduced[[length(H_list_qt25.manhattan.reduced)+1]] <- res25.manhattan.reduced$H
}
plot(lambda_candidate25.manhattan.reduced, SIC_list25.manhattan.reduced, type="l")
plot(lambda_candidate25.manhattan.reduced, GACV_list25.manhattan.reduced, type="l")
plot(lambda_candidate25.manhattan.reduced, imputed_CV_list25.manhattan.reduced, type="l")
opt.lambda.ind.SIC25.manhattan.reduced <- which.min(SIC_list25.manhattan.reduced)
opt.lambda.SIC25.manhattan.reduced <- lambda_candidate25.manhattan.reduced[opt.lambda.ind.SIC25.manhattan.reduced]
opt.lambda.ind.GACV25.manhattan.reduced <- which.min(GACV_list25.manhattan.reduced)
opt.lambda.GACV25.manhattan.reduced <- lambda_candidate25.manhattan.reduced[opt.lambda.ind.GACV25.manhattan.reduced] 
opt.lambda.ind.imputedCV25.manhattan.reduced <- which.min(imputed_CV_list25.manhattan.reduced)
opt.lambda.imputedCV25.manhattan.reduced <- lambda_candidate25.manhattan.reduced[opt.lambda.ind.imputedCV25.manhattan.reduced]

# tau=0.1
SIC_list10.manhattan.reduced <- c()
GACV_list10.manhattan.reduced <- c()
imputed_CV_list10.manhattan.reduced <- c()
H_list_qt10.manhattan.reduced <- list()

lambda_candidate10.manhattan.reduced <- exp(seq(-27, -4, by=0.3))

for(lambda in lambda_candidate10.manhattan.reduced){
  res10.manhattan.reduced <- PIRLS(lambda=lambda, y=M.taxi.reduced$f, tau=0.1, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list10.manhattan.reduced <- c(SIC_list10.manhattan.reduced, res10.manhattan.reduced$SIC)
  GACV_list10.manhattan.reduced <- c(GACV_list10.manhattan.reduced, res10.manhattan.reduced$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=M.taxi.reduced$f, tau=0.1, alpha=0.1, B=diag(n.manhattan.reduced), 
                            P=L.manhattan.reduced%*%L.manhattan.reduced, gamma_init=rep(mean(M.taxi.reduced$f), n.manhattan.reduced), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.manhattan.reduced, k=1, M=NULL, method="quantile")
  
  imputed_CV_list10.manhattan.reduced <- c(imputed_CV_list10.manhattan.reduced, imputed_CV)
  H_list_qt10.manhattan.reduced[[length(H_list_qt10.manhattan.reduced)+1]] <- res10.manhattan.reduced$H
}
plot(lambda_candidate10.manhattan.reduced, SIC_list10.manhattan.reduced, type="l")
plot(lambda_candidate10.manhattan.reduced, GACV_list10.manhattan.reduced, type="l")
plot(lambda_candidate10.manhattan.reduced, imputed_CV_list10.manhattan.reduced, type="l")
opt.lambda.ind.SIC10.manhattan.reduced <- which.min(SIC_list10.manhattan.reduced)
opt.lambda.SIC10.manhattan.reduced <- lambda_candidate10.manhattan.reduced[opt.lambda.ind.SIC10.manhattan.reduced] 
opt.lambda.ind.GACV10.manhattan.reduced <- which.min(GACV_list10.manhattan.reduced)
opt.lambda.GACV10.manhattan.reduced <- lambda_candidate10.manhattan.reduced[opt.lambda.ind.GACV10.manhattan.reduced] 
opt.lambda.ind.imputedCV10.manhattan.reduced <- which.min(imputed_CV_list10.manhattan.reduced)
opt.lambda.imputedCV10.manhattan.reduced <- lambda_candidate10.manhattan.reduced[opt.lambda.ind.imputedCV10.manhattan.reduced] 


# lambda selection result visualization
par(mfrow=c(2,3))
plot(lambda_candidate90.manhattan.reduced, SIC_list90.manhattan.reduced, type="l", main = "SIC (tau=0.9)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC90.manhattan.reduced, SIC_list90.manhattan.reduced[opt.lambda.ind.SIC90.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate90.manhattan.reduced, GACV_list90.manhattan.reduced, type="l", main = "GACV (tau=0.9)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV90.manhattan.reduced, GACV_list90.manhattan.reduced[opt.lambda.ind.GACV90.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate90.manhattan.reduced, imputed_CV_list90.manhattan.reduced, type="l", main = "GIQCV (tau=0.9)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV90.manhattan.reduced, imputed_CV_list90.manhattan.reduced[opt.lambda.ind.imputedCV90.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV_list.manhattan.reduced, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.manhattan.reduced, GCV_list.manhattan.reduced[opt.lambda.ind.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced[opt.lambda.ind.imputedCV.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV.tps.manhattan.reduced, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.manhattan.reduced, GCV.tps.manhattan.reduced[opt.lambda.ind.tps.manhattan.reduced], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate75.manhattan.reduced, SIC_list75.manhattan.reduced, type="l", main = "SIC (tau=0.75)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC75.manhattan.reduced, SIC_list75.manhattan.reduced[opt.lambda.ind.SIC75.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate75.manhattan.reduced, GACV_list75.manhattan.reduced, type="l", main = "GACV (tau=0.75)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV75.manhattan.reduced, GACV_list75.manhattan.reduced[opt.lambda.ind.GACV75.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate75.manhattan.reduced, imputed_CV_list75.manhattan.reduced, type="l", main = "GIQCV (tau=0.75)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV75.manhattan.reduced, imputed_CV_list75.manhattan.reduced[opt.lambda.ind.imputedCV75.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV_list.manhattan.reduced, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.manhattan.reduced, GCV_list.manhattan.reduced[opt.lambda.ind.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced[opt.lambda.ind.imputedCV.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV.tps.manhattan.reduced, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.manhattan.reduced, GCV.tps.manhattan.reduced[opt.lambda.ind.tps.manhattan.reduced], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate50.manhattan.reduced, SIC_list50.manhattan.reduced, type="l", main = "SIC (tau=0.5)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC50.manhattan.reduced, SIC_list50.manhattan.reduced[opt.lambda.ind.SIC50.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate50.manhattan.reduced, GACV_list50.manhattan.reduced, type="l", main = "GACV (tau=0.5)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV50.manhattan.reduced, GACV_list50.manhattan.reduced[opt.lambda.ind.GACV50.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate50.manhattan.reduced, imputed_CV_list50.manhattan.reduced, type="l", main = "GIQCV (tau=0.5)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV50.manhattan.reduced, imputed_CV_list50.manhattan.reduced[opt.lambda.ind.imputedCV50.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV_list.manhattan.reduced, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.manhattan.reduced, GCV_list.manhattan.reduced[opt.lambda.ind.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced[opt.lambda.ind.imputedCV.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV.tps.manhattan.reduced, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.manhattan.reduced, GCV.tps.manhattan.reduced[opt.lambda.ind.tps.manhattan.reduced], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate25.manhattan.reduced, SIC_list25.manhattan.reduced, type="l", main = "SIC (tau=0.25)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC25.manhattan.reduced, SIC_list25.manhattan.reduced[opt.lambda.ind.SIC25.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate25.manhattan.reduced, GACV_list25.manhattan.reduced, type="l", main = "GACV (tau=0.25)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV25.manhattan.reduced, GACV_list25.manhattan.reduced[opt.lambda.ind.GACV25.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate25.manhattan.reduced, imputed_CV_list25.manhattan.reduced, type="l", main = "GIQCV (tau=0.25)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV25.manhattan.reduced, imputed_CV_list25.manhattan.reduced[opt.lambda.ind.imputedCV25.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV_list.manhattan.reduced, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.manhattan.reduced, GCV_list.manhattan.reduced[opt.lambda.ind.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced[opt.lambda.ind.imputedCV.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV.tps.manhattan.reduced, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.manhattan.reduced, GCV.tps.manhattan.reduced[opt.lambda.ind.tps.manhattan.reduced], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate10.manhattan.reduced, SIC_list10.manhattan.reduced, type="l", main = "SIC (tau=0.1)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC10.manhattan.reduced, SIC_list10.manhattan.reduced[opt.lambda.ind.SIC10.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate10.manhattan.reduced, GACV_list10.manhattan.reduced, type="l", main = "GACV (tau=0.1)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV10.manhattan.reduced, GACV_list10.manhattan.reduced[opt.lambda.ind.GACV10.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate10.manhattan.reduced, imputed_CV_list10.manhattan.reduced, type="l", main = "GIQCV (tau=0.1)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV10.manhattan.reduced, imputed_CV_list10.manhattan.reduced[opt.lambda.ind.imputedCV10.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV_list.manhattan.reduced, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.manhattan.reduced, GCV_list.manhattan.reduced[opt.lambda.ind.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.manhattan.reduced, imputed_CV_list.mean.manhattan.reduced[opt.lambda.ind.imputedCV.mean.manhattan.reduced], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.manhattan.reduced, GCV.tps.manhattan.reduced, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.manhattan.reduced, GCV.tps.manhattan.reduced[opt.lambda.ind.tps.manhattan.reduced], col="red", pch=19, cex=2)

fit.tps.manhattan.reduced <- Tps(M.taxi.reduced$xy, M.taxi.reduced$f, lambda = opt.lambda.tps.manhattan.reduced)

# bivariate QSS
rqss50.manhattan.reduced <- rqss(z ~ qss(cbind(x, y)), tau=0.5,
              data=data.frame(x=M.taxi.reduced$xy[,1], y=M.taxi.reduced$xy[,2], z=M.taxi.reduced$f))
rqss90.manhattan.reduced <- rqss(z ~ qss(cbind(x, y)), tau=0.9,
               data=data.frame(x=M.taxi.reduced$xy[,1], y=M.taxi.reduced$xy[,2], z=M.taxi.reduced$f))
rqss10.manhattan.reduced <- rqss(z ~ qss(cbind(x, y)), tau=0.1,
               data=data.frame(x=M.taxi.reduced$xy[,1], y=M.taxi.reduced$xy[,2], z=M.taxi.reduced$f))


# 90% - 10%  +- q(alpha)*sigma
q1 <-  plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = 
                            as.numeric(H_list_mean.manhattan.reduced[[opt.lambda.ind.imputedCV.mean.manhattan.reduced]]%*%M.taxi.reduced$f) - qnorm(0.1)*sd(M.taxi.reduced$f), 
                          min=2.4, max=8.6, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle(expression(bold(paste("mean + ", z[0.1]*hat(sigma)))))
q2 <- plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_mean.manhattan.reduced[[opt.lambda.ind.imputedCV.mean.manhattan.reduced]]%*%M.taxi.reduced$f),
                         min=2.4, max=8.6, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle("mean")
q3 <-  plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = 
                            as.numeric(H_list_mean.manhattan.reduced[[opt.lambda.ind.imputedCV.mean.manhattan.reduced]]%*%M.taxi.reduced$f) - qnorm(0.9)*sd(M.taxi.reduced$f), 
                          min=2.4, max=8.6, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle(expression(bold(paste("mean - ", z[0.1]*hat(sigma)))))
q4 <- plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt90.manhattan.reduced[[opt.lambda.ind.imputedCV90.manhattan.reduced]]%*%M.taxi.reduced$f),
                         min=2.4, max=8.6, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle("tau = 0.9")
q5 <- plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt50.manhattan.reduced[[opt.lambda.ind.imputedCV50.manhattan.reduced]]%*%M.taxi.reduced$f),
                         min=2.4, max=8.6, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle("tau = 0.5")
q6 <- plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt10.manhattan.reduced[[opt.lambda.ind.imputedCV10.manhattan.reduced]]%*%M.taxi.reduced$f),
                         min=2.4, max=8.6, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle("tau = 0.1")
grid.arrange(q1,q2,q3,q4,q5,q6, nrow=2)


## minus mean / median
q1 <-  plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = 
                            as.numeric(H_list_mean.manhattan.reduced[[opt.lambda.ind.imputedCV.mean.manhattan.reduced]]%*%M.taxi.reduced$f) - qnorm(0.1)*sd(M.taxi.reduced$f) -
                            as.numeric(H_list_mean.manhattan.reduced[[opt.lambda.ind.imputedCV.mean.manhattan.reduced]]%*%M.taxi.reduced$f), 
                          min=-1.4, max=1.4, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle(expression(bold(paste(z[0.1]))))
q2 <-  plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = 
                            as.numeric(H_list_mean.manhattan.reduced[[opt.lambda.ind.imputedCV.mean.manhattan.reduced]]%*%M.taxi.reduced$f) - qnorm(0.9)*sd(M.taxi.reduced$f) -
                            as.numeric(H_list_mean.manhattan.reduced[[opt.lambda.ind.imputedCV.mean.manhattan.reduced]]%*%M.taxi.reduced$f), 
                          min=-1.4, max=1.4, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle(expression(bold(paste(-z[0.1]))))
q3 <- plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = 
                           as.numeric(H_list_qt90.manhattan.reduced[[opt.lambda.ind.imputedCV90.manhattan.reduced]]%*%M.taxi.reduced$f) - 
                           as.numeric(H_list_qt50.manhattan.reduced[[opt.lambda.ind.imputedCV50.manhattan.reduced]]%*%M.taxi.reduced$f),
                         min=-1.4, max=1.4, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle(expression(bold(hat(s)[0.9] - hat(s)[0.5])))
q4 <- plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt10.manhattan.reduced[[opt.lambda.ind.imputedCV10.manhattan.reduced]]%*%M.taxi.reduced$f) -
                           as.numeric(H_list_qt50.manhattan.reduced[[opt.lambda.ind.imputedCV50.manhattan.reduced]]%*%M.taxi.reduced$f),
                         min=-1.4, max=1.4, value="Total Amount", label.title.size=19, label.text.size = 17, ratio=1) + ggtitle(expression(bold(hat(s)[0.1] - hat(s)[0.5])))

grid.arrange(q1,q2,q3,q4, nrow=1)


###################################
### 02. US hourly temperature data 
###################################
data_ustemp <- readMat("US_Hourly_2010_August_1st.mat") # hourly data on 2010.08.01
UShourlytemp <- list()
UShourlytemp$xy <- cbind(data_ustemp$x, data_ustemp$y)
rownames(UShourlytemp$xy) <- 1:nrow(UShourlytemp$xy) 
UShourlytemp$f <- data_ustemp$signal[,15]

distmat.htemp <- distm(UShourlytemp$xy, fun = distHaversine) / 1000
A.htemp <- c()
for(i in 1:(nrow(distmat.htemp)-1)){
  for(j in (i+1):ncol(distmat.htemp)){
    val <- distmat.htemp[i,j]
    A.htemp <- rbind(A.htemp, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.htemp, k=7), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.htemp[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=length(UShourlytemp$f), ncol=length(UShourlytemp$f))

colnames(wmat) <- 1:length(UShourlytemp$f)  
rownames(wmat) <- 1:length(UShourlytemp$f)

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}

# weight matrix
UShourlytemp$A <- wmat

# sparse weight matrix
UShourlytemp$sA <- sp.wmat

UShourlytemp$dist <- distmat.htemp
UShourlytemp$sdist <- A.htemp

# visualize
plot_graph(UShourlytemp)
plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=6, vertex_color = UShourlytemp$f, value="Temperature")

L.ush <- gasper::laplacian_mat(UShourlytemp$A)
n.ush <- nrow(UShourlytemp$xy)

## mean-based
H_list_mean.ush <- list()
GCV_list.ush <- c()
imputed_CV_list.mean.ush <- c()
lambda_candidate.mean.ush <- exp(seq(-17, -8, by=0.3)) 

for(lambda in lambda_candidate.mean.ush){
  H <- solve(diag(1,n.ush)+n.ush*lambda*L.ush%*%L.ush)
  H_list_mean.ush[[length(H_list_mean.ush)+1]] <- H
  GCV_list.ush <- c(GCV_list.ush, GCV(lambda, UShourlytemp$f, H))
  imputed_CV <- imputedCV(lambda=lambda, y=UShourlytemp$f, tau=0.5, alpha=0.1, B=diag(n.ush), 
                          P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.ush, k=1, M=NULL, method="mean")
  
  imputed_CV_list.mean.ush <- c(imputed_CV_list.mean.ush, imputed_CV)
}

plot(lambda_candidate.mean.ush, GCV_list.ush, type="l")
plot(lambda_candidate.mean.ush, imputed_CV_list.mean.ush, type="l")
opt.lambda.ind.mean.ush <- which.min(GCV_list.ush)
opt.lambda.mean.ush <- lambda_candidate.mean.ush[opt.lambda.ind.mean.ush] 
opt.lambda.ind.imputedCV.mean.ush <- which.min(imputed_CV_list.mean.ush)
opt.lambda.imputedCV.mean.ush <- lambda_candidate.mean.ush[opt.lambda.ind.imputedCV.mean.ush] 

# thin plate spline
par(mfrow=c(1,1))
GCV.tps.ush <- sapply(lambda_candidate.mean.ush, function(l){
  tmp <- Tps(UShourlytemp$xy, UShourlytemp$f, lambda=l)
  return(unlist(summary(tmp))$sum.gcv.lambda.GCV)})

plot(lambda_candidate.mean.ush, GCV.tps.ush, type="l")
opt.lambda.ind.tps.ush <- which.min(GCV.tps.ush)
opt.lambda.tps.ush <- lambda_candidate.mean.ush[opt.lambda.ind.tps.ush]


# quantile-based

# tau=0.9
SIC_list90.ush <- c()
GACV_list90.ush <- c()
imputed_CV_list90.ush <- c()
H_list_qt90.ush <- list()

lambda_candidate90.ush <- exp(seq(-8, -5, by=0.3))

for(lambda in lambda_candidate90.ush){
  res90.ush <- PIRLS(lambda=lambda, y=UShourlytemp$f, tau=0.9, alpha=0.1, B=diag(n.ush), 
                     P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list90.ush <- c(SIC_list90.ush, res90.ush$SIC)
  GACV_list90.ush <- c(GACV_list90.ush, res90.ush$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=UShourlytemp$f, tau=0.9, alpha=0.1, B=diag(n.ush), 
                             P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, 
                             check="Oh", L=L.ush, k=1, M=NULL, method="quantile")
  imputed_CV_list90.ush <- c(imputed_CV_list90.ush, imputed_CV)
  H_list_qt90.ush[[length(H_list_qt90.ush)+1]] <- res90.ush$H
}

plot(lambda_candidate90.ush, SIC_list90.ush, type="l")
plot(lambda_candidate90.ush, GACV_list90.ush, type="l")
plot(lambda_candidate90.ush, imputed_CV_list90.ush, type="l")
opt.lambda.ind.SIC90.ush <- which.min(SIC_list90.ush)
opt.lambda.SIC90.ush <- lambda_candidate90.ush[opt.lambda.ind.SIC90.ush] 
opt.lambda.ind.GACV90.ush <- which.min(GACV_list90.ush)
opt.lambda.GACV90.ush <- lambda_candidate90.ush[opt.lambda.ind.GACV90.ush] 
opt.lambda.ind.imputedCV90.ush <- which.min(imputed_CV_list90.ush)
opt.lambda.imputedCV90.ush <- lambda_candidate90.ush[opt.lambda.ind.imputedCV90.ush] 


# tau=0.75
SIC_list75.ush <- c()
GACV_list75.ush <- c()
imputed_CV_list75.ush <- c()
H_list_qt75.ush <- list()

lambda_candidate75.ush <- exp(seq(-10, -6, by=0.3))

for(lambda in lambda_candidate75.ush){
  res75.ush <- PIRLS(lambda=lambda, y=UShourlytemp$f, tau=0.75, alpha=0.1, B=diag(n.ush), 
                     P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list75.ush <- c(SIC_list75.ush, res75.ush$SIC)
  GACV_list75.ush <- c(GACV_list75.ush, res75.ush$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=UShourlytemp$f, tau=0.75, alpha=0.1, B=diag(n.ush), 
                            P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.ush, k=1, M=NULL, method="quantile")
  imputed_CV_list75.ush <- c(imputed_CV_list75.ush, imputed_CV)
  H_list_qt75.ush[[length(H_list_qt75.ush)+1]] <- res75.ush$H
}

plot(lambda_candidate75.ush, SIC_list75.ush, type="l")
plot(lambda_candidate75.ush, GACV_list75.ush, type="l")
plot(lambda_candidate75.ush, imputed_CV_list75.ush, type="l")
opt.lambda.ind.SIC75.ush <- which.min(SIC_list75.ush)
opt.lambda.SIC75.ush <- lambda_candidate75.ush[opt.lambda.ind.SIC75.ush]
opt.lambda.ind.GACV75.ush <- which.min(GACV_list75.ush)
opt.lambda.GACV75.ush <- lambda_candidate75.ush[opt.lambda.ind.GACV75.ush]
opt.lambda.ind.imputedCV75.ush <- which.min(imputed_CV_list75.ush)
opt.lambda.imputedCV75.ush <- lambda_candidate75.ush[opt.lambda.ind.imputedCV75.ush] 

# tau=0.5
SIC_list50.ush <- c()
GACV_list50.ush <- c()
imputed_CV_list50.ush <- c()
H_list_qt50.ush <- list()

lambda_candidate50.ush <- exp(seq(-11, -7, by=0.3))

for(lambda in lambda_candidate50.ush){
  res50.ush <- PIRLS(lambda=lambda, y=UShourlytemp$f, tau=0.5, alpha=0.1, B=diag(n.ush), 
                     P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list50.ush <- c(SIC_list50.ush, res50.ush$SIC)
  GACV_list50.ush <- c(GACV_list50.ush, res50.ush$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=UShourlytemp$f, tau=0.5, alpha=0.1, B=diag(n.ush), 
                            P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.ush, k=1, M=NULL, method="quantile")
  imputed_CV_list50.ush <- c(imputed_CV_list50.ush, imputed_CV)
  H_list_qt50.ush[[length(H_list_qt50.ush)+1]] <- res50.ush$H
}

plot(lambda_candidate50.ush, SIC_list50.ush, type="l")
plot(lambda_candidate50.ush, GACV_list50.ush, type="l")
plot(lambda_candidate50.ush, imputed_CV_list50.ush, type="l")
opt.lambda.ind.SIC50.ush <- which.min(SIC_list50.ush)
opt.lambda.SIC50.ush <- lambda_candidate50.ush[opt.lambda.ind.SIC50.ush]
opt.lambda.ind.GACV50.ush <- which.min(GACV_list50.ush)
opt.lambda.GACV50.ush <- lambda_candidate50.ush[opt.lambda.ind.GACV50.ush]
opt.lambda.ind.imputedCV50.ush <- which.min(imputed_CV_list50.ush)
opt.lambda.imputedCV50.ush <- lambda_candidate50.ush[opt.lambda.ind.imputedCV50.ush] 

# tau=0.25
SIC_list25.ush <- c()
GACV_list25.ush <- c()
imputed_CV_list25.ush <- c()
H_list_qt25.ush <- list()

lambda_candidate25.ush <- exp(seq(-8, -5, by=0.3))

for(lambda in lambda_candidate25.ush){
  res25.ush <- PIRLS(lambda=lambda, y=UShourlytemp$f, tau=0.25, alpha=0.1, B=diag(n.ush), 
                     P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list25.ush <- c(SIC_list25.ush, res25.ush$SIC)
  GACV_list25.ush <- c(GACV_list25.ush, res25.ush$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=UShourlytemp$f, tau=0.25, alpha=0.1, B=diag(n.ush), 
                            P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.ush, k=1, M=NULL, method="quantile")
  imputed_CV_list25.ush <- c(imputed_CV_list25.ush, imputed_CV)
  H_list_qt25.ush[[length(H_list_qt25.ush)+1]] <- res25.ush$H
}

plot(lambda_candidate25.ush, SIC_list25.ush, type="l")
plot(lambda_candidate25.ush, GACV_list25.ush, type="l")
plot(lambda_candidate25.ush, imputed_CV_list25.ush, type="l")
opt.lambda.ind.SIC25.ush <- which.min(SIC_list25.ush)
opt.lambda.SIC25.ush <- lambda_candidate25.ush[opt.lambda.ind.SIC25.ush]
opt.lambda.ind.GACV25.ush <- which.min(GACV_list25.ush)
opt.lambda.GACV25.ush <- lambda_candidate25.ush[opt.lambda.ind.GACV25.ush]
opt.lambda.ind.imputedCV25.ush <- which.min(imputed_CV_list25.ush)
opt.lambda.imputedCV25.ush <- lambda_candidate25.ush[opt.lambda.ind.imputedCV25.ush] 


# tau=0.1
SIC_list10.ush <- c()
GACV_list10.ush <- c()
imputed_CV_list10.ush <- c()
H_list_qt10.ush <- list()

lambda_candidate10.ush <- exp(seq(-11, -7, by=0.3))

for(lambda in lambda_candidate10.ush){
  res10.ush <- PIRLS(lambda=lambda, y=UShourlytemp$f, tau=0.1, alpha=0.1, B=diag(n.ush), 
                     P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list10.ush <- c(SIC_list10.ush, res10.ush$SIC)
  GACV_list10.ush <- c(GACV_list10.ush, res10.ush$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=UShourlytemp$f, tau=0.1, alpha=0.1, B=diag(n.ush), 
                            P=L.ush%*%L.ush, gamma_init=rep(mean(UShourlytemp$f), n.ush), max_iterations=10000, tol=1e-6, 
                            check="Oh", L=L.ush, k=1, M=NULL, method="quantile")
  imputed_CV_list10.ush <- c(imputed_CV_list10.ush, imputed_CV)
  H_list_qt10.ush[[length(H_list_qt10.ush)+1]] <- res10.ush$H
}

plot(lambda_candidate10.ush, SIC_list10.ush, type="l")
plot(lambda_candidate10.ush, GACV_list10.ush, type="l")
plot(lambda_candidate10.ush, imputed_CV_list10.ush, type="l")
opt.lambda.ind.SIC10.ush <- which.min(SIC_list10.ush)
opt.lambda.SIC10.ush <- lambda_candidate10.ush[[opt.lambda.ind.SIC10.ush]] # 0.001230912
opt.lambda.ind.GACV10.ush <- which.min(GACV_list10.ush)
opt.lambda.GACV10.ush <- lambda_candidate10.ush[[opt.lambda.ind.GACV10.ush]] # 0.000911882
opt.lambda.ind.imputedCV10.ush <- which.min(imputed_CV_list10.ush)
opt.lambda.imputedCV10.ush <- lambda_candidate10.ush[opt.lambda.ind.imputedCV10.ush] 



# lambda selection result visualization
par(mfrow=c(2,3))
plot(lambda_candidate90.ush, SIC_list90.ush, type="l", main = "SIC (tau=0.9)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC90.ush, SIC_list90.ush[opt.lambda.ind.SIC90.ush], col="red", pch=19, cex=2)
plot(lambda_candidate90.ush, GACV_list90.ush, type="l", main = "GACV (tau=0.9)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV90.ush, GACV_list90.ush[opt.lambda.ind.GACV90.ush], col="red", pch=19, cex=2)
plot(lambda_candidate90.ush, imputed_CV_list90.ush, type="l", main = "GIQCV (tau=0.9)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV90.ush, imputed_CV_list90.ush[opt.lambda.ind.imputedCV90.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV_list.ush, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.ush, GCV_list.ush[opt.lambda.ind.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, imputed_CV_list.mean.ush, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.ush, imputed_CV_list.mean.ush[opt.lambda.ind.imputedCV.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV.tps.ush, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.ush, GCV.tps.ush[opt.lambda.ind.tps.ush], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate75.ush, SIC_list75.ush, type="l", main = "SIC (tau=0.75)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC75.ush, SIC_list75.ush[opt.lambda.ind.SIC75.ush], col="red", pch=19, cex=2)
plot(lambda_candidate75.ush, GACV_list75.ush, type="l", main = "GACV (tau=0.75)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV75.ush, GACV_list75.ush[opt.lambda.ind.GACV75.ush], col="red", pch=19, cex=2)
plot(lambda_candidate75.ush, imputed_CV_list75.ush, type="l", main = "GIQCV (tau=0.75)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV75.ush, imputed_CV_list75.ush[opt.lambda.ind.imputedCV75.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV_list.ush, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.ush, GCV_list.ush[opt.lambda.ind.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, imputed_CV_list.mean.ush, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.ush, imputed_CV_list.mean.ush[opt.lambda.ind.imputedCV.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV.tps.ush, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.ush, GCV.tps.ush[opt.lambda.ind.tps.ush], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate50.ush, SIC_list50.ush, type="l", main = "SIC (tau=0.5)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC50.ush, SIC_list50.ush[opt.lambda.ind.SIC50.ush], col="red", pch=19, cex=2)
plot(lambda_candidate50.ush, GACV_list50.ush, type="l", main = "GACV (tau=0.5)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV50.ush, GACV_list50.ush[opt.lambda.ind.GACV50.ush], col="red", pch=19, cex=2)
plot(lambda_candidate50.ush, imputed_CV_list50.ush, type="l", main = "GIQCV (tau=0.5)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV50.ush, imputed_CV_list50.ush[opt.lambda.ind.imputedCV50.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV_list.ush, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.ush, GCV_list.ush[opt.lambda.ind.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, imputed_CV_list.mean.ush, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.ush, imputed_CV_list.mean.ush[opt.lambda.ind.imputedCV.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV.tps.ush, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.ush, GCV.tps.ush[opt.lambda.ind.tps.ush], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate25.ush, SIC_list25.ush, type="l", main = "SIC (tau=0.25)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC25.ush, SIC_list25.ush[opt.lambda.ind.SIC25.ush], col="red", pch=19, cex=2)
plot(lambda_candidate25.ush, GACV_list25.ush, type="l", main = "GACV (tau=0.25)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV25.ush, GACV_list25.ush[opt.lambda.ind.GACV25.ush], col="red", pch=19, cex=2)
plot(lambda_candidate25.ush, imputed_CV_list25.ush, type="l", main = "GIQCV (tau=0.25)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV25.ush, imputed_CV_list25.ush[opt.lambda.ind.imputedCV25.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV_list.ush, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.ush, GCV_list.ush[opt.lambda.ind.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, imputed_CV_list.mean.ush, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.ush, imputed_CV_list.mean.ush[opt.lambda.ind.imputedCV.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV.tps.ush, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.ush, GCV.tps.ush[opt.lambda.ind.tps.ush], col="red", pch=19, cex=2)

par(mfrow=c(2,3))
plot(lambda_candidate10.ush, SIC_list10.ush, type="l", main = "SIC (tau=0.1)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.SIC10.ush, SIC_list10.ush[opt.lambda.ind.SIC10.ush], col="red", pch=19, cex=2)
plot(lambda_candidate10.ush, GACV_list10.ush, type="l", main = "GACV (tau=0.1)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.GACV10.ush, GACV_list10.ush[opt.lambda.ind.GACV10.ush], col="red", pch=19, cex=2)
plot(lambda_candidate10.ush, imputed_CV_list10.ush, type="l", main = "GIQCV (tau=0.1)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV10.ush, imputed_CV_list10.ush[opt.lambda.ind.imputedCV10.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV_list.ush, type="l", main = "GCV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.mean.ush, GCV_list.ush[opt.lambda.ind.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, imputed_CV_list.mean.ush, type="l", main = "GICV (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.imputedCV.mean.ush, imputed_CV_list.mean.ush[opt.lambda.ind.imputedCV.mean.ush], col="red", pch=19, cex=2)
plot(lambda_candidate.mean.ush, GCV.tps.ush, type="l", main = "TPS (mean)", cex.main=2,
     ylab="score", xlab="lambda", cex.axis=1.2, cex.lab=1.5)
points(opt.lambda.tps.ush, GCV.tps.ush[opt.lambda.ind.tps.ush], col="red", pch=19, cex=2)

# bivariate QSS
rqss50.ush <- rqss(z ~ qss(cbind(x, y)), tau=0.5,
                                 data=data.frame(x=UShourlytemp$xy[,1], y=UShourlytemp$xy[,2], z=UShourlytemp$f))
rqss90.ush <- rqss(z ~ qss(cbind(x, y)), tau=0.9,
                                 data=data.frame(x=UShourlytemp$xy[,1], y=UShourlytemp$xy[,2], z=UShourlytemp$f))
rqss10.ush <- rqss(z ~ qss(cbind(x, y)), tau=0.1,
                                 data=data.frame(x=UShourlytemp$xy[,1], y=UShourlytemp$xy[,2], z=UShourlytemp$f))

# 90% - 10%  +- q(alpha)*sigma
q1 <-  plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = 
                            as.numeric(H_list_mean.ush[[opt.lambda.ind.imputedCV.mean.ush]]%*%UShourlytemp$f) - qnorm(0.1)*sd(UShourlytemp$f), 
                          min=59, max=109, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle(expression(bold(paste("mean + ", z[0.1]*hat(sigma)))))
q2 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_mean.ush[[opt.lambda.ind.imputedCV.mean.ush]]%*%UShourlytemp$f),
                         min=59, max=109, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle("mean")
q3 <-  plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = 
                            as.numeric(H_list_mean.ush[[opt.lambda.ind.imputedCV.mean.ush]]%*%UShourlytemp$f) - qnorm(0.9)*sd(UShourlytemp$f), 
                          min=59, max=109, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle(expression(bold(paste("mean - ", z[0.1]*hat(sigma)))))
q4 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt90.ush[[opt.lambda.ind.imputedCV90.ush]]%*%UShourlytemp$f),
                         min=59, max=109, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle("tau = 0.9")
q5 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt50.ush[[opt.lambda.ind.imputedCV50.ush]]%*%UShourlytemp$f),
                         min=59, max=109, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle("tau = 0.5")
q6 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt10.ush[[opt.lambda.ind.imputedCV10.ush]]%*%UShourlytemp$f),
                         min=59, max=109, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle("tau = 0.1")

grid.arrange(q1,q2,q3,q4,q5,q6, nrow=2)


## minus mean / median
q1 <-  plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = 
                            as.numeric(H_list_mean.ush[[opt.lambda.ind.imputedCV.mean.ush]]%*%UShourlytemp$f) - qnorm(0.1)*sd(UShourlytemp$f) -
                            as.numeric(H_list_mean.ush[[opt.lambda.ind.imputedCV.mean.ush]]%*%UShourlytemp$f), 
                          min=-21.5, max=22.5, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle(expression(bold(paste(z[0.1]))))
q2 <-  plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color =
                            as.numeric(H_list_mean.ush[[opt.lambda.ind.imputedCV.mean.ush]]%*%UShourlytemp$f) - qnorm(0.9)*sd(UShourlytemp$f) -
                            as.numeric(H_list_mean.ush[[opt.lambda.ind.imputedCV.mean.ush]]%*%UShourlytemp$f),
                          min=-21.5, max=22.5, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle(expression(bold(paste(-z[0.1]))))
q3 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = 
                           as.numeric(H_list_qt90.ush[[opt.lambda.ind.imputedCV90.ush]]%*%UShourlytemp$f) - 
                           as.numeric(H_list_qt50.ush[[opt.lambda.ind.imputedCV50.ush]]%*%UShourlytemp$f),
                         min=-21.5, max=22.5, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle(expression(bold(hat(s)[0.9] - hat(s)[0.5])))
q4 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = as.numeric(H_list_qt10.ush[[opt.lambda.ind.imputedCV10.ush]]%*%UShourlytemp$f) -
                           as.numeric(H_list_qt50.ush[[opt.lambda.ind.imputedCV50.ush]]%*%UShourlytemp$f),
                         min=-21.5, max=22.5, value="Temperature", label.title.size=19, label.text.size = 17, ratio=1/1.3) + ggtitle(expression(bold(hat(s)[0.1] - hat(s)[0.5])))

grid.arrange(q1,q2,q3,q4, nrow=1)


# visualize graph and original graph signal
# graph
q1 <- plot_graph_custom3(M.taxi.reduced, v.size=1.5, ratio=1, signal=FALSE)
q2 <- plot_graph_custom3(UShourlytemp, v.size=1.5, ratio=1/1.3, signal=FALSE)
grid.arrange(q1,q2, ncol=2)

# signal
q3 <- plot_graph_custom3(M.taxi.reduced, e.size=1.3, v.size=8, vertex_color = M.taxi.reduced$f,
                         min=2.4, max=8.6, value="Total Amount", label.title.size=17, label.text.size = 15, ratio=1)
q4 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=8, vertex_color = UShourlytemp$f,
                         min=59, max=109, value="Temperature", label.title.size=17, label.text.size = 15, ratio=1/1.3)

grid.arrange(q3,q4, ncol=2)
