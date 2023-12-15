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
library(quantreg)

source("utils.R")

######################################################################
### Motivation 01: Graph-based vs Euclidean-based (spiral lattice)
######################################################################

# Define the dimensions of the grid
n <- 15

# Initialize an empty adjacency matrix
adj_matrix <- matrix(0, ncol = n * n, nrow = n * n)

# Define clockwise and counterclockwise movement (you mentioned counterclockwise)
# These vectors represent the changes in x and y coordinates for each step
dx <- c(1, 0, -1, 0)
dy <- c(0, 1, 0, -1)

# Initialize starting position and direction
x <- 0
y <- 0
direction <- 1  # 0 = right, 1 = up, 2 = left, 3 = down

# Helper function to convert (x, y) coordinates to node index
coord_to_index <- function(x, y) {
  return(y * n + x + 1)
}

x_start <- 0
y_start <- 1
x_end <- n-1
y_end <- n-1
node.order <- c()
# Iterate through all nodes in a counterclockwise order
for (i in 1:(n * n-1)) {
  node.order <- c(node.order,coord_to_index(x, y))
  if(i==(n * n-1)){
    node.order <- c(node.order, coord_to_index(x + dx[direction], y + dy[direction]))
  }
  # cat("x", x, "y", y)
  # Mark the current node as connected to the next node in the current direction
  adj_matrix[coord_to_index(x, y), coord_to_index(x + dx[direction], y + dy[direction])] <- 1
  
  # Calculate the next position
  x_next <- x + dx[direction]
  y_next <- y + dy[direction]
  
  # Check if the next position is within bounds and not yet visited
  if(direction==1){
    if(x_next < x_end){
      x <- x_next
      y <- y_next
    } else {
      direction <- (direction + 1) %% 4
      x <- x_next
      y <- y_next
      x_end <- x_end - 1
    }
  } else if(direction==2){
    if(y_next < y_end){
      x <- x_next
      y <- y_next
    } else {
      direction <- (direction + 1) %% 4
      x <- x_next
      y <- y_next
      y_end <- y_end - 1
    }
  } else if(direction==3){
    if(x_next > x_start){
      x <- x_next
      y <- y_next
    } else {
      direction <- 4
      x <- x_next
      y <- y_next
      x_start <- x_start + 1
    } 
  } else if(direction==4){
    if(y_next > y_start){
      x <- x_next
      y <- y_next
    } else {
      direction <- (direction + 1) %% 4
      x <- x_next
      y <- y_next
      y_start <- y_start + 1
    }
  }
}

adj_matrix <- adj_matrix + t(adj_matrix)

g2 = make_lattice(c(n,n)) #make 30 x 30 grid graph
g2 <- delete_edges(g2, E(g2))  # Remove all existing edges
x = rep(0:(n-1), n) #define coordinates of vertices
y = rep(0:(n-1), rep(n,n))
graph_attr(g2, 'xy') = data.frame(x,y) #define 'xy' attributes of graph g which represents coordinates
g2.ad_mat = Matrix(adj_matrix, sparse=TRUE) #get unweighted adjacency matrix of graph g
graph_attr(g2, 'sA') = g2.ad_mat
for (i in 1:(nrow(adj_matrix)-1)) {
  for (j in (i+1):ncol(adj_matrix)) {
    if (adj_matrix[i, j] == 1) {
      # Add an edge between nodes i and j
      g2 <- add_edges(g2, c(i, j))
    }
  }
}
lattice02 <- list()
lattice02$xy <- g2$xy
n.lattice02 <- nrow(lattice02$xy)
rownames(lattice02$xy) <- 1:n.lattice02
N <- 100

# weight matrix
lattice02$A <- g2$sA
edge.wt <- igraph::as_data_frame(g2, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              g2$sA[edge.wt[i,1], edge.wt[i,2]]))
}

# sparse weight matrix
lattice02$sA <- sp.wmat
L.lattice02 <- gasper::laplacian_mat(lattice02$A)

# signal assignment
# lattice02$f.true[node.order] <- sin(2*pi*c(0:(n*n-1))/7) # true sine signal
lattice02$f.true[node.order] <- c(1:225)/20 # true line signal

plot_graph(lattice02) # plot graph structure
plot_graph_custom3(lattice02, e.size=1.3, v.size=6, vertex_color = lattice02$f.true, 
                   min=min(lattice02$f.true), max=max(lattice02$f.true), value="value") + ggtitle("original") # plot graph signal

set.seed(100)

lattice02$f <- lattice02$f.true + rnorm(n*n,0,0.5)
plot(lattice02$f[node.order], type="l")

## mean-based
H_list_mean.lattice02 <- list()
GCV_list.lattice02 <- c()
imputed_CV_list.mean.lattice02 <- c()
lambda_candidate.mean.lattice02 <- exp(seq(2, 5, by=0.3))
for(lambda in lambda_candidate.mean.lattice02){
  H <- solve(diag(1,n.lattice02)+n.lattice02*lambda*L.lattice02%*%L.lattice02)
  H_list_mean.lattice02[[length(H_list_mean.lattice02)+1]] <- H
  GCV_list.lattice02 <- c(GCV_list.lattice02, GCV(lambda, lattice02$f, H))
  imputed_CV <- imputedCV(lambda=lambda, y=lattice02$f, tau=0.5, alpha=0.1, B=diag(n.lattice02), 
                          P=L.lattice02%*%L.lattice02, gamma_init=rep(mean(lattice02$f), n.lattice02), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice02, k=1, M=NULL, method="mean")
  imputed_CV_list.mean.lattice02 <- c(imputed_CV_list.mean.lattice02, imputed_CV)
}
plot(lambda_candidate.mean.lattice02, GCV_list.lattice02, type="l")
plot(lambda_candidate.mean.lattice02, imputed_CV_list.mean.lattice02, type="l")
opt.lambda.ind.mean.lattice02 <- which.min(GCV_list.lattice02)
opt.lambda.mean.lattice02 <- lambda_candidate.mean.lattice02[opt.lambda.ind.mean.lattice02] 
opt.lambda.ind.imputedCV.mean.lattice02 <- which.min(imputed_CV_list.mean.lattice02)
opt.lambda.imputedCV.mean.lattice02 <- lambda_candidate.mean.lattice02[opt.lambda.ind.imputedCV.mean.lattice02] 

# thin plate spline
lambda_candidate.mean.lattice02 <- exp(seq(-17, -4, by=0.3))
GCV.tps.lattice02 <- sapply(lambda_candidate.mean.lattice02, function(l){
  tmp <- Tps(lattice02$xy, lattice02$f, lambda=l)
  return(unlist(summary(tmp))$sum.gcv.lambda.GCV)})
plot(lambda_candidate.mean.lattice02, GCV.tps.lattice02, type="l")
opt.lambda.ind.tps.lattice02 <- which.min(GCV.tps.lattice02)
opt.lambda.tps.lattice02 <- lambda_candidate.mean.lattice02[opt.lambda.ind.tps.lattice02]
fit.tps.lattice02 <- Tps(lattice02$xy, lattice02$f, lambda = opt.lambda.tps.lattice02)

## quantile-based
# bivariate QSS
rqss50.lattice02 <- rqss(z ~ qss(cbind(x, y)), tau=0.5,
                         data=data.frame(x=lattice02$xy[,1], y=lattice02$xy[,2], z=lattice02$f))
rqss90.lattice02 <- rqss(z ~ qss(cbind(x, y)), tau=0.9,
                         data=data.frame(x=lattice02$xy[,1], y=lattice02$xy[,2], z=lattice02$f))
rqss10.lattice02 <- rqss(z ~ qss(cbind(x, y)), tau=0.1,
                         data=data.frame(x=lattice02$xy[,1], y=lattice02$xy[,2], z=lattice02$f))
rqss75.lattice02 <- rqss(z ~ qss(cbind(x, y)), tau=0.75,
                         data=data.frame(x=lattice02$xy[,1], y=lattice02$xy[,2], z=lattice02$f))
rqss25.lattice02 <- rqss(z ~ qss(cbind(x, y)), tau=0.25,
                         data=data.frame(x=lattice02$xy[,1], y=lattice02$xy[,2], z=lattice02$f))

# tau=0.5
SIC_list50.lattice02 <- c()
GACV_list50.lattice02 <- c()
imputed_CV_list50.lattice02 <- c()
H_list_qt50.lattice02 <- list()

lambda_candidate50.lattice02 <- exp(seq(0.5,3.5, by=0.3))
for(lambda in lambda_candidate50.lattice02){
  cat("####### lambda = ", lambda, "####### \n")
  res50.lattice02 <- PIRLS(lambda=lambda, y=lattice02$f, tau=0.5, alpha=0.1, B=diag(n.lattice02), 
                           P=L.lattice02%*%L.lattice02, gamma_init=rep(mean(lattice02$f), n.lattice02), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list50.lattice02 <- c(SIC_list50.lattice02, res50.lattice02$SIC)
  GACV_list50.lattice02 <- c(GACV_list50.lattice02, res50.lattice02$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice02$f, tau=0.5, alpha=0.1, B=diag(n.lattice02), 
                          P=L.lattice02%*%L.lattice02, gamma_init=rep(mean(lattice02$f), n.lattice02), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice02, k=1, M=NULL, method="quantile")
  imputed_CV_list50.lattice02 <- c(imputed_CV_list50.lattice02, imputed_CV)
  H_list_qt50.lattice02[[length(H_list_qt50.lattice02)+1]] <- res50.lattice02$H
}
plot(lambda_candidate50.lattice02, SIC_list50.lattice02, type="l")
plot(lambda_candidate50.lattice02, GACV_list50.lattice02, type="l")
plot(lambda_candidate50.lattice02, imputed_CV_list50.lattice02, type="l")
opt.lambda.ind.SIC50.lattice02 <- which.min(SIC_list50.lattice02)
opt.lambda.SIC50.lattice02 <- lambda_candidate50.lattice02[opt.lambda.ind.SIC50.lattice02] 
opt.lambda.ind.GACV50.lattice02 <- which.min(GACV_list50.lattice02)
opt.lambda.GACV50.lattice02 <- lambda_candidate50.lattice02[opt.lambda.ind.GACV50.lattice02] 
opt.lambda.ind.imputedCV50.lattice02 <- which.min(imputed_CV_list50.lattice02)
opt.lambda.imputedCV50.lattice02 <- lambda_candidate50.lattice02[opt.lambda.ind.imputedCV50.lattice02] 

plot(lattice02$f[node.order], type="l")
lines(fitted(rqss50.lattice02)[node.order], col="red")
lines(as.numeric(H_list_mean.lattice02[[opt.lambda.ind.imputedCV.mean.lattice02]]%*%lattice02$f)[node.order], col="blue")

# fitting results
q1 <- plot_graph_custom3(lattice02, e.size=1.3, v.size=11, vertex_color = lattice02$f.true, 
                         min=-1, max=12.1, value="value", label.title.size=19, label.text.size = 17) + ggtitle("true signal") # plot graph signal
q2 <- plot_graph_custom3(lattice02, e.size=1.3, v.size=11, vertex_color = lattice02$f, 
                         min=-1, max=12.1, value="value", label.title.size=19, label.text.size = 17) + ggtitle("noisy signal") # plot graph signal
q3 <- plot_graph_custom3(lattice02, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_mean.lattice02[[opt.lambda.ind.imputedCV.mean.lattice02]]%*%lattice02$f), 
                         min=-1, max=12.1, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("mean")
q4 <- plot_graph_custom3(lattice02, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt50.lattice02[[opt.lambda.ind.imputedCV50.lattice02]]%*%lattice02$f), 
                         min=-1, max=12.1, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.5")
q5 <- plot_graph_custom3(lattice02, e.size=1.3, v.size=11, vertex_color = fitted(fit.tps.lattice02), 
                         min=-1, max=12.1, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("TPS")
q6 <- plot_graph_custom3(lattice02, e.size=1.3, v.size=11, vertex_color = fitted(rqss50.lattice02), 
                         min=-1, max=12.1, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("QSS (tau = 0.5)")

grid.arrange(q1,q3,q5,q2,q4,q6, ncol=3)


######################################################################
### Motivation 02: Existence of outliers (spiral lattice)
######################################################################

# different signal (outlier)
lattice03 <- lattice02
lattice03$f <- rep(1, n*n)
lattice03$f[c(18:24)] <- 10
n.lattice03 <- n.lattice02
L.lattice03 <- L.lattice02
q1 <- plot_graph_custom3(lattice03, e.size=1.3, v.size=11, vertex_color = lattice03$f, 
                         min=min(lattice03$f), max=max(lattice03$f), value="value") + ggtitle("noisy signal") # plot graph signal

## mean-based
H_list_mean.lattice03 <- list()
GCV_list.lattice03 <- c()
imputed_CV_list.mean.lattice03 <- c()
lambda_candidate.mean.lattice03 <- exp(seq(-5,10, by=1))
for(lambda in lambda_candidate.mean.lattice03){
  H <- solve(diag(1,n.lattice03)+n.lattice03*lambda*L.lattice03%*%L.lattice03)
  H_list_mean.lattice03[[length(H_list_mean.lattice03)+1]] <- H
  GCV_list.lattice03 <- c(GCV_list.lattice03, GCV(lambda, lattice03$f, H))
  imputed_CV <- imputedCV(lambda=lambda, y=lattice03$f, tau=0.5, alpha=0.1, B=diag(n.lattice03), 
                          P=L.lattice03%*%L.lattice03, gamma_init=rep(mean(lattice03$f), n.lattice03), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice03, k=1, M=NULL, method="mean")
  imputed_CV_list.mean.lattice03 <- c(imputed_CV_list.mean.lattice03, imputed_CV)
}
plot(lambda_candidate.mean.lattice03, GCV_list.lattice03, type="l")
plot(lambda_candidate.mean.lattice03, imputed_CV_list.mean.lattice03, type="l")
opt.lambda.ind.mean.lattice03 <- which.min(GCV_list.lattice03)
opt.lambda.mean.lattice03 <- lambda_candidate.mean.lattice03[opt.lambda.ind.mean.lattice03] 
opt.lambda.ind.imputedCV.mean.lattice03 <- which.min(imputed_CV_list.mean.lattice03)
opt.lambda.imputedCV.mean.lattice03 <- lambda_candidate.mean.lattice03[opt.lambda.ind.imputedCV.mean.lattice03] 


# tau=0.5
SIC_list50.lattice03 <- c()
GACV_list50.lattice03 <- c()
imputed_CV_list50.lattice03 <- c()
H_list_qt50.lattice03 <- list()

lambda_candidate50.lattice03 <- exp(seq(0,5, by=1))
for(lambda in lambda_candidate50.lattice03){
  cat("####### lambda = ", lambda, "####### \n")
  res50.lattice03 <- PIRLS(lambda=lambda, y=lattice03$f, tau=0.5, alpha=0.1, B=diag(n.lattice03), 
                           P=L.lattice03%*%L.lattice03, gamma_init=rep(mean(lattice03$f), n.lattice03), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list50.lattice03 <- c(SIC_list50.lattice03, res50.lattice03$SIC)
  GACV_list50.lattice03 <- c(GACV_list50.lattice03, res50.lattice03$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice03$f, tau=0.5, alpha=0.1, B=diag(n.lattice03), 
                          P=L.lattice03%*%L.lattice03, gamma_init=rep(mean(lattice03$f), n.lattice03), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice03, k=1, M=NULL, method="quantile")
  imputed_CV_list50.lattice03 <- c(imputed_CV_list50.lattice03, imputed_CV)
  H_list_qt50.lattice03[[length(H_list_qt50.lattice03)+1]] <- res50.lattice03$H
}
plot(lambda_candidate50.lattice03, SIC_list50.lattice03, type="l")
plot(lambda_candidate50.lattice03, GACV_list50.lattice03, type="l")
plot(lambda_candidate50.lattice03, imputed_CV_list50.lattice03, type="l")
opt.lambda.ind.SIC50.lattice03 <- which.min(SIC_list50.lattice03)
opt.lambda.SIC50.lattice03 <- lambda_candidate50.lattice03[opt.lambda.ind.SIC50.lattice03] 
opt.lambda.ind.GACV50.lattice03 <- which.min(GACV_list50.lattice03)
opt.lambda.GACV50.lattice03 <- lambda_candidate50.lattice03[opt.lambda.ind.GACV50.lattice03] 
opt.lambda.ind.imputedCV50.lattice03 <- which.min(imputed_CV_list50.lattice03)
opt.lambda.imputedCV50.lattice03 <- lambda_candidate50.lattice03[opt.lambda.ind.imputedCV50.lattice03] 

# fitting results
q1 <- plot_graph_custom3(lattice03, e.size=1.3, v.size=11, vertex_color = lattice03$f, 
                         min=0, max=11, value="value", label.title.size=19, label.text.size = 17) + ggtitle("noisy signal") # plot graph signal
q2 <- plot_graph_custom3(lattice03, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_mean.lattice03[[opt.lambda.ind.imputedCV.mean.lattice03]]%*%lattice03$f), 
                         min=0, max=11, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("mean") # plot graph signal
q3 <- plot_graph_custom3(lattice03, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt50.lattice03[[opt.lambda.ind.imputedCV50.lattice03]]%*%lattice03$f), 
                         min=0, max=11, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.5") # plot graph signal

grid.arrange(q1,q2,q3, ncol=3)


######################################################################
### Motivation 03: Varying structure (spiral lattice)
######################################################################

# different signal (varying structure)
lattice04 <- lattice02
lattice04$f.true <- rep(5, n*n)
lattice04$f <- rep(5, n*n)
set.seed(10)
lattice04$f[node.order[1:110]] <- lattice04$f[node.order[1:110]] + rnorm(110,0,0.5)
lattice04$f[node.order[111:225]] <- lattice04$f[node.order[111:225]] + rnorm(115,0,2)
n.lattice04 <- n.lattice02
L.lattice04 <- L.lattice02
q1 <- plot_graph_custom3(lattice04, e.size=1.3, v.size=11, vertex_color = lattice04$f, 
                         min=min(lattice04$f), max=max(lattice04$f), value="value") + ggtitle("noisy signal") # plot graph signal


## mean-based
H_list_mean.lattice04 <- list()
GCV_list.lattice04 <- c()
imputed_CV_list.mean.lattice04 <- c()
lambda_candidate.mean.lattice04 <- exp(seq(-10,10, by=1))
for(lambda in lambda_candidate.mean.lattice04){
  H <- solve(diag(1,n.lattice04)+n.lattice04*lambda*L.lattice04%*%L.lattice04)
  H_list_mean.lattice04[[length(H_list_mean.lattice04)+1]] <- H
  GCV_list.lattice04 <- c(GCV_list.lattice04, GCV(lambda, lattice04$f, H))
  imputed_CV <- imputedCV(lambda=lambda, y=lattice04$f, tau=0.5, alpha=0.1, B=diag(n.lattice04), 
                          P=L.lattice04%*%L.lattice04, gamma_init=rep(mean(lattice04$f), n.lattice04), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice04, k=1, M=NULL, method="mean")
  imputed_CV_list.mean.lattice04 <- c(imputed_CV_list.mean.lattice04, imputed_CV)
}
plot(lambda_candidate.mean.lattice04, GCV_list.lattice04, type="l")
plot(lambda_candidate.mean.lattice04, imputed_CV_list.mean.lattice04, type="l")
opt.lambda.ind.mean.lattice04 <- which.min(GCV_list.lattice04)
opt.lambda.mean.lattice04 <- lambda_candidate.mean.lattice04[opt.lambda.ind.mean.lattice04] 
opt.lambda.ind.imputedCV.mean.lattice04 <- which.min(imputed_CV_list.mean.lattice04)
opt.lambda.imputedCV.mean.lattice04 <- lambda_candidate.mean.lattice04[opt.lambda.ind.imputedCV.mean.lattice04] 

# tau=0.5
SIC_list50.lattice04 <- c()
GACV_list50.lattice04 <- c()
imputed_CV_list50.lattice04 <- c()
H_list_qt50.lattice04 <- list()

lambda_candidate50.lattice04 <- exp(seq(0,5, by=1))
for(lambda in lambda_candidate50.lattice04){
  cat("####### lambda = ", lambda, "####### \n")
  res50.lattice04 <- PIRLS(lambda=lambda, y=lattice04$f, tau=0.5, alpha=0.1, B=diag(n.lattice04), 
                           P=L.lattice04%*%L.lattice04, gamma_init=rep(mean(lattice04$f), n.lattice04), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list50.lattice04 <- c(SIC_list50.lattice04, res50.lattice04$SIC)
  GACV_list50.lattice04 <- c(GACV_list50.lattice04, res50.lattice04$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice04$f, tau=0.5, alpha=0.1, B=diag(n.lattice04), 
                          P=L.lattice04%*%L.lattice04, gamma_init=rep(mean(lattice04$f), n.lattice04), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice04, k=1, M=NULL, method="quantile")
  imputed_CV_list50.lattice04 <- c(imputed_CV_list50.lattice04, imputed_CV)
  H_list_qt50.lattice04[[length(H_list_qt50.lattice04)+1]] <- res50.lattice04$H
}
plot(lambda_candidate50.lattice04, SIC_list50.lattice04, type="l")
plot(lambda_candidate50.lattice04, GACV_list50.lattice04, type="l")
plot(lambda_candidate50.lattice04, imputed_CV_list50.lattice04, type="l")
opt.lambda.ind.SIC50.lattice04 <- which.min(SIC_list50.lattice04)
opt.lambda.SIC50.lattice04 <- lambda_candidate50.lattice04[opt.lambda.ind.SIC50.lattice04] 
opt.lambda.ind.GACV50.lattice04 <- which.min(GACV_list50.lattice04)
opt.lambda.GACV50.lattice04 <- lambda_candidate50.lattice04[opt.lambda.ind.GACV50.lattice04] 
opt.lambda.ind.imputedCV50.lattice04 <- which.min(imputed_CV_list50.lattice04)
opt.lambda.imputedCV50.lattice04 <- lambda_candidate50.lattice04[opt.lambda.ind.imputedCV50.lattice04] 

# tau=0.9
SIC_list90.lattice04 <- c()
GACV_list90.lattice04 <- c()
imputed_CV_list90.lattice04 <- c()
H_list_qt90.lattice04 <- list()

lambda_candidate90.lattice04 <- exp(seq(0,5, by=1))
for(lambda in lambda_candidate90.lattice04){
  cat("####### lambda = ", lambda, "####### \n")
  res90.lattice04 <- PIRLS(lambda=lambda, y=lattice04$f, tau=0.9, alpha=0.1, B=diag(n.lattice04), 
                           P=L.lattice04%*%L.lattice04, gamma_init=rep(mean(lattice04$f), n.lattice04), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list90.lattice04 <- c(SIC_list90.lattice04, res90.lattice04$SIC)
  GACV_list90.lattice04 <- c(GACV_list90.lattice04, res90.lattice04$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice04$f, tau=0.9, alpha=0.1, B=diag(n.lattice04), 
                          P=L.lattice04%*%L.lattice04, gamma_init=rep(mean(lattice04$f), n.lattice04), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice04, k=1, M=NULL, method="quantile")
  imputed_CV_list90.lattice04 <- c(imputed_CV_list90.lattice04, imputed_CV)
  H_list_qt90.lattice04[[length(H_list_qt90.lattice04)+1]] <- res90.lattice04$H
}
plot(lambda_candidate90.lattice04, SIC_list90.lattice04, type="l")
plot(lambda_candidate90.lattice04, GACV_list90.lattice04, type="l")
plot(lambda_candidate90.lattice04, imputed_CV_list90.lattice04, type="l")
opt.lambda.ind.SIC90.lattice04 <- which.min(SIC_list90.lattice04)
opt.lambda.SIC90.lattice04 <- lambda_candidate90.lattice04[opt.lambda.ind.SIC90.lattice04] 
opt.lambda.ind.GACV90.lattice04 <- which.min(GACV_list90.lattice04)
opt.lambda.GACV90.lattice04 <- lambda_candidate90.lattice04[opt.lambda.ind.GACV90.lattice04] 
opt.lambda.ind.imputedCV90.lattice04 <- which.min(imputed_CV_list90.lattice04)
opt.lambda.imputedCV90.lattice04 <- lambda_candidate90.lattice04[opt.lambda.ind.imputedCV90.lattice04] 


# tau=0.1
SIC_list10.lattice04 <- c()
GACV_list10.lattice04 <- c()
imputed_CV_list10.lattice04 <- c()
H_list_qt10.lattice04 <- list()

lambda_candidate10.lattice04 <- exp(seq(0,5, by=1))
for(lambda in lambda_candidate10.lattice04){
  cat("####### lambda = ", lambda, "####### \n")
  res10.lattice04 <- PIRLS(lambda=lambda, y=lattice04$f, tau=0.1, alpha=0.1, B=diag(n.lattice04), 
                           P=L.lattice04%*%L.lattice04, gamma_init=rep(mean(lattice04$f), n.lattice04), max_iterations=10000, tol=1e-6, check="Oh")
  SIC_list10.lattice04 <- c(SIC_list10.lattice04, res10.lattice04$SIC)
  GACV_list10.lattice04 <- c(GACV_list10.lattice04, res10.lattice04$GACV)
  imputed_CV <- imputedCV(lambda=lambda, y=lattice04$f, tau=0.1, alpha=0.1, B=diag(n.lattice04), 
                          P=L.lattice04%*%L.lattice04, gamma_init=rep(mean(lattice04$f), n.lattice04), max_iterations=10000, tol=1e-6, 
                          check="Oh", L=L.lattice04, k=1, M=NULL, method="quantile")
  imputed_CV_list10.lattice04 <- c(imputed_CV_list10.lattice04, imputed_CV)
  H_list_qt10.lattice04[[length(H_list_qt10.lattice04)+1]] <- res10.lattice04$H
}
plot(lambda_candidate10.lattice04, SIC_list10.lattice04, type="l")
plot(lambda_candidate10.lattice04, GACV_list10.lattice04, type="l")
plot(lambda_candidate10.lattice04, imputed_CV_list10.lattice04, type="l")
opt.lambda.ind.SIC10.lattice04 <- which.min(SIC_list10.lattice04)
opt.lambda.SIC10.lattice04 <- lambda_candidate10.lattice04[opt.lambda.ind.SIC10.lattice04] 
opt.lambda.ind.GACV10.lattice04 <- which.min(GACV_list10.lattice04)
opt.lambda.GACV10.lattice04 <- lambda_candidate10.lattice04[opt.lambda.ind.GACV10.lattice04] 
opt.lambda.ind.imputedCV10.lattice04 <- which.min(imputed_CV_list10.lattice04)
opt.lambda.imputedCV10.lattice04 <- lambda_candidate10.lattice04[opt.lambda.ind.imputedCV10.lattice04] 

# fitting results
# q1 <- plot_graph_custom3(lattice04, e.size=1.3, v.size=11, vertex_color = lattice04$f.true, 
#                          min=0, max=11, value="value",
#                          ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("true signal") # plot graph signal
q2 <- plot_graph_custom3(lattice04, e.size=1.3, v.size=11, vertex_color = lattice04$f, 
                         min=0, max=11, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("noisy signal") # plot graph signal
q3 <- plot_graph_custom3(lattice04, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_mean.lattice04[[opt.lambda.ind.imputedCV.mean.lattice04]]%*%lattice04$f), 
                         min=0, max=11, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("mean") # plot graph signal
q4 <- plot_graph_custom3(lattice04, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt90.lattice04[[opt.lambda.ind.imputedCV90.lattice04]]%*%lattice04$f), 
                         min=0, max=11, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.9") # plot graph signal
# q5 <- plot_graph_custom3(lattice04, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt50.lattice04[[opt.lambda.ind.imputedCV50.lattice04]]%*%lattice04$f), 
#                          min=0, max=11, value="value",
#                          ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.5") # plot graph signal
q6 <- plot_graph_custom3(lattice04, e.size=1.3, v.size=11, vertex_color = as.numeric(H_list_qt10.lattice04[[opt.lambda.ind.imputedCV10.lattice04]]%*%lattice04$f), 
                         min=0, max=11, value="value",
                         ratio=1, label.title.size=19, label.text.size = 17) + ggtitle("tau = 0.1") # plot graph signal

grid.arrange(q2,q3,q4,q6, ncol=4)

