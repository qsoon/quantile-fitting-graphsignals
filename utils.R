library(gasper)
library(igraph)
library(ggplot2)

checkft <- function(x, tau){
  return(x*(tau-as.numeric(x<0)))
}

checkft_by_Zheng <- function(x, tau, alpha){
  return(tau*x + alpha*log(1+exp(-x/alpha)))
}

grad_checkft_by_Zheng <- function(x, tau, alpha){
  return(tau - 1/(1+exp(x/alpha)))
}

checkft_by_Oh <- Vectorize(function(x, tau, c){
  if((-c<=x) & (x<c)){
    return(0.5*(1-tau)*(x^2)/c + 0.5*(2*tau-1)*(x^2)/c*as.numeric(x>=0))
  }
  else{
    return((tau-1)*(x+0.5*c) + (x-(tau-0.5)*c)*as.numeric(x>=c))
  }
}, vectorize.args = "x")

grad_checkft_by_Oh <- Vectorize(function(x, tau, c){
  if((-c<=x) & (x<c)){
    return(tau*x/c*as.numeric(x >= 0) + (1-tau)*x/c*as.numeric(x < 0))
  }
  else{
    return(tau - as.numeric(x < -c))
  }
}, vectorize.args = "x")

checkft_by_Nychka <- Vectorize(function(x, tau, delta){
  if((-delta<x) & (x<delta)){
    return((tau*as.numeric(x>0) + (1-tau)*as.numeric(x<0))*x^2/delta)
  }
  else{
    return(checkft(x, tau))
  }
}, vectorize.args = "x")


# Define the objective function to minimize
objective_function <- function(y, x, tau, alpha, L, lambda, type="L2", check) {
  if (type=="global"){
    penalty <- as.numeric(t(x)%*%L%*%x)
  } else if (type=="L1"){
    penalty <- sum(abs(as.numeric(L%*%x)))
  } else{
    penalty <- as.numeric(t(x)%*%L%*%L%*%x)
  }
  
  if(check=="Zheng"){
    return(mean(checkft_by_Zheng(y-x, tau, alpha)) + lambda*penalty)
  } else if(check=="Oh"){
    return(mean(checkft_by_Oh(y-x, tau, alpha)) + lambda*penalty)
  }
}

# Gradient of the objective function
gradient <- function(y, x, tau, alpha, L, lambda, type="L2", check) {
  n <- length(y)
  if (type=="global"){
    grad_penalty <- as.numeric(2*L%*%x)
  } else if (type=="L2"){
    grad_penalty <- as.numeric(2*L%*%L%*%x)
  } 
  
  if(check=="Zheng"){
    return(-grad_checkft_by_Zheng(y-x, tau, alpha) / n + lambda*grad_penalty)
  } else if(check=="Oh"){
    return(-grad_checkft_by_Oh(y-x, tau, alpha) / n + lambda*grad_penalty)
  }
  
}

# Gradient Descent Function
gradient_descent <- function(learning_rate, max_iterations, tol=1e-6, y, initial_x, tau, alpha, L, lambda, type="L2", check) {
  x <- initial_x
  val <- c()
  for (iteration in 1:max_iterations) {
    gradient_value <- as.numeric(gradient(y, x, tau, alpha, L, lambda, type, check))
    x <- x - learning_rate * gradient_value
    val0 <- objective_function(y, x, tau, alpha, L, lambda, type, check)
    val <- c(val, val0)
    cat("Iteration:", iteration, "Value:", val0, 
        "Gradient val:", norm(gradient_value, type="2"), "\n")
    
    if (norm(gradient_value, type="2") < tol) {
      cat("Converged to a minimum.\n")
      break
    }
  }
  
  return(list(x, val))
}

plot_gs <- function(G=NULL, W=NULL, loc=NULL, Y=NULL, igraph_layout=2, igraph_params=list(), 
                    vsize=0.75, esize=3, alpha=1, limits = NULL, nplot=4, signal=FALSE,
                    wtoption="alpha", nodecol = "Spectral", edgecol = "gray", label=FALSE,
                    nudge_x = 0.2, nudge_y = 0){
  if(is(G, "igraph")){
    if(is.null(W)){
      if(is_weighted(G)){
        W <- igraph::as_adjacency_matrix(G, attr="weight") # extract edge weight matrix
      } else{
        W <- igraph::as_adjacency_matrix(G) # extract adjacency matrix
      }
    }
    if(is.null(loc)){
      if(igraph_layout == 1){
        loc <- igraph::layout.reingold.tilford(G, params=igraph_params)
      } else if(igraph_layout == 2){
        loc <- igraph::layout.circle(G, params=igraph_params)
      } else if(igraph_layout == 3){
        loc <- igraph::layout.sphere(G, params=igraph_params)
      } else if(igraph_layout == 4){
        loc <- igraph::layout.random(G, params=igraph_params)
      } else if(igraph_layout == 5){
        loc <- igraph::layout.fruchterman.reingold(G, params=igraph_params)
      } else if(igraph_layout == 6){
        loc <- igraph::layout.kamada.kawai(G, params=igraph_params)
      }
    }
  } else if(is.null(W)){
    stop("'W' should be given unless igraph object G is given")
  } else if((!is(W, "matrix")) & (!is(W, "sparseMatrix"))){
    stop("'W' must be matrix or sparseMatrix")
  } else if(is.null(loc)){
    G <- graph.adjacency(W,mode="undirected",weighted=TRUE) # create igraph object
    if(igraph_layout == 1){
      loc <- igraph::layout.reingold.tilford(G, params=igraph_params)
    } else if(igraph_layout == 2){
      loc <- igraph::layout.circle(G, params=igraph_params)
    } else if(igraph_layout == 3){
      loc <- igraph::layout.sphere(G, params=igraph_params)
    } else if(igraph_layout == 4){
      loc <- igraph::layout.random(G, params=igraph_params)
    } else if(igraph_layout == 5){
      loc <- igraph::layout.fruchterman.reingold(G, params=igraph_params)
    } else if(igraph_layout == 6){
      loc <- igraph::layout.kamada.kawai(G, params=igraph_params)
    }
  }
  
  if(!is(W, "sparseMatrix")){
    W <- Matrix::Matrix(W, sparse = TRUE)
  }
  
  res <- list()
  res$W <- W
  res$loc <- loc
  res$Y <- Y
  
  triplet <- Matrix::mat2triplet(W)
  W <- data.frame(i=triplet$i, j=triplet$j, x=triplet$x)
  
  x <- loc[, 1]
  y <- loc[, 2]
  ind_i <- W[, 1]
  ind_j <- W[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df1$label <- 1:nrow(df1)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  
  if(signal==FALSE){
    if(wtoption=="alpha"){
      p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y,
                                                      xend = y1, yend = y2), alpha=alpha*W[,3]/max(W[,3]), linewidth=esize,
                                                  color = edgecol, data = df2, show.legend = FALSE) +
        geom_point(size = vsize) +theme_void() # size aesthetic for lines was deprecated in ggplot2 3.4.0, use linewidth instead
      
      if(label==TRUE){
        p1 <- p1 + geom_text(label=df1$label,nudge_x = nudge_x, nudge_y = nudge_y,
                             check_overlap = T)
      }
      res$graph <- p1
    } else if(wtoption=="width"){
      p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y,
                                                      xend = y1, yend = y2), alpha=alpha, linewidth=esize*W[,3]/max(W[,3]),
                                                  color = edgecol, data = df2, show.legend = FALSE) +
        geom_point(size = vsize) + theme_void()
      
      if(label==TRUE){
        p1 <- p1 + geom_text(label=df1$label,nudge_x = nudge_x, nudge_y = nudge_y,
                             check_overlap = T)
      }
      res$graph <- p1
    } 
  }
  
  
  if(signal==TRUE){
    # tau <- nrow(Y)-1
    # plot.ind <- round(seq(0, tau, length.out=nplot))
    # plt.list <- list()
    # legend.titles <- paste("tau=",plot.ind, sep="")
    if(is.null(limits)){
      limits <- range(Y)
    } 
    if(wtoption=="alpha"){
      f <- Y
      p2<-ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2), alpha=alpha*W[,3]/max(W[,3]), linewidth=esize,
                                                color = edgecol, data = df2, show.legend = FALSE) + 
        geom_point(size = vsize, aes(colour = f)) +
        scale_colour_distiller(palette = nodecol, limits = limits) + 
        labs(colour = "") +
        theme_void() + theme(legend.text = element_text(size = 8), 
                             legend.key.size = unit(0.8, "line"), 
                             legend.margin = margin(0,0, 0, 0), plot.margin = unit(c(0.01, 0.01, 0.01, 
                                                                                     0.01), "cm")) 
      
      res$graphsignal <- p2
    } else if(wtoption=="width"){
      f <- Y
      p2<-ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2), alpha=alpha, linewidth=esize*W[,3]/max(W[,3]),
                                                color = edgecol, data = df2, show.legend = FALSE) + 
        geom_point(size = vsize, aes(colour = f)) +
        scale_colour_distiller(palette = nodecol, limits = limits) + 
        labs(colour = "") +
        theme_void() + theme(legend.text = element_text(size = 8), 
                             legend.key.size = unit(0.8, "line"), 
                             legend.margin = margin(0,0, 0, 0), plot.margin = unit(c(0.01, 0.01, 0.01, 
                                                                                     0.01), "cm")) 
      
      res$graphsignal <- p2
    }
    
  }
  
  
  # if(signal){
  #   print(p2)
  # } else{
  #   print(p1)
  # }     
  return(res)
}

#######################
#### visualization ####
#######################
plot_graph_custom2 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                  xend = y1, yend = y2, color=w), linewidth=e.size, 
                                              data = df2) +
    scale_color_gradient(low="grey", high="black", name="Edge Weight")+
    geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
    # scale_fill_gradient(low="blue", high="red", na.value = "gray", name = "Temperature") +
    scale_fill_gradientn(colors = rev(rainbow(7)[-7]), limits=c(min, max), na.value = "gray", name = value) + 
    geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
    theme_void() +
    guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) +
    theme(legend.margin = margin(10,10,10,10), plot.margin = margin(10,10,10,10), plot.title=element_text(size=30, face="bold", hjust = 0.5))
  print(p1)
}

plot_graph_custom3 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2), color = "gray", data = df2) + 
      geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
    print(p1)
  }

  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) +
      theme(legend.margin = margin(10,10,10,10), plot.margin = margin(10,10,10,10), 
            plot.title=element_text(size=30, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

PIRLS <- function(lambda, y, tau, alpha, B, P, gamma_init, max_iterations, tol, check){
  n <- length(y)
  gamma <- gamma_init
  for (iteration in 1:max_iterations){
    gamma.old <- gamma
    if(check=="Zheng"){
      W <- diag(grad_checkft_by_Zheng(as.numeric(y-B%*%gamma), tau, alpha) / as.numeric((2*(y-B%*%gamma))))
    } else if(check=="Oh"){
      W <- diag(grad_checkft_by_Oh(as.numeric(y-B%*%gamma), tau, alpha) / as.numeric((2*(y-B%*%gamma))))
    }

    tmp <- solve(t(B)%*%W%*%B + n*lambda*P)%*%t(B)%*%W
    gamma <- tmp%*%y
    H <- B%*%tmp
    yhat <- as.numeric(H%*%y)
    
    if (norm(gamma-gamma.old, type="2") < tol) {
      cat("Iteration:", iteration, "\n")
      cat("Converged!\n")
      break
    }
  }
  
  res <- list(yhat = yhat, H=H, SIC=SIC(lambda, y, tau, alpha, H, check), GACV=GACV(lambda, y, tau, alpha, H, check))
  
  return(res)
}

SIC <- function(lambda, y, tau, alpha, H, check){
  n <- length(y)
  yhat <- as.numeric(H%*%y)
  df <- sum(diag(H))
  # return(log(mean(checkft(y-yhat, tau))) + log(n)/(2*n)*df)
  if(check=="Zheng"){
    return(log(mean(checkft_by_Zheng(y-yhat, tau, alpha))) + log(n)/(2*n)*df)
  } else if(check=="Oh"){
    return(log(mean(checkft_by_Oh(y-yhat, tau, alpha))) + log(n)/(2*n)*df)
  }
}

GACV <- function(lambda, y, tau, alpha, H, check){
  n <- length(y)
  yhat <- as.numeric(H%*%y)
  # return(sum(checkft(y-yhat, tau) / (n-sum(diag(H)))))
  if(check=="Zheng"){
    return(sum(checkft_by_Zheng(y-yhat, tau, alpha) / (n-sum(diag(H)))))
  } else if(check=="Oh"){
    return(sum(checkft_by_Oh(y-yhat, tau, alpha) / (n-sum(diag(H)))))
  }
}

GCV <- function(lambda, y, H){
  n <- length(y)
  yhat <- as.numeric(H%*%y)
  gcv <- mean((y-yhat)^2) / (1-sum(diag(H))/n)^2
  return(gcv)
}

imputegs <- function(Laplacian, signal, vertices){
  n <- dim(Laplacian)[1]
  s.known <- vertices
  s.unknown <- setdiff(1:n, s.known)
  
  ordered_Laplacian <- Laplacian[c(s.known, s.unknown), c(s.known, s.unknown)]
  ordered_signal <- signal[c(s.known, s.unknown)]
  
  RT <- ordered_Laplacian[(length(s.known)+1):n,1:length(s.known)]
  Lu <- ordered_Laplacian[(length(s.known)+1):n,(length(s.known)+1):n]
  sb <- ordered_signal[1:length(s.known)] #signal_known
  
  su <- solve(Lu, -RT %*% sb)
  
  signal[s.unknown] <- su
  
  return(signal)
}

imputedCV <- function(lambda, y, tau, alpha, B, P, gamma_init, max_iterations, tol, check, L, k, M=NULL, method){
  if(method=="mean"){ # for mean denoising
    n <- length(y)
    if(k==1){ # leave one out
      tmp <- rep(0,n)
      for(i in 1:n){
        z <- y
        imputed <- imputegs(L, y, c(1:n)[-i])
        z[i] <- imputed[i]
        H <- solve(diag(1,n)+n*lambda*L%*%L)
        tmp[i] <- (y[i]-as.numeric(H%*%z)[i])^2
      }
      res <- mean(tmp)
    } else{ # randomly chosen k
      tmp <- rep(0,M)
      for(i in 1:M){
        z <- y
        removed <- sample(n,k)
        imputed <- imputegs(L, y, c(1:n)[-removed])
        z[removed] <- imputed[removed]
        H <- solve(diag(1,n)+n*lambda*L%*%L)
        tmp[i] <- mean((y[removed]-as.numeric(H%*%z)[removed])^2)
      }
      res <- mean(tmp)
    }
  } else if(method=="quantile"){ # for quantile analysis
    n <- length(y)
    if(k==1){ # leave one out
      tmp <- rep(0,n)
      for(i in 1:n){
        z <- y
        imputed <- imputegs(L, y, c(1:n)[-i])
        z[i] <- imputed[i]
        res <- PIRLS(lambda, z, tau, alpha, B, P, gamma_init, max_iterations, tol, check)
        H <- res$H
        if(check=="Zheng"){
          tmp[i] <- checkft_by_Zheng(y[i]-as.numeric(H%*%z)[i], tau, alpha)
        } else if(check=="Oh"){
          tmp[i] <- checkft_by_Oh(y[i]-as.numeric(H%*%z)[i], tau, alpha)
        }
      }
      res <- mean(tmp)
    } else{ # randomly chosen k
      tmp <- rep(0,M)
      for(i in 1:M){
        z <- y
        removed <- sample(n,k)
        imputed <- imputegs(L, y, c(1:n)[-removed])
        z[removed] <- imputed[removed]
        res <- PIRLS(lambda, z, tau, alpha, B, P, gamma_init, max_iterations, tol, check)
        H <- res$H
        if(check=="Zheng"){
          tmp[i] <- mean(checkft_by_Zheng(y[removed]-as.numeric(H%*%z)[removed], tau, alpha))
        } else if(check=="Oh"){
          tmp[i] <- mean(checkft_by_Oh(y[removed]-as.numeric(H%*%z)[removed], tau, alpha))
        }
      }
      res <- mean(tmp)
    }
  }
  
  return(res)
}

pseudoimputedCV <- function(lambda, y, tau, alpha, B, P, gamma_init, max_iterations, tol, check, L, k, M=NULL, method){
  if(method=="mean"){ # for mean denoising
    res <- PIRLS(lambda, y, tau, alpha, B, P, gamma_init, max_iterations, tol, check)
    H <- res$H
    emp.pseudodata <- y
    n <- length(y)
    if(k==1){ # leave one out
      tmp <- rep(0,n)
      for(i in 1:n){
        z <- emp.pseudodata
        imputed <- imputegs(L, emp.pseudodata, c(1:n)[-i])
        z[i] <- imputed[i]
        H <- solve(diag(1,n)+n*lambda*L%*%L)
        tmp[i] <- (y[i]-as.numeric(H%*%z)[i])^2
      }
      res <- mean(tmp)
    } else{ # randomly chosen k
      tmp <- rep(0,M)
      for(i in 1:M){
        z <- emp.pseudodata
        removed <- sample(n,k)
        imputed <- imputegs(L, emp.pseudodata, c(1:n)[-removed])
        z[removed] <- imputed[removed]
        H <- solve(diag(1,n)+n*lambda*L%*%L)
        tmp[i] <- mean((y[removed]-as.numeric(H%*%z)[removed])^2)
      }
      res <- mean(tmp)
    }
  } else if(method=="quantile"){ # for quantile analysis
    res <- PIRLS(lambda, y, tau, alpha, B, P, gamma_init, max_iterations, tol, check)
    H <- res$H
    fhat <- as.numeric(H%*%y)
    if(check=="Zheng"){
      emp.pseudodata <- fhat + grad_checkft_by_Zheng(y-fhat, tau, alpha) / 2
    } else if(check=="Oh"){
      emp.pseudodata <- fhat + grad_checkft_by_Oh(y-fhat, tau, alpha) / 2
    }
    n <- length(y)
    if(k==1){ # leave one out
      tmp <- rep(0,n)
      for(i in 1:n){
        z <- emp.pseudodata
        imputed <- imputegs(L, emp.pseudodata, c(1:n)[-i])
        z[i] <- imputed[i]
        res <- PIRLS(lambda, z, tau, alpha, B, P, gamma_init, max_iterations, tol, check)
        H <- res$H
        if(check=="Zheng"){
          tmp[i] <- checkft_by_Zheng(y[i]-as.numeric(H%*%z)[i], tau, alpha)
        } else if(check=="Oh"){
          tmp[i] <- checkft_by_Oh(y[i]-as.numeric(H%*%z)[i], tau, alpha)
        }
      }
      res <- mean(tmp)
    } else{ # randomly chosen k
      tmp <- rep(0,M)
      for(i in 1:M){
        z <- emp.pseudodata
        removed <- sample(n,k)
        imputed <- imputegs(L, emp.pseudodata, c(1:n)[-removed])
        z[removed] <- imputed[removed]
        res <- PIRLS(lambda, z, tau, alpha, B, P, gamma_init, max_iterations, tol, check)
        H <- res$H
        if(check=="Zheng"){
          tmp[i] <- mean(checkft_by_Zheng(y[removed]-as.numeric(H%*%z)[removed], tau, alpha))
        } else if(check=="Oh"){
          tmp[i] <- mean(checkft_by_Oh(y[removed]-as.numeric(H%*%z)[removed], tau, alpha))
        }
      }
      res <- mean(tmp)
    }
  }
}


mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}
