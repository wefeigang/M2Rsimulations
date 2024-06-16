library(glasso) # glasso package
library(MASS) # contains mvrnorm
library(igraph) # plotting library

set.seed(10)

# Make sparse precision matrix
sparse_matrix <- function(p, k){ # k is the number from paper by Luan AR(k)
  true_prec <- diag(1, p)
  if (k == 1){
    true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.5
  }
  if (k == 2){
    true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.5
    true_prec[abs(row(true_prec) - col(true_prec)) == 2] <- 0.25
  }
  if (k == 3){
    true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.4
    true_prec[abs(row(true_prec) - col(true_prec)) == 2] <- 0.2
    true_prec[abs(row(true_prec) - col(true_prec)) == 3] <- 0.2
  }
  return(true_prec)
}

# Simulation function
glasso_simulation <- function(n, p, lambda, k){
  # Forming a tridiagonal precision matrix
  true_prec <- sparse_matrix(p = p, k = k)
  true_cov <- solve(true_prec)
  data <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
  data <- scale(data)
  # should I scale the data?
  #est_prec <- CVglasso(X = data)$Omega
  #est_cov <- CVglasso(X = data)$Sigma
  est_prec <- glasso(s = cov(data), rho = lambda)$wi
  return(est_prec)
}

CV_glasso_simulation <- function(n, p){
  # Forming a tridiagonal precision matrix
  true_prec <- sparse_matrix(p)
  true_cov <- round(solve(true_prec), 5) # round to 5 d.p. to get rid of slight inaccuracies when inverting matrix
  data <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
  # should I scale the data?
  est_prec <- CVglasso(X = data)$Omega
  est_cov <- CVglasso(X = data)$Sigma
  return(est_prec)
}

n = 1000
p = 10
k = 1

true_prec <- sparse_matrix(p, k)

#Plotting
par(mfrow=c(3, 7), mar=c(3,1,3,1)) # set up graphs (mar is the margins)

# To make labels outside circle
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

lab.locs <- radian.rescale(x=1:p, direction=-1, start=0)

# True network
# adjacency_matrix <- (true_prec != 0) * 1
# network <- graph_from_adjacency_matrix(adjacency_matrix , mode = 'undirected', diag = F )
# plot(network, layout = layout.circle, 
#      vertex.label.degree = lab.locs, vertex.label.dist = 2, vertex.label.cex = 2,
#      vertex.size = 10, vertex.color = "white", vertex.frame.color = "black", 
#      edge.color = "black", vertex.label.color = "black", edge.width = 2, vertex.frame.width = 2)

# CV network
# true_cov <- solve(true_prec)
# data <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
# est_lambda <- cv_glasso(X = data, lambda_seq = seq(0.001, 1, 0.001))$best_lambda
# est_prec <- glasso_simulation(n = n, p = p, lambda = CVglasso(X = data)$Tuning[2])
# adjacency_matrix <- (est_prec != 0) * 1
# network <- graph_from_adjacency_matrix(adjacency_matrix , mode = 'undirected', diag = F )
# plot(network, layout = layout.circle,
#      vertex.label.degree = lab.locs, vertex.label.dist = 2.8, vertex.label.cex = 1,
#      vertex.size = 10, vertex.color = "white", vertex.frame.color = "black",
#      edge.color = "black", vertex.label.color = "black")

# Plotting networks for varying lambda
set.seed(123)

for (k in 1:3){
  true_prec <- sparse_matrix(p = p, k = k)
  true_cov <- solve(true_prec)
  data <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
  S <- cov(data)
  diag(S) <- 0
  lam.max <- max(abs(S))
  lam.min <- ifelse(nrow(X) > ncol(X), 10^-4, 0.01)
  lambda_seq <- 10^seq(log10(lam.min), log10(lam.max), length = 56)
  for (i in 41:47){
    lambda <- lambda_seq[i]
    est_prec <- glasso_simulation(n = n, p = p, lambda = lambda, k = k)
    adjacency_matrix <- (est_prec != 0) * 1
    network <- graph_from_adjacency_matrix(adjacency_matrix , mode = 'undirected', diag = F )
    if (i == 44){
      plot(network, layout = layout.circle,
           vertex.label.degree = lab.locs, vertex.label.dist = 4, vertex.label.cex = 2, 
           vertex.size = 20, vertex.color = "white", vertex.frame.color = "black",
           edge.color = "black", vertex.label.color = "black", edge.width = 2, vertex.frame.width = 2)
      title(paste0("Model ", k), cex.main = 2)
      mtext(paste0('lambda = ', round(lambda, 4)), side=1, line=1.7, cex=1) # text below plot
    }else{
      plot(network, layout = layout.circle,
           vertex.label.degree = lab.locs, vertex.label.dist = 4, vertex.label.cex = 2, 
           vertex.size = 20, vertex.color = "white", vertex.frame.color = "black",
           edge.color = "black", vertex.label.color = "black", edge.width = 2, vertex.frame.width = 2)
      mtext(paste0('lambda = ', round(lambda, 4)), side=1, line=1.7, cex=1) # text below plot
    }
  }
}

# par(mfrow=c(1, 3), mar=c(3,3,3,3))
# 
# for (k in 1:3){
#   true_prec <- sparse_matrix(p, k)
#   
#   # True network
#   adjacency_matrix <- (true_prec != 0) * 1
#   network <- graph_from_adjacency_matrix(adjacency_matrix , mode = 'undirected', diag = F )
#   plot(network, layout = layout.circle, 
#        vertex.label.degree = lab.locs, vertex.label.dist = 2.5, vertex.label.cex = 1.5,
#        vertex.size = 10, vertex.color = "white", vertex.frame.color = "black", 
#        edge.color = "black", vertex.label.color = "black", edge.width = 2, vertex.frame.width = 2,
#        main = paste0("Model", k))
# }