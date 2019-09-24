library(igraph)

get_GMRF_network <- function(p = 100, v = 0.3, u = 0.1){
  p<-100
  if(!exists('n')){
    n = p * 100 # if number of samples is not provided, the number of samples will be 100 times the number of nodes
  }
  g = sample_pa(p, directed = T)
  master_edgelist <- as_edgelist(g)
  sorted_p <- as.numeric(topo_sort(g))
  neighbors_sorted_p <- sapply(sorted_p, function(x) as.numeric(neighbors(g, x, mode = "in")))
  nodes_with_zero_neighbors <- sorted_p[sapply(neighbors_sorted_p, length)==0]
  dat.sim <- matrix(NA, nrow = n, ncol = p)
  dat.sim[,nodes_with_zero_neighbors] <- MASS::mvrnorm(n,
                                                       mu = rep(0, length(nodes_with_zero_neighbors)),
                                                       Sigma = diag(length(nodes_with_zero_neighbors)),
                                                       empirical = T)

  for(i in 1:p){
    if(sorted_p[i] %in% nodes_with_zero_neighbors){
      next
    }else{
      neighbor_dat <- scale(dat.sim[,neighbors_sorted_p[[i]]])
      dat.sim[,sorted_p[i]] <- rowSums(sapply(as.data.frame(neighbor_dat), .random_sign))+ (rnorm(n))
      # dat.sim[,sorted_p[i]] <- rowSums(as.data.frame(neighbor_dat)) + (rnorm(n))
    }
  }
  theta <- cor(dat.sim)
  diag(theta) = 0
  omega = theta * v
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = cov2cor(solve(omega))
  return(list(precision = omega, covariance = sigma, graph = g))
}
