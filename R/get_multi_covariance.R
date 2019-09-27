# Code to simulate undirected graphs with shared and specific edges

library(igraph)

# ################# functions ####################
# random_sign<- function(x){
#   x * sample(c(-1,1), 1)
# }
#
# add <- function(x) Reduce("+", x)
###############################################
get_multi_covariance <- function(p = 50, num_nets = 3, prop_sample = 0.7, v = 0.3, u = 0.1, ...){
  # p: number of nodes in the graphs
  # num_nets: number of related networks to generate
  # prop_sample: sample this proportion of total edges in the network
  if(!exists('n')){
    n = p * 10 # if number of samples is not provided, the number of samples will be 10 times the number of nodes
  }

  ## temporory --> update with more efficient version
  pg = igraph::union(sample_pa(p, directed = T))

  master_edgelist <- as_edgelist(pg)

  multi_g <- lapply(1:num_nets, function(each){
    subgraph.edges(pg, sample(1:ecount(pg), round(prop_sample * ecount(pg))), delete.vertices = FALSE)
  })

  # print(sapply(multi_g, function(x) fit_power_law(degree(x))))
  theta_cov_mat <- lapply(multi_g, function(g, n, p){
    sorted_p <- as.numeric(topo_sort(g))
    neighbors_sorted_p <- sapply(sorted_p, function(x) as.numeric(neighbors(g, x, mode = "in")))

    ## simulate uncorrelated data from mvn for nodes with zero neighbors
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
        # dat.sim[,sorted_p[i]] <- rowSums(sapply(as.data.frame(neighbor_dat), .random_sign)) + (rnorm(n))
        dat.sim[,sorted_p[i]] <- rowSums(as.data.frame(neighbor_dat)) + (rnorm(n))
      }
    }
    theta <- cor(dat.sim)
    diag(theta) = 0
    omega = theta * v
    diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
    sigma = cov2cor(solve(omega))
    list(precision = omega, covariance = sigma, graph = as.undirected(g))
  }, n = n, p = p)

  theta_cov_mat
}

