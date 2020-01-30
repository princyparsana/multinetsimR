## construct a sparse inverse covariance matrix
sample_ic <- function(p = 10, type = "scale-free" ,min_support = c(-0.4, 0.1), max_support = c(0.1, 0.4), seed = 1){
  set.seed(seed)
  if(type == "scale-free"){
    pg = igraph::as.undirected(igraph::sample_pa(p, directed = T))
  }else{
    pg = igraph::sample_gnm(p, m = p*3)
  }
  n_edges = igraph::ecount(pg)
  g <- igraph::as_adj(pg)
  g <- Matrix::triu(g)
  e_vals = runif(n_edges, min = c(-0.4, 0.1), max = c(0.1, 0.4))
  g[g!=0] <- e_vals
  g <- convert_psd(g)
  if(!all(eigen(g)$values >=0)){
    cat("PSD check unsuccessful - negative eigenvalues")
  }else{
    cat("PSD check was successful - all eigenvalues are greater than or equal to 0")
  }
  list(graph = pg, theta = g)
}
