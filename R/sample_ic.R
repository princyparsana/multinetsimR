## construct a sparse inverse covariance matrix
sample_ic <- function(p = 10, type = "scale-free"){
  if(type == "scale-free"){
    pg = igraph::as.undirected(igraph::sample_pa(p, directed = T))
  }else{
    pg = igraph::sample_gnm(p, m = p*3)
  }
  n_edges = igraph::ecount(pg)
  g <- igraph::as_adj(pg)
  g <- Matrix::triu(g)
  e_vals = runif(n_edges, min = c(-0.9, 0.5), max = c(-0.5, 0.9))
  g[g!=0] <- e_vals
  g <- convert_psd(g)
  list(graph = pg, theta = g)
}
