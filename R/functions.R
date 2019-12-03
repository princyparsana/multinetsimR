

convert_psd <- function(thisgraph){
  thisgraph = apply(thisgraph, 1, function(x){
    if(sum(abs(x))!=0){
      scaled_val = x/(1.5*sum(abs(x)))
    }else{
      x
    }
  })
  diag(thisgraph) <- 1
  thisgraph = (thisgraph + t(thisgraph))/2
  thisgraph
}

## construct a master sparse inverse covariance matrix

generate_master <- function(p = 10){
  pg = as.undirected(sample_pa(p, directed = T))
  n_edges = 2*ecount(pg)
  g <- as_adj(pg)
  e_vals = runif(n_edges, min = c(-1, 0.5), max = c(-0.5, 1))
  g[g!=0] <- e_vals
  g <- convert_psd(g)
  list(graph = pg, theta = g)
}

generate_multi_ic <- function(parent_graph, parent_theta, s = 0.5, num_nets = 3){
  full_g = graph_from_edgelist(t(combn(seq(1, vcount(parent_graph)),2)), directed = F)
  diff_parent_graph = difference(full_g, parent_graph)
  n_edges = 2*ecount(parent_graph)
  subsample_count = round(s*(n_edges/2))
  lapply(1:num_nets, function(x,y){
    g_sub_delete = subgraph.edges(y, sample(1:(n_edges/2), subsample_count), delete.vertices = F)
    update_gg = difference(y, g_sub_delete)
    g_sub_new = subgraph.edges(diff_parent_graph, sample(1:(n_edges/2), subsample_count), delete.vertices = F)
    g_sub = igraph::union(update_gg, g_sub_new)
    g_new = as_adj(g_sub)
    this_n_edges <- ecount(g_sub)*2
    e_vals = runif(this_n_edges, min = c(-1, 0.5), max = c(-0.5, 1))
    g_new[g_new!=0] = e_vals
    g_task = parent_theta * as_adj(update_gg)
    g_task= g_new + g_task
    convert_psd(g_task)
  }, y = parent_graph)
}

