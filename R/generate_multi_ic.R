## construct multiple sparse inverse covariance matrix

generate_multi_ic <- function(parent_graph, parent_theta, s = 0.5, num_nets = 3){
  full_g = graph_from_edgelist(t(combn(seq(1, vcount(parent_graph)),2)), directed = F)
  diff_parent_graph = difference(full_g, parent_graph, delete.vertices = F)
  n_edges = ecount(parent_graph)
  subsample_count = round(s*(n_edges))
  lapply(1:num_nets, function(x,y){
    g_sub_delete = subgraph.edges(y, sample(1:n_edges, subsample_count), delete.vertices = F)
    update_gg = difference(y, g_sub_delete)
    g_sub_new = subgraph.edges(diff_parent_graph, sample(1:ecount(diff_parent_graph), subsample_count), delete.vertices = F)
    g_sub = igraph::union(update_gg, g_sub_new)
    g_new = as_adj(g_sub)
    g_new = triu(g_new)
    this_n_edges <- ecount(g_sub)
    e_vals = runif(this_n_edges, min = c(-0.5, 0.1), max = c(-0.1, 0.5))
    g_new[g_new!=0] = e_vals
    parent_theta = triu(parent_theta)
    g_task = parent_theta * as_adj(update_gg)
    g_task= g_new + g_task
    convert_psd(g_task)
  }, y = parent_graph)
}

