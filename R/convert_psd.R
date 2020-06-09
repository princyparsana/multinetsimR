library(Matrix)
library(igraph)

convert_psd <- function(thisgraph, psd_type = "jgl", v  = 0.3, u = 0.1){
  thisgraph <- as.matrix(thisgraph)
  if(psd_type == "jgl"){
    rowsum_graph = rowSums(abs(as.matrix(thisgraph))) * 1.5
    for(i in 1:length(rowsum_graph)){
      thisgraph[i,] <- thisgraph[i,]/rowsum_graph[i]
    }
    #   if(sum(abs(x))!=0){
    #     scaled_val = x/(1.5*sum(abs(x)))
    #   }else{
    #     x
    #   }
    # })
    thisgraph[is.nan(thisgraph)]<-0

    diag(thisgraph) <- 1.00
    thisgraph = (thisgraph + t(thisgraph))/2
    inv_thisgraph <- solve(thisgraph)
    denominator_cov <- tcrossprod(diag(inv_thisgraph))
    d_mat <- matrix(0.6, nrow = nrow(inv_thisgraph), ncol = ncol(inv_thisgraph))
    diag(d_mat) <- 1
    cov_mat <- d_mat/sqrt(denominator_cov)
    thisgrpah = solve(cov_mat)
  }

  if(psd_type == "huge"){
    ## adapted from huge package, allows to control the magnitude of PSD precision matrix
    thisgraph  = thisgraph + t(thisgraph)
    thisgraph[thisgraph!=0] <- 1
    diag(thisgraph) <- 0
    theta = thisgraph*v
    diag(theta) = abs(min(eigen(theta)$values)) + 0.1 + u
    sigma = cov2cor(solve(theta))
    thisgraph = solve(sigma)
  }

  thisgraph
}
