library(igraph)
library(netsimulatR)
set.seed(10)
dat_sim <- get_multi_covariance(p=100, num_nets = 3, n= 300)

covariance_mat <- lapply(dat_sim, function(x){
  x$covariance
})

mixture_weights <- list(d1 = c(0.6, 0.3, 0.1),
                        d2 = c(0.7, 0.1, 0.2),
                        d3 = c(0.4, 0.5, 0.1),
                        d4 = c(0.1, 0.7, 0.2),
                        d5 = c(0.3, 0.1, 0.6))

gen_empirical_corr = lapply(mixture_weights, function(eachmix){
  dat_tmp = (eachmix[1] *  mvtnorm::rmvnorm(300, sigma = covariance_mat[[1]])) +
    (eachmix[2] *  mvtnorm::rmvnorm(300, sigma = covariance_mat[[2]])) +
    (eachmix[3] * mvtnorm::rmvnorm(300, sigma = covariance_mat[[3]]))
  cov(scale(dat_tmp))
})


dist_norm2 = sapply(gen_empirical_corr, function(x){
  sapply(gen_empirical_corr, function(a,b){
    norm(a-b, type = 'f')
  }, x)
})

diag(dist_norm2) <- NA
gplots::heatmap.2(dist_norm2, trace = "none")
dist_norm2

################################################
covariance_mat <- lapply(dat_sim, function(x){
  x$precision
})


dat_cov <- lapply(mixture_weights, function(eachmix, covlist){
  mix_list = mapply(function(x,y){
    x*y
  }, eachmix, covlist, SIMPLIFY = F)
  solve(netsimulatR:::.add(mix_list))
}, covlist = covariance_mat)

gen_empirical_covariance = lapply(dat_cov, function(x){
  dat_tmp = mvtnorm::rmvnorm(300, sigma = x) #+ (mvtnorm::rmvnorm(100, sigma = diag(500)))
  solve(cov(scale(dat_tmp)))
})


dist_norm2 = sapply(gen_empirical_covariance, function(x){
  sapply(gen_empirical_covariance, function(a,b){
    a[abs(a)<0.05] <- 0
    b[abs(b)<0.05] <-0
    k <- a-b
    norm(k, type = 'f')
  }, x)
})

diag(dist_norm2) <- NA
gplots::heatmap.2(dist_norm2, trace = "none")
dist_norm2
