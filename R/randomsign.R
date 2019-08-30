## assign random sign to correlations

.random_sign<- function(x){
  x * sample(c(-1,1), 1)
}

.add <- function(x) Reduce("+", x)
