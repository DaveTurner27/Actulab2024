sinistre <- read.csv("nath_perte")$x
boost <- read.csv("nath_boosting")$x
glm <- read.csv("nath_glm")$x

primeChoisi <- function(v){
  p1 <- min(v)
  p2 <- min(v[which(v =! p1)])
  if(min(1, -0.504244*(exp(19.11*(p1 - p2)/p2)-1) +0.5) > runif(1)) {
    return(p1)
  }
  p2
}
moy <- rep(mean(sinistre), 747628)
pi <- data.frame("boost" = moy*1.15,
                 "b-g" = glm*1.20,
                 "glm" = boost*1.20)


risoto <- t(apply(pi, MARGIN = 1, FUN = function(x) {
  idx <- sample(1:3, 2)
  x == primeChoisi(x[idx])
}))

colMeans(risoto)

sum(pi[risoto[,1], 1] - sinistre[risoto[,1]])/1000000
sum(pi[risoto[,2], 2] - sinistre[risoto[,2]])/1000000
sum(pi[risoto[,3], 3] - sinistre[risoto[,3]])/1000000

mean(pi[risoto[,1], 1] - sinistre[risoto[,1]])
mean(pi[risoto[,2], 2] - sinistre[risoto[,2]])
mean(pi[risoto[,3], 3] - sinistre[risoto[,3]])
