Dis <- function (cl, x) {
  # Dis : Total separation between clusters
  x <- as.matrix(x)
  k <- max(cl)
  centers.matrix <- centers(cl, x)
  Distance.centers <- dist(centers.matrix)
  Dmin <- min(Distance.centers)
  Dmax <- max(Distance.centers)
  Distance.centers <- as.matrix(Distance.centers)

  s2 <- 0
  for (u in 1:k) {
    s1 <- 0
    for (j in 1:ncol(Distance.centers)) {
      s1 <- s1 + Distance.centers[u, j]
    }
    s2 <- s2 + 1 / s1
  }
  Dis <- (Dmax / Dmin) * s2
  return(Dis)
}
