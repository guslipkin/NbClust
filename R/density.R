density.clusters <- function(cl, x) {
  x <- as.matrix(x)
  k <- max(cl)
  n <- length(cl)

  distance <- matrix(0, ncol = 1, nrow = n)
  density <-  matrix(0, ncol = 1, nrow = k)

  average_scattering <- Average.scattering(cl, x)
  centers.matrix <- average_scattering$centers
  stdev <- average_scattering$stdev

  distance <-
    cl |>
    unique() |>
    lapply(\(i) {
      centers.matrix |>
        ncol() |>
        seq_len() |>
        sapply(\(j) { (x[cl == i,j] - centers.matrix[i,j])^2 }) |>
        vector_to_matrix() |>
        apply(1, sum) |>
        sqrt()
    }) |>
    unlist() |>
    matrix(ncol = 1) |>
    cbind(matrix(cl, ncol = 1)) |>
    data.frame()
  distance <-
    distance[order(distance$X2),1] |>
    as.matrix()

  density <-
    as.numeric(distance <= stdev)[unique(cl)] |>
    as.matrix(ncol = 1)

  dens <- list(
    "distance" = distance,
    "density" = density
  )
  return(dens)
}


density.bw<-function(cl, x) {
  x <- as.matrix(x)
  k <- max(cl)
  n <- length(cl)

  average_scattering <- Average.scattering(cl, x)
  centers.matrix <- average_scattering$centers
  stdev <- average_scattering$stdev
  u <- 1

  density.bw <-
    expand.grid("u" = 1:k, "v" = 1:k) |>
    apply(1, \(d) {
      u <- d["u"]; v <- d["v"];
      if (u == v) { return(0) }
      moy <- (centers.matrix[u,] + centers.matrix[v,]) / 2
      dist <-
        moy |>
        seq_along() |>
        sapply(\(j) { (x[cl == v | cl == u,j] - moy[j])^2 }) |>
        matrix(ncol = length(moy)) |>
        apply(1, sum) |>
        sqrt()
      sum(dist <= stdev)
    }) |>
    matrix(ncol = k)

  density.clust <- density.clusters(cl,x)$density

  density.bw <-
    expand.grid("u" = 1:k, "v" = 1:k) |>
    apply(1, \(d) {
      u <- d["u"]; v <- d["v"];
      uv_max <- max(density.clust[u], density.clust[v])
      d <- if (uv_max != 0) { density.bw[u,v] / uv_max } else { 0 }
      d / (k * (k - 1))
    }) |>
    sum()

  return(density.bw)
}
