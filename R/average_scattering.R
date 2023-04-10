Average.scattering <- function (cl, x) {
  x <- as.matrix(x)
  n <- length(cl)
  k <- max(cl)
  centers.matrix <- centers(cl,x)

  cluster.size <-
    cl |>
    table() |>
    as.numeric()

  var.clusters <-
    cl |>
    unique() |>
    sort() |>
    lapply(\(i) {
      centers.matrix |>
        ncol() |>
        seq_len() |>
        sapply(\(j) { sum((x[cl == i,j] - centers.matrix[i,j])^2) }) /
          cluster.size[i]
    }) |>
    unlist() |>
    matrix(nrow = max(cl), byrow = TRUE)

  var.matrix <- apply(x, 2, \(x) { var(x)*(n-1)/n })

  Somme.var.clusters <-
    var.clusters |>
    apply(1, \(x) (x %*% x)) |>
    sqrt() |>
    sum()

  # Standard deviation
  stdev <- (1 / k) * sqrt(Somme.var.clusters)

  #Average scattering for clusters
  scat <-
    (1 / k) * (Somme.var.clusters / sqrt(var.matrix %*% var.matrix))

  scat <- list(
    "stdev" = stdev,
    "centers" = centers.matrix,
    "variance.intraclusters" = var.clusters,
    "scatt" = scat
  )
  return(scat)
}
