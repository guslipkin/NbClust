centers <- function(cl, x) {
  x <- as.matrix(x)

  centers <-
    cl |>
    unique() |>
    sort() |>
    lapply(\(i) {
      x[cl == i,] |>
        vector_to_matrix() |>
        apply(2, mean)
    }) |>
    unlist() |>
    matrix(ncol = length(cl), byrow = TRUE)

  return(centers)
}

vector_to_matrix <- function(x) {
  if (!is.matrix(x)) { x <- matrix(x, nrow = 1) }
  return(x)
}
