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
    matrix(nrow = max(cl), byrow = TRUE)

  return(centers)
}

vector_to_matrix <- function(x) {
  if (!is.matrix(x)) { x <- matrix(x, nrow = 1) }
  return(x)
}

vector_to_col_matrix <- function(x, name = NULL) {
  y <- x
  if (!is.matrix(x)) {
    y <- matrix(x, ncol = 1)
    if (!is.null(name)) {
      colnames(y) <- name
      rownames(y) <- names(x)
    }
  }
  return(y)
}
