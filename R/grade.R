grade <- function(res, resCritical, master_index, row_names, alphaBeale, DiffLev) {
  row_names <- as.numeric(row_names)
  mapply(\(x, name) {
    if (name %in% master_index[c(3:6, 15, 19, 21, 25, 27:28)]) {
      x <- grade_max(x, row_names)
    } else if (name %in% master_index[c(13:14, 24, 26, 30, 32)]) {
      x <- grade_min(x, row_names)
    } else if (name %in% master_index[c(16:17, 22:23)]) {
      row_names <- row_names - 1
      if (name == "gap") {
        y <- which(resCritical[row_names, "CritValue_Gap"] >= 0)[1]
      } else if (name == "duda") {
        y <- which(x >= resCritical[row_names, "CritValue_Duda"])[1]
      } else if (name == "pseudot2") {
        y <- which(x <= resCritical[row_names, "CritValue_PseudoT2"])[1]
      } else if (name == "frey") {
        y <- which(x < 1)[1]
      }
      x <- grade_inequality(x, y)
    } else if (name == "beale") {
      y <- sort(x)
      y <- y[which(y >= alphaBeale)[1]]
      x <- matrix(c(row_names[order(x)[1]], y), ncol = 1)
    } else if (name %in% master_index[c(7:11, 20)]) {
      x <- grade_max(DiffLev[,name], row_names)
    } else if (name == "rubin") {
      x <- grade_min(DiffLev[,name], row_names)
    } else if (name == "hubert") {
      x <- grade_hubert(x, DiffLev[,"hubert"], row_names)
    } else if (name == "dindex") {
      x <- grade_dindex(x, DiffLev[,"dindex"], row_names)
    }
    colnames(x) <- name
    return(x)
  }, data.frame(res), colnames(res))
}

grade_max <- function(x, row_names) {
  clusters <- row_names[which.max(x)]
  index <- max(x, na.rm = TRUE)
  matrix(c(clusters, index), ncol = 1)
}

grade_min <- function(x, row_names) {
  clusters <- row_names[which.min(x)]
  index <- min(x, na.rm = TRUE)
  matrix(c(clusters, index), ncol = 1)
}

grade_inequality <- function(x, clusters) {
  index <- x[clusters]
  clusters <-
    clusters |>
    names() |>
    as.numeric()
  if (length(clusters) == 0) { clusters <- NA_real_ }
  if (is.na(index)) { index <- NA_real_ }
  matrix(c(clusters, index), ncol = 1)
}

grade_hubert <- function(x, y, row_names) {
  par(mfrow = c(1, 2))
  plot(
    x_axis,
    x,
    tck = 0,
    type = "b",
    col = "red",
    xlab = expression(paste("Number of clusters ")),
    ylab = expression(paste("Hubert Statistic values"))
  )
  plot(
    row_names,
    y,
    tck = 0,
    type = "b",
    col = "blue",
    xlab = expression(paste("Number of clusters ")),
    ylab = expression(paste("Hubert statistic second differences"))
  )
  cat("*** : The Hubert index is a graphical method of determining the number of clusters.
                In the plot of Hubert index, we seek a significant knee that corresponds to a
                significant increase of the value of the measure i.e the significant peak in Hubert
                index second differences plot.\n\n")
  matrix(c(0, 0), ncol = 1)
}

grade_dindex <- function(x, y, row_names) {
  par(mfrow = c(1, 2))
  plot(
    x_axis,
    x,
    tck = 0,
    type = "b",
    col = "red",
    xlab = expression(paste("Number of clusters ")),
    ylab = expression(paste("Dindex Values"))
  )
  plot(
    row_names,
    y,
    tck = 0,
    type = "b",
    col = "blue",
    xlab = expression(paste("Number of clusters ")),
    ylab = expression(paste("Second differences Dindex Values"))
  )
  cat("*** : The D index is a graphical method of determining the number of clusters.
                In the plot of D index, we seek a significant knee (the significant peak in Dindex
                second differences plot) that corresponds to a significant increase of the value of
                the measure.\n\n")
  matrix(c(0, 0), ncol = 1)
}

create_difflev <- function(res, row_names) {
  DiffLev <- matrix(0, nrow = length(row_names), ncol = 12)
  DiffLev[,1] <- as.numeric(row_names)

  # top row
  DiffLev[1,c(2:10, 12)] <- NA_real_

  # middle rows
  for (n in 2:(nrow(DiffLev) - 1)) {
    cols <- c(3,5,7,9,18,27)
    DiffLev[n, c(2,3,5,7,9,10)] <- abs(res[n, cols] - res[n - 1, cols])

    cols <- c(6,8,10,29)
    DiffLev[n, c(4,6,8,12)] <-
      ((res[n+1, cols] - res[n, cols]) - (res[n, cols] - res[n-1, cols]))
  }

  # bottom row
  max_row <- nrow(DiffLev)
  DiffLev[max_row, c(4, 6, 8, 10, 12)] <- NA_real_

  DiffLev[max_row, c(2:3,5,7,9)] <-
    abs(res[max_row, c(3,5,7,9,18)] - res[max_row - 1, c(3,5,7,9,18)])

  colnames(DiffLev) <- c(
    "clusters", "hartigan", "scott", "marriot", "trcovw", "tracew",
    "friedman", "rubin", "ball", "hubert", "skipped", "dindex"
  )
  return(DiffLev)
}
