indice.kl <- function (x, clall, d = NULL, nb1.cl1) {
  if (nb1.cl1 > 0) {
    KL <- NA
  }
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  m <- ncol(x)
  g <- k <- max(clall[, 2])
  KL <-
    abs((g - 1) ^ (2 / m) * gss(x, clall[, 1], d)$wgss - g ^ (2 / m) * gss(x, clall[, 2], d)$wgss) / abs((g) ^ (2 / m) * gss(x, clall[, 2], d)$wgss - (g + 1) ^ (2 / m) * gss(x, clall[, 3], d)$wgss)
  return(KL)
}

indice.ch <- function (x, cl, d = NULL, nb1.cl1){
  if (nb1.cl1 > 0) {
    CH <- NA
  }
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  n <- length(cl)
  k <- max(cl)
  CH <- (gss(x, cl, d)$bgss / (k - 1)) / (gss(x, cl, d)$wgss / (n - k))
  return(CH)
}

indice.hart <- function(x, clall, d = NULL, centrotypes = "centroids"){
  if (sum(c("centroids", "medoids") == centrotypes) == 0)
    stop("Wrong centrotypes argument")
  if ("medoids" == centrotypes && is.null(d))
    stop("For argument centrotypes = 'medoids' d cannot be null")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  n <- nrow(x)
  g <- max(clall[, 1])
  HART <-
    (gss(x, clall[, 2], d)$wgss / gss(x, clall[, 3], d)$wgss - 1) * (n - g - 1)
  return(HART)
}

indice.ball <- function(x, cl, d = NULL, centrotypes = "centroids"){
  wgssB <- gss(x, cl, d)$wgss
  qq <- max(cl)
  ball <- wgssB / qq
  return(ball)
}

indice.ratkowsky <- function(x, cl, d, centrotypes = "centroids"){
  qq <- max(cl)
  clsize <- table(cl)
  centers <- gss(x, cl, d)$centers
  varwithins <- varwithinss(x, centers, cl)
  zvargss <- vargss(x, clsize, varwithins)
  ratio <- mean(sqrt(zvargss$varbgss / zvargss$vartss))
  ratkowsky <- ratio / sqrt(qq)
  return(ratkowsky)
}
