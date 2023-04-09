wss <- function(x) {
  x <- as.matrix(x)
  n <- length(x)
  centers <- matrix(nrow = 1, ncol = ncol(x))

  if (ncol(x) == 1) {	centers[1, ] <- mean(x) 	}
  if (is.null(dim(x))) {
    bb <- matrix(x,byrow=FALSE,nrow=1,ncol=ncol(x))
    centers[1, ] <- apply(bb, 2, mean)
  } else {
    centers[1, ] <- apply(x, 2, mean)
  }

  x.2 <- sweep(x,2,centers[1,],"-")
  withins <- sum(x.2^2)
  wss <- sum(withins)
  return(wss)
}

gss <- function(x, cl, d) {
  n <- length(cl)
  k <- max(cl)
  centers <- matrix(nrow = k, ncol = ncol(x))
  for (i in 1:k) {
    if (ncol(x) == 1) {
      centers[i,] <- mean(x[cl == i,])
    }
    if (is.null(dim(x[cl == i,]))) {
      bb <- matrix(x[cl == i,],
                   byrow = FALSE,
                   nrow = 1,
                   ncol = ncol(x))
      centers[i,] <- apply(bb, 2, mean)
    } else {
      centers[i,] <- apply(x[cl == i,], 2, mean)
    }
  }
  allmean <- apply(x, 2, mean)
  dmean <- sweep(x, 2, allmean, "-")
  allmeandist <- sum(dmean ^ 2)
  withins <- rep(0, k)
  x.2 <- (x - centers[cl,]) ^ 2
  for (i in 1:k) {
    withins[i] <- sum(x.2[cl == i,])
  }
  wgss <- sum(withins)
  bgss <- allmeandist - wgss

  results <- list(wgss = wgss,
                  bgss = bgss,
                  centers = centers)
  return(results)
}

vargss <- function(x, clsize, varwithins) {
  nvar <- dim(x)[2]
  n <- sum(clsize)
  k <- length(clsize)
  varallmean <- rep(0, nvar)
  varallmeandist <- rep(0, nvar)
  varwgss <- rep(0, nvar)
  for (l in 1:nvar)
    varallmean[l] <- mean(x[, l])
  vardmean <- sweep(x, 2, varallmean, "-")
  for (l in 1:nvar) {
    varallmeandist[l] <- sum((vardmean[, l]) ^ 2)
    varwgss[l] <- sum(varwithins[, l])
  }
  varbgss <- varallmeandist - varwgss
  vartss <- varbgss + varwgss
  zvargss <- list(vartss = vartss, varbgss = varbgss)
  return(zvargss)
}

varwithinss <- function(x, centers, cluster) {
  nrow <- dim(centers)[1]
  nvar <- dim(x)[2]
  varwithins <- matrix(0, nrow, nvar)
  x <- (x - centers[cluster,]) ^ 2
  for (l in 1:nvar) {
    for (k in 1:nrow) {
      varwithins[k, l] <- sum(x[cluster == k, l])
    }
  }
  return(varwithins)
}

biserial.cor <- function (x, y, use = c("all.obs", "complete.obs"), level = 1) {
  if (!is.numeric(x))
    stop("'x' must be a numeric variable.\n")
  y <- as.factor(y)
  if (length(levs <- levels(y)) > 2)
    stop("'y' must be a dichotomous variable.\n")
  if (length(x) != length(y))
    stop("'x' and 'y' do not have the same length")
  use <- match.arg(use)
  if (use == "complete.obs") {
    cc.ind <- complete.cases(x, y)
    x <- x[cc.ind]
    y <- y[cc.ind]
  }
  ind <- y == levs[level]
  diff.mu <- mean(x[ind]) - mean(x[!ind])
  prob <- mean(ind)
  diff.mu * sqrt(prob * (1 - prob)) / sd(x)
}
