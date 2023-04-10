simgap <- function(Xvec) {
  ma <- max(Xvec)
  mi <- min(Xvec)
  set.seed(1)
  Xout <- runif(length(Xvec), min = mi, max = ma)
  return(Xout)
}

GAP <- function(X, cl, method) {
  B <- 10
  set.seed(1)
  if (is.null(dim(x))) {
    dim(x) <- c(length(x), 1)
  }
  ClassNr <- max(cl)
  Wk0 <- 0
  WkB <- matrix(0, 1, B)
  for (bb in (1:B)) {
    Xnew <- apply(X, 2, simgap)
    if (bb == 1) {
      pp <- cl
      if (ClassNr == length(cl))
        pp2 <- 1:ClassNr
      else if (method == "k-means") {
        set.seed(1)
        pp2 <- kmeans(Xnew, ClassNr, 100)$cluster
      }
      else if (method == "single" || method == "complete" ||
               method == "average" || method == "ward.D2" ||
               method == "mcquitty" || method == "median" ||
               method == "centroid" || method == "ward.D")
        pp2 <- cutree(hclust(dist(Xnew), method = method), ClassNr)
      else
        stop("Wrong clustering method")
      if (ClassNr > 1) {
        for (zz in (1:ClassNr)) {
          Xuse <- X[pp == zz, ]
          Wk0 <- Wk0 + sum(diag(var(Xuse))) * (length(pp[pp == zz]) - 1) / (dim(X)[1] - ClassNr)
          Xuse2 <- Xnew[pp2 == zz, ]
          WkB[1, bb] <- WkB[1, bb] + sum(diag(var(Xuse2))) * (length(pp2[pp2 == zz]) - 1) / (dim(X)[1] - ClassNr)
        }
      }
      if (ClassNr == 1) {
        Wk0 <- sum(diag(var(X)))
        WkB[1, bb] <- sum(diag(var(Xnew)))
      }
    }
    if (bb > 1) {
      if (ClassNr == length(cl))
        pp2 <- 1:ClassNr
      else if (method == "k-means") {
        set.seed(1)
        pp2 <- kmeans(Xnew, ClassNr, 100)$cluster
      }
      else if (method == "single" || method == "complete" ||
               method == "average" || method == "ward.D2" ||
               method == "mcquitty" || method == "median" ||
               method == "centroid" || method == "ward.D")
        pp2 <- cutree(hclust(dist(Xnew), method = method), ClassNr)
      else
        stop("Wrong clustering method")
      if (ClassNr > 1) {
        for (zz in (1:ClassNr)) {
          Xuse2 <- Xnew[pp2 == zz, ]
          WkB[1, bb] <- WkB[1, bb] + sum(diag(var(Xuse2))) *
            length(pp2[pp2 == zz]) / (dim(X)[1] - ClassNr)
        }
      }
      if (ClassNr == 1) {
        WkB[1, bb] <- sum(diag(var(Xnew)))
      }
    }
  }
  Sgap <- mean(log(WkB[1, ])) - log(Wk0)
  Sdgap <- sqrt(1 + 1 / B) * sqrt(var(log(WkB[1, ]))) * sqrt((B - 1) / B)
  resul <- list(Sgap = Sgap, Sdgap = Sdgap)
  resul
}
