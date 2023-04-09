Index.Hubert <- function(x, cl) {
  k <- max(cl)
  n <- dim(x)[1]
  y <- matrix(0, ncol = dim(x)[2], nrow = n)
  P <- as.matrix(md)
  meanP <- mean(P)
  variance.matrix <- numeric(0)
  md <- dist(x, method = "euclidean")
  for (j in 1:n)
    variance.matrix[j] = var(P[, j]) * (n - 1) / n
  varP <- sqrt(variance.matrix %*% variance.matrix)

  centers.clusters <- centers(cl, x)
  for (i in 1:n) {
    for (u in 1:k) {
      if (cl[i] == u)
        y[i,] <- centers.clusters[u,]
    }
  }

  Q <- as.matrix(dist(y, method = "euclidean"))
  meanQ <- mean(Q)
  for (j in 1:n)
    variance.matrix[j] = var(Q[, j]) * (n - 1) / n
  varQ <- sqrt(variance.matrix %*% variance.matrix)

  M <- n * (n - 1) / 2
  S <- 0
  n1 <- n - 1

  for (i in 1:n1) {
    j <- i + 1
    while (j <= n) {
      S <- S + (P[i, j] - meanP) * (Q[i, j] - meanQ)
      j <- j + 1
    }
  }
  gamma <- S / (M * varP * varQ)

  return(gamma)
}

Index.sPlussMoins <- function (cl1, md) {
  cn1 <- max(cl1)
  n1 <- length(cl1)
  dmat <- as.matrix(md)
  average.distance <-
    median.distance <-
    separation <-
    cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
  separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
  di <- list()
  for (u in 1:cn1) {
    cluster.size[u] <- sum(cl1 == u)
    du <- as.dist(dmat[cl1 == u, cl1 == u])
    within.dist1 <- c(within.dist1, du)
    average.distance[u] <- mean(du)
    median.distance[u] <- median(du)
    bv <- numeric(0)
    for (v in 1:cn1) {
      if (v != u) {
        suv <- dmat[cl1 == u, cl1 == v]
        bv <- c(bv, suv)
        if (u < v) {
          separation.matrix[u, v] <- separation.matrix[v, u] <- min(suv)
          between.dist1 <- c(between.dist1, suv)
        }
      }
    }
  }

  nwithin1 <- length(within.dist1)
  nbetween1 <- length(between.dist1)
  meanwithin1 <- mean(within.dist1)
  meanbetween1 <- mean(between.dist1)

  s.plus <- s.moins <- 0
  for (k in 1:nwithin1) {
    s.plus <-
      s.plus + (colSums(outer(between.dist1, within.dist1[k], ">")))
    s.moins <-
      s.moins + (colSums(outer(between.dist1, within.dist1[k], "<")))
  }

  Index.Gamma <- (s.plus - s.moins) / (s.plus + s.moins)
  Index.Gplus <- (2 * s.moins) / (n1 * (n1 - 1))
  t.tau  <- (nwithin1 * nbetween1) - (s.plus + s.moins)
  Index.Tau <-
    (s.plus - s.moins) / (((n1 * (n1 - 1) / 2 - t.tau) * (n1 * (n1 - 1) / 2)) ^
                            (1 / 2))

  results <- c(
    "gamma" = Index.Gamma,
    "gplus" = Index.Gplus,
    "tau" = Index.Tau
    )
  return(results)
}

Index.15and28  <- function (cl1, cl2, md) {
  cn1 <- max(cl1)
  n1 <- length(cl1)
  dmat <- as.matrix(md)
  average.distance <-
    median.distance <-
    separation <-
    cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
  separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
  di <- list()
  for (u in 1:cn1) {
    cluster.size[u] <- sum(cl1 == u)
    du <- as.dist(dmat[cl1 == u, cl1 == u])
    within.dist1 <- c(within.dist1, du)
    for (v in 1:cn1) {
      if (v != u) {
        suv <- dmat[cl1 == u, cl1 == v]
        if (u < v) {
          separation.matrix[u, v] <- separation.matrix[v, u] <- min(suv)
          between.dist1 <- c(between.dist1, suv)
        }
      }
    }
  }
  cn2 <- max(cl2)
  n2 <- length(cl2)
  dmat <- as.matrix(md)
  average.distance <-
    median.distance <-
    separation <-
    cluster.size <- within.dist2 <- between.dist2 <- numeric(0)
  separation.matrix <- matrix(0, ncol = cn2, nrow = cn2)
  di <- list()
  for (w in 1:cn2) {
    cluster.size[w] <- sum(cl2 == w)
    dw <- as.dist(dmat[cl2 == w, cl2 == w])
    within.dist2 <- c(within.dist2, dw)
    bx <- numeric(0)
    for (x in 1:cn2) {
      if (x != w) {
        swx <- dmat[cl2 == w, cl2 == x]
        bx <- c(bx, swx)
        if (w < x) {
          separation.matrix[w, x] <- separation.matrix[x, w] <- min(swx)
          between.dist2 <- c(between.dist2, swx)
        }
      }
    }
  }
  nwithin1 <- length(within.dist1)
  nbetween1 <- length(between.dist1)
  meanwithin1 <- mean(within.dist1)
  meanbetween1 <- mean(between.dist1)
  meanwithin2 <- mean(within.dist2)
  meanbetween2 <- mean(between.dist2)
  Index.15 <-
    (meanbetween2 - meanbetween1) / (meanwithin2 - meanwithin1)
  Index.28 <- (meanwithin1 / nwithin1) / (meanbetween1 / nbetween1)

  results <- c("frey" = Index.15, "mcclain" = Index.28)
  return(results)
}

Indices.WKWL <- function (x, cl1 = cl1, cl2 = cl2) {
  dim2 <- dim(x)[2]

  ncg1 <- 1
  ncg1max <- max(cl1)
  while ((sum(cl1 == ncg1) == sum(cl2 == ncg1)) && ncg1 <= ncg1max) {
    ncg1 <- ncg1 + 1
  }
  g1 <- ncg1

  ncg2 <- max(cl2)
  nc2g2 <- ncg2 - 1
  while ((sum(cl1 == nc2g2) == sum(cl2 == ncg2)) && nc2g2 >= 1) {
    ncg2 <- ncg2 - 1
    nc2g2 <- nc2g2 - 1
  }
  g2 <- ncg2

  NK <- sum(cl2 == g1)
  WK.x <- x[cl2 == g1, ]
  WK <- wss(x = WK.x)

  NL <- sum(cl2 == g2)
  WL.x <- x[cl2 == g2, ]
  WL <- wss(x = WL.x)

  NM <- sum(cl1 == g1)
  WM.x <- x[cl1 == g1, ]
  WM <- wss(x = WM.x)

  duda <- (WK + WL) / WM

  BKL <- WM - WK - WL
  pseudot2 <- BKL / ((WK + WL) / (NK + NL - 2))

  beale <- (BKL / (WK + WL)) / (((NM - 1) / (NM - 2)) * (2 ^ (2 / dim2) -
                                                           1))

  results <- c(
      "duda" = duda,
      "pseudot2" = pseudot2,
      "NM" = NM,
      "NK" = NK,
      "NL" = NL,
      "beale" = beale
    )
  return(results)
}

Indices.WBT <- function(x, cl, P, s, vv) {
  n <- dim(x)[1]
  pp <- dim(x)[2]
  qq <- max(cl)
  z <- matrix(0, ncol = qq, nrow = n)
  clX <- as.matrix(cl)

  for (i in 1:n)
    for (j in 1:qq) {
      z[i, j] == 0
      if (clX[i, 1] == j) {
        z[i, j] = 1
      }
    }

  xbar <- solve(t(z) %*% z) %*% t(z) %*% x
  B <- t(xbar) %*% t(z) %*% z %*% xbar
  W <- P - B
  marriot <- (qq ^ 2) * det(W)
  trcovw <- sum(diag(cov(W)))
  tracew <- sum(diag(W))
  if (det(W) != 0)
    scott <- n * log(det(P) / det(W))
  else {
    cat("Error: division by zero!")
  }
  friedman <- sum(diag(solve(W) * B))
  rubin <- sum(diag(P)) / sum(diag(W))


  R2 <- 1 - sum(diag(W)) / sum(diag(P))
  v1 <- 1
  u <- rep(0, pp)
  c <- (vv / (qq)) ^ (1 / pp)
  u <- s / c
  k1 <- sum((u >= 1) == TRUE)
  p1 <- min(k1, qq - 1)
  if (all(p1 > 0, p1 < pp)) {
    for (i in 1:p1)
      v1 <- v1 * s[i]
    c <- (v1 / (qq)) ^ (1 / p1)
    u <- s / c
    b1 <- sum(1 / (n + u[1:p1]))
    b2 <- sum(u[p1 + 1:pp] ^ 2 / (n + u[p1 + 1:pp]), na.rm = TRUE)
    E_R2 <- 1 - ((b1 + b2) / sum(u ^ 2)) * ((n - qq) ^ 2 / n) * (1 + 4 /
                                                                   n)
    ccc <- log((1 - E_R2) / (1 - R2)) * (sqrt(n * p1 / 2) / ((0.001 + E_R2) ^
                                                               1.2))
  } else {
    b1 <- sum(1 / (n + u))
    E_R2 <- 1 - (b1 / sum(u ^ 2)) * ((n - qq) ^ 2 / n) * (1 + 4 / n)
    ccc <- log((1 - E_R2) / (1 - R2)) * (sqrt(n * pp / 2) / ((0.001 + E_R2) ^
                                                               1.2))
  }
  results <- c(
      "ccc" = ccc,
      "scott" = scott,
      "marriot" = marriot,
      "trcovw" = trcovw,
      "tracew" = tracew,
      "friedman" = friedman,
      "rubin" = rubin
  )
  return(results)
}

Indices.Traces <- function(data, d, clall, index = "all") {
  x <- data
  cl0 <- clall[, 1]
  cl1 <- clall[, 2]
  cl2 <- clall[, 3]
  clall <- clall
  nb.cl0 <- table(cl0)
  nb.cl1 <- table(cl1)
  nb.cl2 <- table(cl2)
  nb1.cl0 <- sum(nb.cl0 == 1)
  nb1.cl1 <- sum(nb.cl1 == 1)
  nb1.cl2 <- sum(nb.cl2 == 1)

  indice <-
    pmatch(index, c("kl", "ch", "hartigan", "ratkowsky", "ball", "all"))
  if (is.na(indice))
    stop("invalid clustering index")
  if (indice == -1)
    stop("ambiguous index")
  vecallindex <- numeric(5)
  if (any(indice == 1) || (indice == 6))
    vecallindex[1] <- indice.kl(x, clall, d, nb1.cl1)
  if (any(indice == 2) || (indice == 6))
    vecallindex[2] <- indice.ch(x, cl = clall[, 2], d, nb1.cl1)
  if (any(indice == 3) || (indice == 6))
    vecallindex[3] <- indice.hart(x, clall, d)
  if (any(indice == 4) || (indice == 6))
    vecallindex[4] <- indice.ratkowsky(x, cl = cl1, d)
  if (any(indice == 5) || (indice == 6))
    vecallindex[5] <- indice.ball(x, cl = cl1, d)
  names(vecallindex) <- c("kl", "ch", "hartigan", "ratkowsky", "ball")
  if (indice < 6)
    vecallindex <- vecallindex[indice]
  return(vecallindex)
}

Indice.ptbiserial <- function (x, md, cl1) {
  nn <- dim(x)[1]
  pp <- dim(x)[2]

  md2 <- as.matrix(md)
  m01 <- array(NA, c(nn, nn))
  nbr <- (nn * (nn - 1)) / 2
  pb <- array(0, c(nbr, 2))

  m3 <- 1
  for (m1 in 2:nn) {
    m12 <- m1 - 1
    for (m2 in 1:m12) {
      if (cl1[m1] == cl1[m2])
        m01[m1, m2] <- 0
      if (cl1[m1] != cl1[m2])
        m01[m1, m2] <- 1
      pb[m3, 1] <- m01[m1, m2]
      pb[m3, 2] <- md2[m1, m2]
      m3 <- m3 + 1
    }
  }

  y <- pb[, 1]
  x <- pb[, 2]

  ptbiserial <- biserial.cor(x = pb[, 2], y = pb[, 1], level = 2)
  return(ptbiserial)
}

Indice.cindex <- function (d, cl) {
  d <- data.matrix(d)
  DU <- 0
  r <- 0
  v_max <- array(1, max(cl))
  v_min <- array(1, max(cl))
  for (i in 1:max(cl)) {
    n <- sum(cl == i)
    if (n > 1) {
      t <- d[cl == i, cl == i]
      DU = DU + sum(t) / 2
      v_max[i] = max(t)
      if (sum(t == 0) == n)
        v_min[i] <- min(t[t != 0])
      else
        v_min[i] <- 0
      r <- r + n * (n - 1) / 2
    }
  }
  Dmin <- min(v_min)
  Dmax <- max(v_max)
  if (Dmin == Dmax)
    result <- NA
  else
    result <- (DU - r * Dmin) / (Dmax * r - Dmin * r)
  result
}

Indice.DB <- function (x, cl, d = NULL, centrotypes = "centroids", p = 2, q = 2) {
  if (sum(c("centroids") == centrotypes) == 0)
    stop("Wrong centrotypes argument")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  if (is.null(dim(x))) {
    dim(x) <- c(length(x), 1)
  }
  x <- as.matrix(x)
  n <- length(cl)
  k <- max(cl)
  dAm <- d
  centers <- matrix(nrow = k, ncol = ncol(x))
  if (centrotypes == "centroids") {
    for (i in 1:k) {
      for (j in 1:ncol(x)) {
        centers[i, j] <- mean(x[cl == i, j])
      }
    }
  }
  else {
    stop("wrong centrotypes argument")
  }
  S <- rep(0, k)
  for (i in 1:k) {
    ind <- (cl == i)
    if (sum(ind) > 1) {
      centerI <- centers[i,]
      centerI <- rep(centerI, sum(ind))
      centerI <- matrix(centerI, nrow = sum(ind), ncol = ncol(x), byrow = TRUE)
      S[i] <- mean(sqrt(apply((x[ind,] - centerI) ^ 2, 1, sum)) ^ q) ^ (1 / q)
    }
    else
      S[i] <- 0
  }
  M <- as.matrix(dist(centers, p = p))
  R <- array(Inf, c(k, k))
  r = rep(0, k)
  for (i in 1:k) {
    for (j in 1:k) {
      R[i, j] = (S[i] + S[j]) / M[i, j]
    }
    r[i] = max(R[i,][is.finite(R[i,])])
  }
  DB = mean(r[is.finite(r)])
  resul <-
    list(
      DB = DB,
      r = r,
      R = R,
      d = M,
      S = S,
      centers = centers
    )
  resul
}

Indice.S <- function (d, cl) {
  d <- as.matrix(d)
  Si <- 0
  for (k in 1:max(cl)) {
    if ((sum(cl == k)) <= 1)
      Sil <- 1
    else {
      Sil <- 0
      for (i in 1:length(cl)) {
        if (cl[i] == k) {
          ai <- sum(d[i, cl == k]) / (sum(cl == k) - 1)
          dips <- NULL
          for (j in 1:max(cl))
            if (cl[i] != j)
              if (sum(cl == j) != 1)
                dips <- cbind(dips, c((sum(d[i, cl == j])) / (sum(cl == j))))
          else
            dips <- cbind(dips, c((sum(d[i, cl == j]))))
          bi <- min(dips)
          Sil <- Sil + (bi - ai) / max(c(ai, bi))
        }
      }
    }
    Si <- Si + Sil
  }
  Si / length(cl)
}

Indice.Gap <- function (x, clall, reference.distribution = "unif", B = 10, method = "ward.D2", d = NULL, centrotypes = "centroids") {
  if (sum(c("centroids", "medoids") == centrotypes) == 0)
    stop("Wrong centrotypes argument")
  if ("medoids" == centrotypes && is.null(d))
    stop("For argument centrotypes = 'medoids' d can not be null")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  X <- as.matrix(x)
  gap1 <- GAP(X, clall[, 1], reference.distribution, B, method, d, centrotypes)
  gap <- gap1$Sgap
  gap2 <- GAP(X, clall[, 2], reference.distribution, B, method, d, centrotypes)
  diffu <- gap - (gap2$Sgap - gap2$Sdgap)
  resul <- list(gap = gap, diffu = diffu)
  resul
}

Index.sdindex <- function(x, clmax, cl) {
  x <- as.matrix(x)
  Alpha <- Dis(clmax, x)
  Scatt <- Average.scattering(cl, x)$scatt
  Dis0 <- Dis(cl, x)
  SD.indice <- Alpha * Scatt + Dis0
  return(SD.indice)
}

Index.SDbw <- function(x, cl) {
  x <- as.matrix(x)
  Scatt <- Average.scattering(cl, x)$scatt
  Dens.bw <- density.bw(cl, x)
  SDbw <- Scatt + Dens.bw
  return(SDbw)
}

Index.Dindex <- function(cl, x) {
  x <- as.matrix(x)
  distance <- density.clusters(cl, x)$distance
  n <- length(distance)
  S <- 0
  for (i in 1:n)
    S <- S + distance[i]
  inertieIntra <- S / n
  return(inertieIntra)
}

Index.dunn <- function(md, clusters, Data=NULL, method="euclidean") {
  distance <- as.matrix(md)
  nc <- max(clusters)
  interClust <- matrix(NA, nc, nc)
  intraClust <- rep(NA, nc)

  for (i in 1:nc) {
    c1 <- which(clusters == i)
    for (j in i:nc) {
      if (j == i)
        intraClust[i] <- max(distance[c1, c1])
      if (j > i) {
        c2 <- which(clusters == j)
        interClust[i, j] <- min(distance[c1, c2])
      }
    }
  }
  dunn <- min(interClust, na.rm = TRUE) / max(intraClust)
  return(dunn)
}
