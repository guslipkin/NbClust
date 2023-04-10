NbClust <- function(
    data = NULL,
    diss = NULL,
    distance = c(
      "diss", "euclidean", "maximum", "manhattan", "canberra", "binary",
      "minkowski"
    ),
    min.nc = 2,
    max.nc = 15,
    method = c(
      "ward.D2", "single", "complete", "average", "mcquitty", "median",
      "centroid", "kmeans","ward.D"
    ),
    index = c("all", "alllong",
              "kl", "ch", "hartigan", "ccc", "scott",
              "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex",
              "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky",
              "ball", "ptbiserial", "gap",  "frey",  "mcclain",   "gamma",
              "gplus",  "tau",  "dunn", "hubert",  "sdindex",  "dindex",
              "sdbw"
    ),
    alphaBeale = 0.1
    ) {

  master_distance = c(
    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  )
  master_method <- c(
    "ward.D2", "single", "complete", "average", "mcquitty", "median",
    "centroid", "kmeans","ward.D"
  )
  master_index <- c(
    "all", "alllong", "kl", "ch", "hartigan", "ccc", "scott", "marriot",
    "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette",
    "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap",
    "frey",  "mcclain",   "gamma", "gplus",  "tau",  "dunn", "hubert",
    "sdindex",  "dindex", "sdbw"
    )

  x <- 0
  min_nc <- min.nc
  max_nc <- max.nc

  distance <- rlang::arg_match0(distance, master_distance)
  method <- rlang::arg_match0(method, master_method)
  indice <- rlang::arg_match0(index, master_index)

  if (
    indice %in% master_index[c(1:2, 5, 7:13, 20, 29, 31)] &
    max.nc - min.nc < 2
    ) {
    stop("The difference between the minimum and the maximum number of clusters must be at least equal to 2")
  } else if (method == "kmeans" & is.null(data)) {
    stop("Method = kmeans, data matrix is needed")
  } else if (
    indice %in% master_index[c(1:12, 14, 16:22, 23:27, 29:32)] &
    is.null(data)
    ) {
    stop("Data matrix is needed. Only frey, mcclain, cindex, sihouette and dunn can be computed.")
  } else if (is.null(data) & is.null(diss)) {
    stop("data matrix and dissimilarity matrix are both null")
  } else if (is.null(data) & !is.null(diss)) {
    message("Only frey, mcclain, cindex, sihouette and dunn can be computed. To compute the other indices, data matrix is needed")
  }

  jeu1 <- as.matrix(data)
  numberObsBefore <- nrow(jeu1)
  jeu <- na.omit(jeu1)
  nn <- numberObsAfter <- nrow(jeu)
  pp <- ncol(jeu)
  TT <- t(jeu) %*% jeu
  eigenValues <- eigen(TT / (nn - 1))$value
  sizeEigenTT <- length(eigenValues)

  if (indice %in% master_index[c(1:2, 6:12)]) {
    if (!all(eigenValues >= 0)) {
      stop("The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated.")
    }
    s1 <- sqrt(eigenValues)
    ss <- rep(1, sizeEigenTT)
    ss[s1 != 0] <- s1[s1 != 0]
    vv <- prod(ss)
  }

  #### Distances ####

  stopifnot(
    !is.null(diss) && distance == "diss"
    || is.null(diss) && distance != "diss"
  )
  md <- if (!is.null(diss)) diss else dist(jeu, method = distance)

  #### Methods ####

  res <-
    matrix(0, max_nc - min_nc + 1, 30) |>
    `rownames<-`(min_nc:max_nc) |>
    `colnames<-`(master_index[-c(1:2)])
  x_axis <- min_nc:max_nc
  resCritical <-
    matrix(0, max_nc - min_nc + 1, 4) |>
    `rownames<-`(min_nc:max_nc) |>
    `colnames<-`(c(
      "CritValue_Duda", "CritValue_PseudoT2",
      "Fvalue_Beale", "CritValue_Gap"
    ))

  hc <- hclust(md, method = method)

  #### Clustering ####

  for (nc in min_nc:max_nc) {

    if (method %in% master_method[c(1:7, 9)]) {
      cl1 <- cutree(hc, k = nc)
      cl2 <- cutree(hc, k = nc + 1)
      clall <- cbind(cl1, cl2)
      clmax <- cutree(hc, k = max_nc)

      if (nc >= 2) {
        cl0 <- cutree(hc, k = nc - 1)
        clall1 <- cbind(cl0, cl1, cl2)
      } else if (nc == 1) {
        cl0 <- rep(NA, nn)
        clall1 <- cbind(cl0, cl1, cl2)
      }
    } else if (method == "kmeans") {
      set.seed(1)
      cl2 <- kmeans(jeu, nc + 1)$cluster
      set.seed(1)
      clmax <- kmeans(jeu, max_nc)$cluster
      if (nc > 2) {
        set.seed(1)
        cl1 <- kmeans(jeu, nc)$cluster
        clall <- cbind(cl1, cl2)
        set.seed(1)
        cl0 <- kmeans(jeu, nc - 1)$cluster
        clall1 <- cbind(cl0, cl1, cl2)
      } else if (nc == 2) {
        set.seed(1)
        cl1 <- kmeans(jeu, nc)$cluster
        clall <- cbind(cl1, cl2)
        cl0 <- rep(1, nn)
        clall1 <- cbind(cl0, cl1, cl2)
      } else if (nc == 1) {
        stop("Number of clusters must be higher than 2")
      }
    }

    # table uses the cross-classifying factors to build a contingency
    #   table of the counts at each combination of factor levels.
    j <- table(cl1)
    s <- sum(j == 1)
    j2 <- table(cl2)
    s2 <- sum(j2 == 1)

    if (indice %in% master_index[c(1:2, 5)]) {
      traces <- Indices.Traces(jeu, md, clall1, index = indice)
      res[nc - min_nc + 1, names(traces)] <- traces
    }

    if (indice %in% master_index[c(1:2, 6:12)]) {
      wbt <- Indices.WBT(x = jeu, cl = cl1, P = TT,s = ss,vv = vv)
      res[nc - min_nc + 1, names(wbt)] <- wbt
    }

    if (indice %in% master_index[c(1:2, 16:18)]) {
      wkwl <- Indices.WKWL(x = jeu, cl1 = cl1, cl2 = cl2)
      res[nc - min_nc + 1, c("duda", "pseudot2", "beale")] <-
        wkwl[c("duda", "pseudot2", "beale")]

      zz <- 3.20 # Best standard score in Milligan and Cooper 1985
      zzz <-
        zz * sqrt(2 * (1 - 8 / ((pi ^ 2) * pp)) / (wkwl["NM"] * pp))

      critValue <- 1 - (2 / (pi * pp)) - zzz
      resCritical[nc - min_nc + 1, "CritValue_Duda"] <- critValue
      resCritical[nc - min_nc + 1, "CritValue_PseudoT2"] <-
        ((1 - critValue) / critValue) * (wkwl["NK"] + wkwl["NL"] - 2)
      resCritical[nc - min_nc + 1, "Fvalue_Beale"] <-
        1 - pf(wkwl["beale"], pp, (wkwl["NM"] - 2) * pp)
    }

    if (indice %in% master_index[c(1:2, 21)]) {
      ptbiserial <- Indice.ptbiserial(x = jeu, md = md, cl1 = cl1)
      res[nc - min_nc + 1, "ptbiserial"] <- ptbiserial
    }

    if (indice %in% master_index[c(1:2, 22)]) {
      resultSGAP <- Indice.Gap(jeu, clall, method)
      res[nc - min_nc + 1, "gap"] <- resultSGAP["gap"]
      resCritical[nc - min_nc + 1, "CritValue_Gap"] <-
        resultSGAP["diffu"]
    }

    if (nc == 1) {
      res[nc - min_nc + 1, c(1:2, 11:13, 17, 21:30)] <- NA_real_
    } else {
      if (indice %in% master_index[c(1:4)]) {
        col_pos <- which(c("kl", "ch") == indice)
        traces <- Indices.Traces(jeu, md, clall1, index = indice)
        res[nc - min_nc + 1, col_pos] <- traces
      }

      if (indice %in% master_index[c(1:2, 13)]) {
        cindex <- Indice.cindex(d = md, cl = cl1)
        res[nc - min_nc + 1, "cindex"] <- cindex
      }

      if (indice %in% master_index[c(1:2, 14)]) {
        db <-  Indice.DB(
          x = jeu,
          cl = cl1,
          d = NULL,
          centrotypes = "centroids",
          p = 2,
          q = 2
        )
        res[nc - min_nc + 1, "db"] <- db$DB
      }

      if (indice %in% master_index[c(1:2, 15)]) {
        silhouette <- Indice.S(d = md, cl = cl1)
        res[nc - min_nc + 1, "silhouette"] <- silhouette
      }

      if (indice %in% master_index[c(1:2, 23:24)]) {
        fifteen_twentyeight <-
          Index.15and28(cl1 = cl1, cl2 = cl2, md = md)
        res[nc - min_nc + 1, names(fifteen_twentyeight)] <-
          fifteen_twentyeight
      }

      if (indice %in% master_index[c(1:2, 25:27)]) {
        splussmoins <- Index.sPlussMoins(cl1 = cl1, md = md)
        res[nc - min_nc + 1, names(splussmoins)] <- splussmoins
      }

      if (indice %in% master_index[c(1:2, 28)]) {
        dunn <- Index.dunn(md, cl1, Data = jeu, method = NULL)
        res[nc - min_nc + 1, "dunn"] <- dunn
      }

      if (indice %in% master_index[c(1:2, 29)]) {
        res[nc - min_nc + 1, "hubert"] <- Index.Hubert(jeu, cl1)
      }

      if (indice %in% master_index[c(1:2, 30)]) {
        res[nc - min_nc + 1, "sdindex"] <-
          Index.sdindex(jeu, clmax, cl1)
      }

      if (indice %in% master_index[c(1:2, 31)]) {
        res[nc - min_nc + 1, "dindex"] <- Index.Dindex(cl1, jeu)
      }

      if (indice %in% master_index[c(1:2, 32)]) {
        res[nc - min_nc + 1, "sdbw"] <- Index.SDbw(jeu, cl1)
      }
    }
}

#### Grading ####

  if (indice %in% master_index[c(1:2, 7:11, 20)]) {
    DiffLev <- create_difflev(res, rownames(res))
  } else {
    DiffLev <- NULL
  }

  col_res <- if (grepl("all", indice)) master_index[-c(1:2)] else indice
  res[,col_res] |>
    vector_to_col_matrix(col_res) |>
    grade(resCritical, master_index, rownames(res), alphaBeale, DiffLev)

######################################################################################################################
########################                Displaying results             #########################################
######################################################################################################################

 if (indice < 31)
 {
     res <- res[,c(indice)]

     if (indice == 14) { resCritical <- resCritical[,1]  }
     if (indice == 15) { resCritical <- resCritical[,2] }
     if (indice == 16) { resCritical <- resCritical[,3] }
     if (indice == 20) { resCritical <- resCritical[,4] }

 }

 if (indice == 31)
  {
      res <- res[,c(1:19,21:22,26:30)]
		  resCritical <- resCritical[,c(1:3)]

  }



 if (any(indice == 20) || (indice == 23) || (indice == 24) || (indice == 25) || (indice == 32))
 {

  results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan, indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
		nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW, nc.TraceW, indice.TraceW, nc.Friedman,
		indice.Friedman, nc.Rubin, indice.Rubin, nc.cindex, indice.cindex, nc.DB, indice.DB, nc.Silhouette,
		indice.Silhouette, nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo, nc.Beale, indice.Beale, nc.Ratkowsky,
		indice.Ratkowsky, nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial, nc.Gap, indice.Gap,
		nc.Frey, indice.Frey, nc.McClain, indice.McClain, nc.Gamma, indice.Gamma, nc.Gplus, indice.Gplus,
		nc.Tau, indice.Tau, nc.Dunn, indice.Dunn, nc.Hubert, indice.Hubert, nc.sdindex, indice.sdindex, nc.Dindex, indice.Dindex, nc.SDbw, indice.SDbw)
  results1<-matrix(c(results),nrow=2,ncol=30)
  resultats <- matrix(c(results),nrow=2,ncol=30,dimnames = list(c("Number_clusters","Value_Index"),
		     c("KL","CH","Hartigan","CCC", "Scott", "Marriot", "TrCovW",
		       "TraceW", "Friedman", "Rubin", "Cindex", "DB", "Silhouette",
			"Duda","PseudoT2", "Beale", "Ratkowsky", "Ball", "PtBiserial",
			"Gap", "Frey", "McClain", "Gamma", "Gplus", "Tau", "Dunn",
      "Hubert", "SDindex", "Dindex", "SDbw")))
   }
 else
  {

    results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan, indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
		nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW, nc.TraceW, indice.TraceW, nc.Friedman, indice.Friedman,
    nc.Rubin, indice.Rubin, nc.cindex, indice.cindex, nc.DB, indice.DB, nc.Silhouette, indice.Silhouette,
    nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo, nc.Beale, indice.Beale, nc.Ratkowsky, indice.Ratkowsky,
    nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial, nc.Frey, indice.Frey, nc.McClain, indice.McClain,
    nc.Dunn, indice.Dunn, nc.Hubert, indice.Hubert, nc.sdindex, indice.sdindex, nc.Dindex, indice.Dindex, nc.SDbw, indice.SDbw
    )
    results1<-matrix(c(results),nrow=2,ncol=26)
    resultats <- matrix(c(results),nrow=2,ncol=26,dimnames = list(c("Number_clusters","Value_Index"),
		c("KL","CH","Hartigan","CCC", "Scott", "Marriot", "TrCovW",
		"TraceW", "Friedman", "Rubin", "Cindex", "DB", "Silhouette",
		 "Duda","PseudoT2", "Beale", "Ratkowsky", "Ball", "PtBiserial",
		"Frey", "McClain", "Dunn", 		"Hubert", "SDindex", "Dindex", "SDbw")))

   }


 if (any(indice <= 20)||(indice == 23)||(indice == 24)||(indice == 25))
 {
   resultats <- resultats[,c(indice)]
 }

 if (any(indice == 21)|| (indice == 22))
 {
  indice3 <-indice-1
  resultats <- resultats[,c(indice3)]
 }

 if (any(indice == 26) || (indice == 27) || (indice == 28) || ( indice == 29)|| ( indice == 30))
 {
  indice4 <- indice-4
  resultats <- resultats[,c(indice4)]
 }


  resultats<-round(resultats, digits=4)
  res<-round(res, digits=4)
  resCritical<-round(resCritical, digits=4)

#  if (numberObsAfter != numberObsBefore)
#  {
#     cat(paste(numberObsAfter,"observations were used out of", numberObsBefore ,"possible observations due to missing values."))
#  }

#  if (numberObsAfter == numberObsBefore)
#  {
#     cat(paste("All", numberObsAfter,"observations were used.", "\n", "\n"))
#  }



    ######################## Summary results #####################################


    if(any(indice == 31) || (indice == 32))
    {
      cat("*******************************************************************", "\n")
      cat("* Among all indices:                                               ", "\n")
      BestCluster<-results1[1,]
      c=0
      for(i in min.nc:max.nc)
      {
        vect<-which(BestCluster==i)
        if(length(vect)>0)
        cat("*",length(vect), "proposed", i,"as the best number of clusters", "\n")

        if(c<length(vect))
        {
          j=i
          c<-length(vect)
        }
      }

        cat("\n","                  ***** Conclusion *****                           ", "\n", "\n")
        cat("* According to the majority rule, the best number of clusters is ",j , "\n", "\n", "\n")
        cat("*******************************************************************", "\n")


      ########################## The Best partition    ###################

        if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) ||
          (method == 5) || (method == 6) || (method == 7)||(method == 9))
            partition<- cutree(hc, k=j)

        else
        {
            set.seed(1)
            partition<-kmeans(jeu,j)$cluster
        }

    }


    if (any(indice==1)||(indice==2)||(indice==3)||(indice==4)||(indice==5)||(indice==6)||(indice==7)
        ||(indice==8)||(indice==9)||(indice==10)||(indice==11)||(indice==12)||(indice==13)||(indice==14)
        ||(indice==15)||(indice==16)||(indice==17)||(indice==18)||(indice==19)||(indice==20)
        ||(indice==21)||(indice==22)||(indice==23)||(indice==24)||(indice==25)||(indice==26)
        ||(indice==28)||(indice==30))
    {
      if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) ||
            (method == 5) || (method == 6) || (method == 7) || (method == 9))

        partition<- cutree(hc, k=best.nc)

      else
      {
        set.seed(1)
        partition<-kmeans(jeu,best.nc)$cluster
      }

    }


    #########################  Summary results   ############################



    if ((indice == 14)|| (indice == 15)|| (indice == 16)|| (indice == 20)|| (indice == 31)|| (indice == 32))
    {
      results.final <- list(All.index=res,All.CriticalValues=resCritical,Best.nc=resultats, Best.partition=partition)
    }

    if ((indice == 27)|| (indice == 29))
       results.final <- list(All.index=res)

    if (any(indice==1)||(indice==2)||(indice==3)||(indice==4)||(indice==5)||(indice==6)||(indice==7)
        ||(indice==8)||(indice==9)||(indice==10)||(indice==11)||(indice==12)||(indice==13)
        ||(indice==17)||(indice==18)||(indice==19)||(indice==21)||(indice==22)||(indice==23)||(indice==24)
        ||(indice==25)||(indice==26)||(indice==28)||(indice==30))

      results.final <- list(All.index=res,Best.nc=resultats, Best.partition=partition)



    return(results.final)


}
