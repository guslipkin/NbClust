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

################


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


  ########### Indices.Traces-hartigan - 3e colonne de res ############
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

    if (indice %in% master_index[c(2, 22)]) {
      resultSGAP <- Indice.Gap(
        x = jeu,
        clall = clall,
        reference.distribution = "unif",
        B = 10,
        method = method,
        d = NULL,
        centrotypes = "centroids"
        )
      res[nc - min_nc + 1, "gap"] <- resultSGAP$gap
      resCritical[nc - min_nc + 1, "CritValue_Gap"] <- resultSGAP$diffu
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
        res[nc - min_nc + 1, "sdindex"] <- Index.sdindex(jeu, clmax, cl1)
      }

      if (indice %in% master_index[c(1:2, 31)]) {
        res[nc - min_nc + 1, "dindex"] <- Index.Dindex(cl1, jeu)
      }

      if (indice %in% master_index[c(1:2, 32)]) {
        res[nc - min_nc + 1, "sdbw"] <- Index.SDbw(jeu, cl1)
      }
    }
}

#########################################################################################################
#######################                      Best Number of Clusters                  ###################
#########################################################################################################



   nc.KL<-indice.KL<-0
   if (any(indice == 1) || (indice == 31) || (indice == 32))
	 {
     # KL - The value of u, which maximizes KL(u), is regarded as specifying the number of clusters [ClusterSim package].
     nc.KL <- (min_nc:max_nc)[which.max(res[,1])]
     indice.KL <- max(res[,1],na.rm = TRUE)
     best.nc<-nc.KL
   }

   nc.CH<-indice.CH<-0
   if (any(indice == 2) || (indice == 31) || (indice == 32))
	 {
     # CH - The value of u, which maximizes CH(u), is regarded as specifying the number of clusters [ClusterSim package].
     nc.CH <- (min_nc:max_nc)[which.max(res[,2])]
     indice.CH <- max(res[,2],na.rm = TRUE)
     best.nc<-nc.CH
   }

  nc.CCC<-indice.CCC<-0
  if (any(indice == 4) || (indice == 31) || (indice == 32))
	{
    # CCC - The maximum value accross the hierarchy levels is used to indicate the optimal number of clusters in data [29].
    nc.CCC <- (min_nc:max_nc)[which.max(res[,4])]
    indice.CCC <- max(res[,4],na.rm = TRUE)
    best.nc<-nc.CCC
  }

  nc.DB<-indice.DB<-0
  if (any(indice == 12) || (indice == 31) || (indice == 32))
	{
    # DB - The value of u, which minimizes DB(u), is regarded as specifying the number of clusters [clusterSim package].
    nc.DB <- (min_nc:max_nc)[which.min(res[,12])]
    indice.DB <- min(res[,12],na.rm = TRUE)
    best.nc<-nc.DB
  }

  nc.Silhouette<-indice.Silhouette<-0
  if (any(indice == 13) || (indice == 31) || (indice == 32))
	{
    # SILHOUETTE - The value of u, which maximizes S(u), is regarded as specifying the number of clusters [ClusterSim package].
    nc.Silhouette <- (min_nc:max_nc)[which.max(res[,13])]
    indice.Silhouette <- max(res[,13],na.rm = TRUE)
    best.nc<-nc.Silhouette
  }

  nc.Gap<-indice.Gap<-0
  # GAP - Choose the number of clusters via finding the smallest q such that: Gap(q)=Gap(q+1)-Sq+1 (q=1,\u{85},n-2) [ClusterSim package].
  if (any(indice == 20) || (indice == 32))
  {
    found <- FALSE
    for (ncG in min_nc:max_nc){
      if ((resCritical[ncG-min_nc+1,4] >=0) && (!found)){
          ncGap <- ncG
     	  indiceGap <- res[ncG-min_nc+1,20]
    	  found <- TRUE
    	  }
     }
     if (found){
  	 nc.Gap <- ncGap
  	 indice.Gap <- indiceGap
     best.nc<-nc.Gap
      }else{
  	    nc.Gap <- NA
  	    indice.Gap <- NA
      }

  }

  nc.Duda<-indice.Duda<-0
  # DUDA - Choose the number of clusters via finding the smallest q such that: duda >= critical_value [Duda and Hart (1973)].


   if (any(indice == 14) || (indice == 31) || (indice == 32))
	{

        foundDuda <- FALSE
        for (ncD in min_nc:max_nc)
        {

           if ((res[ncD-min_nc+1,14]>=resCritical[ncD-min_nc+1,1]) && (!foundDuda))
           {
             ncDuda <- ncD
     	       indiceDuda <- res[ncD-min_nc+1,14]
    	       foundDuda <- TRUE
    	     }
        }
        if (foundDuda)
        {
  	      nc.Duda <- ncDuda
  	      indice.Duda <- indiceDuda
          best.nc<-nc.Duda
        }
        else
        {
  	        nc.Duda <- NA
  	        indice.Duda <- NA
        }


   }

  nc.Pseudo<-indice.Pseudo<-0
  # PSEUDOT2 - Chooses the number of clusters via finding the smallest q such that: pseudot2 <= critical_value [SAS User's guide].
	if (any(indice == 15) || (indice == 31) || (indice == 32))
	{

     foundPseudo <- FALSE
     for (ncP in min_nc:max_nc)
       {

      if ((res[ncP-min_nc+1,15]<=resCritical[ncP-min_nc+1,2]) && (!foundPseudo))
        {
          ncPseudo <- ncP
     	  indicePseudo <- res[ncP-min_nc+1,15]
    	  foundPseudo <- TRUE
    	  }
     }
      if (foundPseudo)
        {
  	 nc.Pseudo <- ncPseudo
  	 indice.Pseudo <- indicePseudo
     best.nc<-nc.Pseudo
      }
     else
        {
  	    nc.Pseudo <- NA
  	    indice.Pseudo <- NA
      }
    }


  nc.Beale<-indice.Beale<-0
  	if (any(indice == 16) || (indice == 31) || (indice == 32))
	{
  # BEALE - Chooses the number of clusters via finding the smallest q such that: Fvalue_beale >= 0.1 [Gordon (1999)].
     foundBeale <- FALSE
     for (ncB in min_nc:max_nc)
       {

      if ((resCritical[ncB-min_nc+1,3]>=alphaBeale) && (!foundBeale)){
          ncBeale <- ncB
     	  indiceBeale <- res[ncB-min_nc+1,16]
    	  foundBeale <- TRUE
    	  }
     }
       if (foundBeale){
  	 nc.Beale <- ncBeale
  	 indice.Beale <- indiceBeale
     best.nc<-nc.Beale
      }
     else
       {
  	    nc.Beale <- NA
  	    indice.Beale <- NA
      }
  }


  nc.ptbiserial<-indice.ptbiserial<-0
  if (any(indice == 19) || (indice == 31) || (indice == 32))
	{
    # POINT-BISERIAL - The maximum value was used to suggest the optimal number of clusters in the data [29].
    nc.ptbiserial <- (min_nc:max_nc)[which.max(res[,19])]
    indice.ptbiserial <- max(res[,19],na.rm = TRUE)
    best.nc<-nc.ptbiserial
  }

   foundNC<-foundIndice<-numeric(0)
   nc.Frey<-indice.Frey<-0
   if (any(indice == 21) || (indice == 31) || (indice == 32))
	 {
  # FREY AND VAN GROENEWOUD - The best results occured when clustering was continued until the ratio fell below 1.00 for the last
  #			      series of times. At this point, the cluster level before this series was taken as the optimal partition.
  #			      If the ration never fell below 1.00, a one cluster solution was assumed [29].

     foundFrey <- FALSE
     i<-1
     for (ncF in min_nc:max_nc)
     {

          if (res[ncF-min_nc+1,21] < 1)
          {

              ncFrey <- ncF-1
     	        indiceFrey <- res[ncF-1-min_nc+1,21]
     	        foundFrey <- TRUE
              foundNC[i]<-ncFrey
              foundIndice[i]<-indiceFrey
              i<-i+1

    	     }

     }
     if (foundFrey)
     {
  	   nc.Frey <- foundNC[1]
  	   indice.Frey <- foundIndice[1]
       best.nc<-nc.Frey
     }
      else
      {
  	    nc.Frey <- NA
  	    indice.Frey <- NA
  	    print(paste("Frey index : No clustering structure in this data set"))
      }


   }


   nc.McClain<-indice.McClain<-0
   if (any(indice == 22) || (indice == 31) || (indice == 32))
	{
  # MCCLAIN AND RAO - The minimum value of the index was found to give the best recovery information [29].
  nc.McClain <- (min_nc:max_nc)[which.min(res[,22])]
  indice.McClain <- min(res[,22],na.rm = TRUE)
  best.nc<-nc.McClain

  }

   nc.Gamma<-indice.Gamma<-0
   if (any(indice == 23) || (indice == 31) || (indice == 32))
	{
       # GAMMA - Maximum values were taken to represent the correct hierarchy level [29].
       nc.Gamma <- (min_nc:max_nc)[which.max(res[,23])]
       indice.Gamma <- max(res[,23],na.rm = TRUE)
       best.nc<-nc.Gamma

  }

   nc.Gplus<-indice.Gplus<-0
   if (any(indice == 24) || (indice == 31) || (indice == 32))
	{
  # GPLUS - Minimum values were used to determine the number of clusters in the data [29].
  nc.Gplus  <- (min_nc:max_nc)[which.min(res[,24])]
  indice.Gplus <- min(res[,24],na.rm = TRUE)
  best.nc<-nc.Gplus
  }

   nc.Tau<-indice.Tau<-0
   if (any(indice == 25) || (indice == 31) || (indice == 32))
	{
  # TAU - The maximum value in the hierarchy sequence was taken as indicating the correct number of clusters [29].
  nc.Tau <- (min_nc:max_nc)[which.max(res[,25])]
  indice.Tau <- max(res[,25],na.rm = TRUE)
  best.nc<-nc.Tau
  }


#Some indices need to compute difference between hierarchy levels to identify relevant number of clusters


  if((indice==3)||(indice == 5)||(indice == 6)||(indice == 7)||(indice == 8)||(indice == 9)||(indice == 10)||(indice == 18)||(indice == 27)||(indice == 11)||(indice == 29)||(indice == 31)||(indice == 32))
  {

    DiffLev <- array(0, c(max_nc-min_nc+1,12))
    DiffLev[,1] <- min_nc:max_nc
       for (nc3 in min_nc:max_nc)
      {
        if (nc3==min_nc)
        {
	       DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-NA)   # Hartigan
	       DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-NA)   #Scott
	       DiffLev[nc3-min_nc+1,4] <- abs(res[nc3-min_nc+1,6]-NA)   # Marriot
	       DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-NA)   #Trcovw
	       DiffLev[nc3-min_nc+1,6] <- abs(res[nc3-min_nc+1,8]-NA)   #Tracew
	       DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-NA)   #Friedman
	       DiffLev[nc3-min_nc+1,8] <- abs(res[nc3-min_nc+1,10]-NA)  #Rubin
	       DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-NA)  # Ball
         DiffLev[nc3-min_nc+1,10] <- abs(res[nc3-min_nc+1,27]-NA) # Hubert
         DiffLev[nc3-min_nc+1,12] <- abs(res[nc3-min_nc+1,29]-NA) # D index


	      }
        else
        {
          if(nc3==max_nc)
          {
            DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-res[nc3-min_nc,3])
            DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-res[nc3-min_nc,5])
            DiffLev[nc3-min_nc+1,4] <- abs(res[nc3-min_nc+1,6]-NA) # Marriot
            DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-res[nc3-min_nc,7])  # trcovw
            DiffLev[nc3-min_nc+1,6] <- abs(res[nc3-min_nc+1,8]-NA) #traceW
            DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-res[nc3-min_nc,9])
            DiffLev[nc3-min_nc+1,8] <- abs(res[nc3-min_nc+1,10]-NA) #Rubin
            DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-res[nc3-min_nc,18])
            DiffLev[nc3-min_nc+1,10] <- abs(res[nc3-min_nc+1,27]-NA)
            DiffLev[nc3-min_nc+1,12] <- abs(res[nc3-min_nc+1,29]-NA) # D index


	         }
          else
          {

           DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-res[nc3-min_nc,3]) # Hartigan
	         DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-res[nc3-min_nc,5])
           DiffLev[nc3-min_nc+1,4] <- ((res[nc3-min_nc+2,6]-res[nc3-min_nc+1,6])-(res[nc3-min_nc+1,6]-res[nc3-min_nc,6]))
           DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-res[nc3-min_nc,7])
           DiffLev[nc3-min_nc+1,6] <- ((res[nc3-min_nc+2,8]-res[nc3-min_nc+1,8])-(res[nc3-min_nc+1,8]-res[nc3-min_nc,8]))
           DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-res[nc3-min_nc,9])
           DiffLev[nc3-min_nc+1,8] <- ((res[nc3-min_nc+2,10]-res[nc3-min_nc+1,10])-(res[nc3-min_nc+1,10]-res[nc3-min_nc,10]))
           DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-res[nc3-min_nc,18])
           DiffLev[nc3-min_nc+1,10] <- abs((res[nc3-min_nc+1,27]-res[nc3-min_nc,27]))
           DiffLev[nc3-min_nc+1,12] <-((res[nc3-min_nc+2,29]-res[nc3-min_nc+1,29])-(res[nc3-min_nc+1,29]-res[nc3-min_nc,29])) #Dindex

          }
       }
    }
   }

  nc.Hartigan<-indice.Hartigan<-0
  if (any(indice == 3) || (indice == 31) || (indice == 32))
	{
	# HARTIGAN - The maximum differences between hierarchy levels were taken as indicating the correct number of clusters in the data [29].
	 nc.Hartigan <- DiffLev[,1][which.max(DiffLev[,2])]
	 indice.Hartigan <- max(DiffLev[,2],na.rm = TRUE)
   best.nc<-nc.Hartigan
  }

  nc.Ratkowsky<-indice.Ratkowsky<-0
  if (any(indice == 17) || (indice == 31) || (indice == 32))
	{
  # RATKOWSKY - The optimal number of groups is taken as the level where this criterion exhibits its maximum value [29].
    nc.Ratkowsky <- (min_nc:max_nc)[which.max(res[,17])]
    indice.Ratkowsky <- max(res[,17],na.rm = TRUE)
    best.nc<-nc.Ratkowsky
  }

    nc.cindex<-indice.cindex<-0
    if (any(indice == 11) || (indice == 31) || (indice == 32))
    {
      #CINDEX - The minimum value across the hierarchy levels was used to indicate the optimal number of clusters [29].
      nc.cindex <- (min_nc:max_nc)[which.min(res[,11])]
      indice.cindex <- min(res[,11],na.rm = TRUE)
      best.nc<-nc.cindex
    }

  nc.Scott<-indice.Scott<-0
  if (any(indice == 5) || (indice == 31) || (indice == 32))
	{
 	# SCOTT - The maximum difference between hierarchy levels was used to suggest the correct number of partitions [29].
	 nc.Scott <- DiffLev[,1][which.max(DiffLev[,3])]
	 indice.Scott <- max(DiffLev[,3],na.rm = TRUE)
   best.nc<-nc.Scott
  }

  nc.Marriot<-indice.Marriot<-0
  if (any(indice == 6) || (indice == 31) || (indice == 32))
	{
	# MARRIOT - The maximum difference between successive levels was used to determine the best partition level [29].
	 nc.Marriot <- DiffLev[,1][which.max(DiffLev[,4])]
	 round(nc.Marriot, digits=1)
	 indice.Marriot <- max(DiffLev[,4],na.rm = TRUE)
   best.nc<-nc.Marriot
  }

  nc.TrCovW<-indice.TrCovW<-0
  if (any(indice == 7) || (indice == 31) || (indice == 32))
	{
	nc.TrCovW <- DiffLev[,1][which.max(DiffLev[,5])]
	indice.TrCovW <- max(DiffLev[,5],na.rm = TRUE)
	best.nc<-nc.TrCovW
  }


  nc.TraceW<-indice.TraceW<-0
  if (any(indice == 8) || (indice == 31) || (indice == 32))
	{
  	# TRACE W - To determine the number of clusters in the data, maximum difference scores were used [29].
	  nc.TraceW <- DiffLev[,1][which.max(DiffLev[,6])]
  	indice.TraceW <- max(DiffLev[,6],na.rm = TRUE)
	  best.nc<-nc.TraceW
   }

  nc.Friedman<-indice.Friedman<-0
  if (any(indice == 9) || (indice == 31) || (indice == 32))
	{
  	# FRIEDMAN - The maximum difference in values of trace W-1B criterion was used to indicate the optimal number of clusters [29].
  	nc.Friedman <- DiffLev[,1][which.max(DiffLev[,7])]
  	indice.Friedman <- max(DiffLev[,7],na.rm = TRUE)
  	best.nc<-nc.Friedman
	}

	nc.Rubin<-indice.Rubin<-0
  if (any(indice == 10) || (indice == 31) || (indice == 32))
	{
	  # RUBIN - The difference between levels was used [29].
	  nc.Rubin <- DiffLev[,1][which.min(DiffLev[,8])]
	  indice.Rubin <- min(DiffLev[,8],na.rm = TRUE)
	  best.nc<-nc.Rubin
  }

  nc.Ball<-indice.Ball<-0
  if (any(indice == 18) || (indice == 31) || (indice == 32))
	{
	  # BALL - The largest difference between levels was used to indicate the optimal solution [29].
	  nc.Ball <- DiffLev[,1][which.max(DiffLev[,9])]
	  indice.Ball <- max(DiffLev[,9],na.rm = TRUE)
	  best.nc<-nc.Ball
  }


   nc.Dunn<-indice.Dunn<-0
   if (any(indice == 26) || (indice == 31) || (indice == 32))
	 {
     # Dunn -
     nc.Dunn <- (min_nc:max_nc)[which.max(res[,26])]
     indice.Dunn <- max(res[,26],na.rm = TRUE)
     best.nc<-nc.Dunn
   }


   nc.Hubert<-indice.Hubert<-0
   if (any(indice == 27) || (indice == 31) || (indice == 32))
	 {
	   # Hubert -
     nc.Hubert  <- 0.00
     indice.Hubert  <- 0.00
     #x11()
     par(mfrow = c(1,2))
     plot(x_axis,res[,27], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert Statistic values")))
     plot(DiffLev[,1],DiffLev[,10], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert statistic second differences")))
     cat(paste ("*** : The Hubert index is a graphical method of determining the number of clusters.
                In the plot of Hubert index, we seek a significant knee that corresponds to a
                significant increase of the value of the measure i.e the significant peak in Hubert
                index second differences plot.", "\n", "\n"))
     }

   nc.sdindex<-indice.sdindex<-0
   if (any(indice == 28) || (indice == 31) || (indice == 32))
	 {
     # SD -
     nc.sdindex <- (min_nc:max_nc)[which.min(res[,28])]
     indice.sdindex<- min(res[,28],na.rm = TRUE)
     best.nc<-nc.sdindex
   }


    nc.Dindex<-indice.Dindex<-0
    if (any(indice == 29) || (indice == 31) || (indice == 32))
	  {

     nc.Dindex <- 0.00
     indice.Dindex<- 0.00
     #x11()
     par(mfrow = c(1,2))
     plot(x_axis,res[,29], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Dindex Values")))
     plot(DiffLev[,1],DiffLev[,12], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Second differences Dindex Values")))
     cat(paste ("*** : The D index is a graphical method of determining the number of clusters.
                In the plot of D index, we seek a significant knee (the significant peak in Dindex
                second differences plot) that corresponds to a significant increase of the value of
                the measure.", "\n", "\n"))
    }

    nc.SDbw<-indice.SDbw<-0
    if (any(indice == 30) || (indice == 31) || (indice == 32))
	   {
       # SDbw -
       nc.SDbw <- (min_nc:max_nc)[which.min(res[,30])]
       indice.SDbw<- min(res[,30],na.rm = TRUE)
       best.nc<-nc.SDbw
     }



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
