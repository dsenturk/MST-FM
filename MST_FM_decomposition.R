MST_FM_decomposition <- function(data      # data.frame in long format with eight labeled columns (described below)
                              # and row length equal to the length of the vectorized observations across all 
                              # regions, facilities and time points (NOTE: all facilities must have the same set of time points)
                              # DATA.FRAME COLUMNS: 
                              # rid: region IDs (vector of length T*sum(Ni)) 
                              # fid: facility IDs (vector of length T*sum(Ni))
                              # y: hospitalization rate data (vector of length T*sum(Ni))
                              # t: follow-up time (vector of length T*sum(Ni)) 
){
  
  #############################################################################
  ## Description: Function for FPCA decomposition (estimation steps 1-5 in Table 1) of MST-FM model described in "Multilevel Modeling of Spatially Nested
  ##              Functional Data: Spatial Variation in Hospitalization Rates among U.S. Dialysis Patients", including estimation 
  ##              of mean function, multilevel eigenfunctions and eigenvalues. 
  ## Definition:  n: number of regions, Ni: number of facilities in region i, T: number of time points, 
  ##              L: number of first-level eigencomponents, M: number of second-level eigencomponents
  ## Args:        see above
  ## Returns:     list()
  ##              mu: estimated mean function (vector of length T)
  ##              psi1: estimated first-level eigenfunctions (matrix of dimension L*T) 
  ##              psi2: estimated second-level eigenfunctions (matrix of dimension M*T) 
  ##              sigma2: estimated measurement error variance, used as an initial value in the MCMC step (scalar)
  ##              lambda: estimated lambda (stablizing component of the second-level eigenvalues) (vector of length M)
  ##              tau2: estimated tau^2 (region-specific component of the second-level eigenvalues) (vector of length n) 
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("MASS", "caTools", "locpol", "KernSmooth", "fANCOVA", "mgcv", "mvtnorm", "spdep")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  
  # Load packages  
  library(MASS)
  library(caTools)
  library(locpol)
  library(KernSmooth)
  library(fANCOVA)
  library(mgcv)
  library(mvtnorm)
  library(spdep)
  
  # Format data
  df <- data
  # Number of regions n
  nregion <- length(unique(df$rid))
  
  # Number of time points T
  ngrid <- length(unique(df$t))
  gridPoints <- unique(df$t)
  
  # Number of facilities per region Ni
  nFac <- aggregate(df$fid, by = list(df$rid), FUN = length)[,2] / ngrid
  
  
  ###########################################################################
  # Implement the estimation algorithm steps 1-5 as described in Section 2.2 Table 1
  ###########################################################################
  
  # Step 1: estimation of the mean function
  rr <- smooth.spline(df$t, df$y)
  muEst <- rr$y
  
  # Step 2: estimation of between and within facility raw covariance function
  # center data
  df$dt <- round(df$t * (ngrid-1)) + 1
  X.c <- df$y - muEst[df$dt]
  
  # Calculate raw between facility covariance 
  Xmat <- matrix(X.c, nrow = ngrid)
  covmat.c <- matrix(0,ngrid,ngrid)
  Fac.Cov <- list()
  for(regionid in 1:nregion){
    X.c1 <- X.c[df$rid==regionid]
    covmat.c1 <- matrix(0,ngrid,ngrid)
    Xmat1 <- matrix(X.c1, nrow = ngrid)
    numFac.region <- nFac[regionid]
    for(i in 1:(numFac.region-2)){
      Xi <- Xmat1[,i]
      XXi <- rowSums(Xmat1[,(i+1):numFac.region])
      covmat.c1 <- covmat.c1 + tcrossprod(Xi,XXi) + tcrossprod(XXi,Xi)
    }
    covmat.c1 <- covmat.c1 + tcrossprod(Xmat1[,numFac.region-1], Xmat1[,numFac.region]) + tcrossprod(Xmat1[,numFac.region], Xmat1[,numFac.region-1])
    Fac.Cov[[regionid]] <- covmat.c1 / numFac.region / (numFac.region-1) # Between facility covariance for regions
    covmat.c <- covmat.c +   Fac.Cov[[regionid]] # Pool data from all regions
  }
  
  # Calculate raw within facility covariance
  G.WithinCov <- list()
  covmat.fac <- matrix(0,ngrid,ngrid)
  for(regionid in 1:nregion){
    X.c1 <- X.c[df$rid==regionid]
    covmat.c1 <- matrix(0,ngrid,ngrid)
    Xmat1 <- matrix(X.c1, nrow = ngrid)
    G.withinFac <- tcrossprod(Xmat1, Xmat1)
    numFac.region <- nFac[regionid]
    G.withinFac <- G.withinFac / numFac.region
    # subtract the between facility covariance from total covariance
    G.Fac <-  G.withinFac - Fac.Cov[[regionid]]
    G.WithinCov[[regionid]] <- G.Fac
    covmat.fac <- covmat.fac + G.Fac # Pool data from all regions
  }
  
  
  # Steps 3 and 4: Obtain estimators of between and within facility covariance functions and employ FPCA to estimate eigenfunctions
  x0 <- rep(gridPoints,each = ngrid)
  x1 <- rep(gridPoints, ngrid)
  
  # Bivariate penalized spline smoothing of between facility covariance function
  cov.mat.c.s <- gam(as.vector(covmat.c / nregion) ~ te(x0, x1, k=10, bs = "ps"))
  cov.c.s <- matrix(cov.mat.c.s$fitted.values, nrow = ngrid)
  
  # FPCA on between facility covariance function
  eigen_temp <- eigen(cov.c.s, symmetric = TRUE)
  eigen_temp$values <- eigen_temp$values[which(eigen_temp$values > 0)]  # Obtain positive eigenvalues
  eigen_temp$vectors <- eigen_temp$vectors[, 1:length(eigen_temp$values)]  # Obtain eigenvectors associated with positive eigenvalues
  # eigenfunctions
  for(e in 1:length(eigen_temp$values)){  # Normalize the eigenvalues over the domain
    normal.factor <- trapz(gridPoints, eigen_temp$vectors[, e]^2)
    eigen_temp$vectors[, e] <- eigen_temp$vectors[, e] / sqrt(normal.factor)
    eigen_temp$values[e] <- eigen_temp$values[e] * normal.factor
  }
  L <- length(which(cumsum(eigen_temp$values) / sum(eigen_temp$values) < .90)) + 1 # Number of first-level eigen components

  # Estimated first-level eigen functions
  eifun1s.2 <- eigen_temp$vectors[,1]
  eifun2s.2 <- eigen_temp$vectors[,2] 
  
  # Bivariate penalized spline smoothing of within facility covariance function
  # Remove diagonal entries
  diag(covmat.fac) <- NA
  cov.mat.fac.s <- gam(as.vector(covmat.fac) ~ te(x0, x1, k=10, bs = "ps"))
  grids2d <- data.frame(x0 = x0, x1 = x1)
  cov.fac.s <- matrix(predict(cov.mat.fac.s, newdata = grids2d), nrow = ngrid)
  cov.fac.s <- (cov.fac.s + t(cov.fac.s)) / 2 #  Symmetrize covariance function
  
  # FPCA on within facility covariance function
  eigen_temp <- eigen(cov.fac.s, symmetric = TRUE)
  eigen_temp$values <- eigen_temp$values[which(eigen_temp$values > 0)]  # Obtain positive eigenvalues
  eigen_temp$vectors <- eigen_temp$vectors[, 1:length(eigen_temp$values)]  # Obtain eigenvectors associated with positive eigenvalues
  # eigenfunctions
  for(e in 1:length(eigen_temp$values)){  # Normalize the eigenvalues over the domain
    normal.factor <- trapz(gridPoints, eigen_temp$vectors[, e]^2)
    eigen_temp$vectors[, e] <- eigen_temp$vectors[, e] / sqrt(normal.factor)
    eigen_temp$values[e] <- eigen_temp$values[e] * normal.factor
  }
  M <- length(which(cumsum(eigen_temp$values) / sum(eigen_temp$values) < .90)) + 1 # Number of second-level eigen components
  
  # Estimated second-level eigen functions
  eifun1s.fac <- eigen_temp$vectors[,1]
  eifun2s.fac <- eigen_temp$vectors[,2]
  
  # Estimated lambda (stablizing component of the second-level eigenvalues)
  lambdaEst <- eigen_temp$values[1:M]
  
  # Step 5: estimation of region-specific eigenvalues
  sigmas <- c()  # Estimate measurement error variance (initial values for MCMC)
  eivalue <- matrix(0, nrow = nregion, M)
  for(regionid in 1:nregion){
    G.withinCov.region <- G.WithinCov[[regionid]]
    diag.within <- diag(G.withinCov.region) # Extract diagonals of within facility covariances
    loess_diag <- suppressWarnings(loess.as(gridPoints, diag.within, degree = 1,
                                            criterion = "gcv", user.span = NULL, plot = F))  # Smooth diagonal entries
    # Remove diagonal entries
    diag(G.withinCov.region) <- NA
    cov.mat.fac.s <- gam(as.vector(G.withinCov.region) ~ te(x0, x1, k=10, bs = "ps"))
    cov.fac.s <- matrix(predict(cov.mat.fac.s, newdata = grids2d), nrow = ngrid)
    cov.fac.s <- (cov.fac.s + t(cov.fac.s)) / 2 #  Symmetrize covariance function
    # Calculate difference between the diagonals of the raw within facility covariances and the smoothed within facility covariances
    sigmas <- c(sigmas, mean(loess_diag$fitted - diag(cov.fac.s)))
    for(e in 1:M){
      eivec1 <- eigen_temp$vectors[,e]
      surf1 <- tcrossprod(eivec1, eivec1) * cov.fac.s
      # Double integral
      l <- apply(surf1, MARGIN = 1, FUN = trapz, x = gridPoints)
      eivalue[regionid,e] <-   trapz(gridPoints,l)
    }
  }
  sigmaEst <- mean(sigmas)  #Estimated measurement error variance (initial values for MCMC)
  
  # Estimate tau_i^2
  lambdaMat <- matrix(rep(lambdaEst, nregion), nrow = nregion, byrow = TRUE)
  taui <- rowMeans(eivalue / lambdaMat) 
  
  # Modify negative tau_i^2
  taui.1 <- eivalue[,1] / lambdaEst[1]
  taui[taui<0] <- taui.1[taui<0]
  eivalue <- matrix(rep(taui,M),nrow = nregion) * lambdaMat
  if(min(taui) < 0){
    taui[taui < 0] <- 1e-6 # If tau_i^2 is estimated to be negative, replace it with 1e-6
    warning("Negative variance")
  }
  
  # Construct output
  psi1.t <- cbind(eifun1s.2, eifun2s.2)
  psi2.t <- cbind(eifun1s.fac, eifun2s.fac)
  out <- list(psi1 = psi1.t, psi2 = psi2.t, mu = muEst, sigma2 = sigmaEst, lambda = lambdaEst, tau2 = taui)
  return(out)
}



