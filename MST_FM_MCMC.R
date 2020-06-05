MST_FM_MCMC <- function(FPCAout, # output from function MST_FM_decomposition, including the follows:
                                # mu: estimated mean function (vector of length T)
                                # psi1: estimated first-level eigenfunctions (matrix of dimension L*T)
                                # psi2: estimated second-level eigenfunctions (matrix of dimension M*T)
                                # sigma2: estimated measurement error variance, used as an initial value in the MCMC step (scalar)
                                # lambda: estimated lambda (stablizing component of the second-level eigenvalues) (vector of length M)
                                # tau2: estimated tau^2 (region-specific component of the second-level eigenvalues) (vector of length n)
                       
                       data, # data.frame in long format with eight labeled columns (described below)
                       # and row length equal to the length of the vectorized observations across all 
                       # regions and facilities
                       # DATA.FRAME COLUMNS: 
                       # rid: region IDs (vector of length T*sum(Ni) )
                       # fid: facility IDs (vector of length T*sum(Ni))
                       # y: hospitalization rate data (vector of length T*sum(Ni))
                       # t: follow-up time (vector of length T*sum(Ni)) 
                       
                       AdjMat  # Adjacency matrix from the map (0-1 matrix of dimension n*n)
){
  
  #############################################################################
  ## Description: Function for MCMC estimation (estimation step 6 in Table 1) of MST-FM model described in "Multilevel Modeling of Spatially Nested
  ##              Functional Data: Spatial Variation in Hospitalization Rates among U.S. Dialysis Patients", including estimation 
  ##              of spatial variance parameters, measurement error variance and region- and facility-specific PC scores. 
  ## Definition:  n: number of regions, Ni: number of facilities in region i, T: number of time points, 
  ##              L: number of first-level eigencomponents, M: number of second-level eigencomponents
  ## Args:        see above
  ## Returns:     list()
  ##              alpha: posterior samples of spatial variance parameter (matrix of dimension L*2000)
  ##              nv: posterior samples of spatial correlation parameter (vector of length 2000)
  ##              sigma2: posterior samples of measurement error variance (vector of length 2000)
  ##              xi: posterior samples of region-specific PC scores (matrix of dimension n*2000)
  ##              zeta: posterior samples of facility-specific PC scores (matrix of dimension sum(Ni)*2000)
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
  
  # Define functions
  logit <- function(x){log(x/(1-x))}
  
  unlogit <- function(x){exp(x)/(1+exp(x))}
  
  repeat.row <- function(m, times){
    return(m[rep(1:nrow(m), times = times),])
  }
  
  # Format data
  df <- data
  # Number of regions n
  nregion <- length(unique(df$rid))
  
  # Number of time points T
  ngrid <- length(unique(df$t))
  gridPoints <- unique(df$t)
  
  # Number of facilities per region Ni
  nFac <- aggregate(df$fid, by = list(df$rid), FUN = length)[,2] / ngrid
  numFac <- sum(nFac) # Total number of facilities
  
  # Extract estimated parameters from function MST_FM_decomposition
  muEst <- FPCAout$mu
  nei1 <- dim(FPCAout$psi1)[2]
  nei2 <- dim(FPCAout$psi2)[2]
  psi1.t <- FPCAout$psi1
  psi2.t <- FPCAout$psi2
  spsi1 <- c(sum(psi1.t[,1]^2), sum(psi1.t[,2]^2))
  spsi2 <- c(sum(psi2.t[,1]^2), sum(psi2.t[,2]^2))
  df$dt <- round(df$t * (ngrid-1)) + 1
  X.c <- df$y - FPCAout$mu[df$dt]
  lambda <- FPCAout$lambda
  tau2 <- FPCAout$tau2
  sigmaEst <- FPCAout$sigma2
  
  # Initial values for MCMC iteration
  xi=rep(0, nei1*nregion)
  zeta = rep(0, nei2*numFac)
  nv=.9
  alpha = c(8,4) 
  n.sample=1 
  v.accept=.4
  qv=.01
  
  # list for saving MCMC results
  mod=list(xi=xi,zeta=zeta,nv=nv,alpha=alpha, sigma = sigmaEst, seed=.Random.seed,
           v.accept=v.accept,qv=qv,n.sample=n.sample,total=n.sample)
  
  # prior for nv
  av <- 9
  bv <- 1
  # prior for alpha
  aa <- c(2,2)
  ba <- c(8,4)
  # prior for sigma2
  as <- 2
  bs <- sigmaEst
  
  # Number of MCMC samples
  S_inc=2500
  n.sample=mod$n.sample+S_inc
  oldT=mod$n.sample #this should always be 1
  
  # Create arrays/vectors for our parameters, plug in current values into the first slot
  xi <- array(dim=c(nei1*nregion, n.sample))
  xi[,1:mod$n.sample] <- mod$xi
  zeta <- array(dim=c(nei2*numFac, n.sample))
  zeta[,1:mod$n.sample] <- mod$zeta
  nv <- array(dim=n.sample)
  nv[1:mod$n.sample] <- mod$nv
  alpha <- array(dim = c(nei1, n.sample))
  alpha[,1:mod$n.sample] <- mod$alpha
  sigma.MC <- NULL
  sigma.MC[1:mod$n.sample] <- mod$sigma
  qv <- mod$qv
  total <- mod$total
  
  # Make the CAR precision matrix and get its determinant
  W <- AdjMat
  D <- rowSums(W)
  DvW <- diag(D) - nv[oldT] * W
  detDvW <- det(DvW) 
  
  ## Fix nv 
  # nv <- rep(0.9, n.sample)
  # skip.v <- TRUE
  
  ## Estimate nv
  skip.v <- FALSE
  
  for(it in (oldT+1):n.sample){
    
    ############# Update region-specific PC score xi's ############# 
    xiMat <- matrix(xi[, it-1], nrow = nregion)
    zetaMat <- matrix(zeta[, it-1], nrow = numFac)
    faceff <- rowSums(repeat.row(psi2.t, numFac) * repeat.row(zetaMat, rep(ngrid, numFac)))
    XX <- X.c - faceff
    
    for(i in 1:nei1){
      if(nei1==2){
        speff <- rep(psi1.t[,-i], numFac) * rep(xiMat[,-i], nFac * ngrid)
      } else{
        speff <- rowSums(repeat.row(psi1.t[,-i], numFac) * repeat.row(xiMat[,-i], nFac * ngrid))
      }
      Y <- (XX - speff) * psi1.t[df$dt,i]
      mu1 <- aggregate(Y, by = list(df$rid), FUN = sum)[,2]
      mu1 <- mu1/sigma.MC[it-1] + nv[it-1] * W %*% xiMat[,i] / alpha[i, it-1]
      v.xi <- alpha[i,it-1] * nFac * spsi1[i] + D * sigma.MC[it-1]
      v.xi <- alpha[i,it-1] * sigma.MC[it-1] / v.xi
      mu.xi <- v.xi * mu1
      xiMat[,i] <- rnorm(nregion, mu.xi, sqrt(v.xi))
      xi[(nregion*(i-1)+1):(nregion*i), it] <- xiMat[,i]
    }
    
    ############# Update facility-specific PC score zeta's #############
    xiMat <- matrix(xi[, it], nrow = nregion)
    zetaMat <- matrix(zeta[, it-1], nrow = numFac)
    speff <- rowSums(repeat.row(psi1.t, numFac) * repeat.row(xiMat, nFac * ngrid))
    XX <- X.c - speff
    
    for(i in 1:nei2){
      if(nei2==2){
        faceff <- rep(psi2.t[,-i], numFac) * rep(zetaMat[,-i], rep(ngrid, numFac))
      } else{
        faceff <- rowSums(repeat.row(psi2.t[,-i], numFac) * repeat.row(zetaMat[,-i], rep(ngrid, numFac)))
      }
      Y <- (XX - faceff) * psi2.t[df$dt,i]
      mu1 <- aggregate(Y, by = list(df$fid), FUN = sum)[,2]
      v.zeta <- lambda[i] * tau2 * spsi2[i] + sigma.MC[it-1]
      v.zeta <- lambda[i] * sigma.MC[it-1] * tau2 / v.zeta
      v.zeta <- rep(v.zeta, nFac)
      mu.zeta <- v.zeta * mu1 / sigma.MC[it-1]
      zeta[(numFac*(i-1)+1):(numFac*i), it] <- rnorm(numFac, mu.zeta, sqrt(v.zeta))
    }
    zetaMat <- matrix(zeta[, it], nrow = numFac)
    
    ############## Update spatial variance parameter alpha's ###########
    for(i in 1:nei1){
      bap <- ba[i] + 1/2 * t(xiMat[,i]) %*% DvW %*% xiMat[,i]
      alpha[i,it] <- 1/rgamma(1,aa[i] + nregion/2, bap)
    }
    #sig2[,it]=sig2[,it-1]
    
    ############## Update measurement error variance sigma^2 ###########
    y <- as.vector(t(repeat.row(xiMat %*% t(psi1.t), nFac))) + as.vector(t(zetaMat %*% t(psi2.t))) + muEst[df$dt]
    b.post <- bs + sum((df$y-y)^2)/2
    sigma.MC[it] <- 1/rgamma(1, as + numFac*ngrid/2, b.post)
    
    
    
    ############## Update spatial correlation parameter nv #############
    if(skip.v==FALSE){ 
      lam <- logit(nv[it-1])
      lams <- rnorm(1,lam,qv) 
      nvs <- unlogit(lams)
      DvWs <- diag(D) - nvs*W
      detDvWs <- det(DvWs)
      r1 <- nei1/2 * log(detDvWs/detDvW) 
      
      r2 <- 0
      for(i in 1:nei1){
        r2 <- r2 + t(xiMat[,i]) %*% W %*% xiMat[,i] / alpha[i,it]
      }
      
      r2 <- (nvs - nv[it-1])/2 * r2
      
      r3 <- exp(lams-lam) * ( (1+exp(lam))/(1+exp(lams)) )^2 *  # this is the jacobian of the logit transformation
        (nvs/nv[it-1])^(av-1) * ((1-nvs)/(1-nv[it-1]) )^(bv-1)
      
      r <- exp(r1+r2)*r3
      
      accept <- ifelse(r>runif(1),1,0)
      if(accept==1){
        nv[it]=nvs
        DvW=DvWs
        detDvW=detDvWs
      }else{
        nv[it] <- nv[it-1]
      }
    }
    
  }
  # drop the first 500 samples as burn-in and save the last 2000 samples for further analysis
  mod=list(xi=xi[,-(1:501)],zeta=zeta[,-(1:501)],nv=nv[-(1:501)],alpha=alpha[,-(1:501)],
           sigma2=sigma.MC[-(1:501)],seed=.Random.seed)
  return(mod)
}


