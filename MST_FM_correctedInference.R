MST_FM_correctedInference <- function(FPCAout, # output from function MST_FM_decomposition, including the follows:
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
                              
                              MCMCout, # output from function MST_FM_MCMC, including the follows:
                              # alpha: posterior samples of spatial variance parameter (matrix of dimension L*2000)
                              # nv: posterior samples of spatial correlation parameter (vector of length 2000)
                              # sigma2: posterior samples of measurement error variance (vector of length 2000)
                              # xi: posterior samples of region-specific PC scores (matrix of dimension n*2000)
                              # zeta: posterior samples of facility-specific PC scores (matrix of dimension sum(Ni)*2000)
                              
                              nboot = 100, # Number of bootstrap samples used in the correction procedure
                              
                              AdjMat  # Adjacency matrix from the map (0-1 matrix of dimension n*n)
                              
){
  
  #############################################################################
  ## Description: Function for obtaining corrected prediction and inference for region-specific trajectories for applications with 
  ##              small number of regions described in Section 2.3 and Web Appendix B Table S1
  ## Definition:  n: number of regions, Ni: number of facilities in region i, T: number of time points, 
  ##              L: number of first-level eigencomponents, M: number of second-level eigencomponents
  ## Args: see above
  ## Returns:     list()
  ##              C.Ests.IV: Corrected region-specific trajectory predictions (matrix of dimension n*T)
  ##              C.sd.IV:  Pointwise standard errors of the Corrected region-specific trajectory predictions (matrix of dimension n*T)
  ##              C.CI.upper: Corrected upper bound of the 95% CI for region-specific trajectory predictions (matrix of dimension n*T)
  ##              C.CI.lower: Corrected lower bound of the 95% CI for region-specific trajectory predictions (matrix of dimension n*T)
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
  numFac <- sum(nFac) # Total number of facilities
  
  # Adjacency matrix from the map 
  W <- AdjMat
  D <- rowSums(W)
  
  # Extract estimated parameters from functions MST_FM_decomposition and MST_FM_MCMC
  muEst <- FPCAout$mu
  psi1.t <- FPCAout$psi1
  psi2.t <- FPCAout$psi2
  lambdaEst <- FPCAout$lambda
  tau2 <- FPCAout$tau2
  nvEst <- mean(MCMCout$nv)
  sigmaEst.MC <- mean(MCMCout$sigma2)
  alpha1Est.MC <- mean(MCMCout$alpha[1,])
  alpha2Est.MC <- mean(MCMCout$alpha[2,])
  
  # List for saving bootstrap results 
  C.Est.BS <- list()
  C.Var.BS <- list()
  for(boot in 1:nboot){
    ###########################################################################
    # Step 0: Construct parametric bootstrap samples 
    ###########################################################################
    # Create data.frame to store dataset
    df.bs <- data.frame(matrix(ncol = 8, nrow = numFac * ngrid))
    colnames(df.bs) <- c("rid","fid","y","t","c.eff1","c.eff2","f.eff1","f.eff2")
    df.bs$rid <- rep(1:nregion, nFac*ngrid)
    df.bs$fid <- rep(1:numFac, each = ngrid)
    df.bs$t <- rep(gridPoints,numFac)
    df.bs$dt <- round(df.bs$t * (ngrid-1)) + 1
    covmat <- solve(diag(D) - nvEst * W)
    
    # First-level region-specific deviation 
    # Generate region-specific PC scores from multivariate normal distribution
    county.eff1 <- mvrnorm(1,rep(0,nregion), covmat * alpha1Est.MC)
    county.eff2 <- mvrnorm(1,rep(0,nregion), covmat * alpha2Est.MC)
    df.bs$c.eff1 <- county.eff1[df.bs$rid] * psi1.t[df.bs$dt,1]
    df.bs$c.eff2 <- county.eff2[df.bs$rid] * psi1.t[df.bs$dt,2]
    
    # Second-level facility-specific deviation
    # Generate facility-specific PC scores from normal distribution
    fac.eff1 <- rnorm(numFac, rep(0, numFac), rep(sqrt(tau2 * lambdaEst[1]), nFac))
    fac.eff2 <- rnorm(numFac, rep(0, numFac), rep(sqrt(tau2 * lambdaEst[2]), nFac))
    df.bs$f.eff1 <- fac.eff1[df.bs$fid] * psi2.t[df.bs$dt,1]
    df.bs$f.eff2 <- fac.eff2[df.bs$fid] * psi2.t[df.bs$dt,2]
    
    # Generate outcome
    measure.err <- rnorm(dim(df.bs)[1], 0, sqrt(sigmaEst.MC))
    df.bs$y <- muEst[df.bs$dt] + df.bs$c.eff1 + df.bs$c.eff2 + df.bs$f.eff1 + df.bs$f.eff2 + measure.err
    
    df.bs <- df.bs[,1:4]
    ###########################################################################
    # Step 1: Estimate model parameters based on the bootstrap sample and obtain region-specific trajectory predictions and variances
    ###########################################################################
    
    # Estimate model parameters based on the bootstrap sample
    FPCAout.bs <- MST_FM_decomposition(data = df.bs)  # MST_FM_decomposition.R
    MCMCout.bs <- MST_FM_MCMC(FPCAout = FPCAout.bs, data = df, AdjMat = W) # MST_FM_MCMC.R
    
    # Obtain region-specific trajectory predictions and variances
    C.traj.Est.bs <- matrix(0,nrow = nregion, ncol = ngrid)
    C.traj.Var.bs <- matrix(0,nrow = nregion, ncol = ngrid)
    for(i in 1:nregion){
      # BLUP
      C.traj.Est.bs[i,] <- FPCAout.bs$mu + FPCAout.bs$psi1[,1] * mean(MCMCout.bs$xi[i,]) + FPCAout.bs$psi1[,2] * mean(MCMCout.bs$xi[i+nregion,])
      Cov.Est <- cov(t(MCMCout.bs$xi[c(i,i+nregion),]))
      C.traj.Var.bs[i,] <- diag(FPCAout.bs$psi1 %*% Cov.Est %*% t(FPCAout.bs$psi1))
    }
    C.Est.BS[[boot]] <- C.traj.Est.bs 
    C.Var.BS[[boot]] <- C.traj.Var.bs
  }
  

  
 ###########################################################################
 # Step 2: Obtain corrected region-specific trajectory predictions and their variance using iterative expections and variances 
 ###########################################################################
  
  #Iterative expectation and variance
  C.Ests.IV <- matrix(0,nrow = nregion, ncol = ngrid)
  C.sd.IV <- matrix(0,nrow = nregion, ncol = ngrid)
  C.CI.upper <- matrix(0,nrow = nregion, ncol = ngrid)
  C.CI.lower <- matrix(0,nrow = nregion, ncol = ngrid)
  for(i in 1:nregion){
    C.Est.b <- c()
    C.Var.b <- c()
    for(b in 1:nboot){
      C.Est.b <- cbind(C.Est.b, C.Est.BS[[b]][i,])
      C.Var.b <- cbind(C.Var.b, C.Var.BS[[b]][i,])
    }
    # Iterative expections
    C.Est.IV <- rowMeans(C.Est.b)
    # Iterative variances
    var1 <- rowMeans(C.Var.b)
    var2 <- apply(C.Est.b,1,var)
    # Corrected variances
    sd.IV <- sqrt(var1 + var2)
    C.Ests.IV[i,] <- C.Est.IV
    C.sd.IV[i,] <- sd.IV
    ###########################################################################
    # Step 3: Obtain corrected pointsiwe 95% CIs for region-specific trajectories
    ###########################################################################
    C.CI.upper[i,] <- C.Est.IV + 1.96 * sd.IV
    C.CI.lower[i,] <- C.Est.IV - 1.96 * sd.IV
  }
  
  # Construct output
  out <- list(C.Ests.IV = C.Ests.IV, C.sd.IV = C.sd.IV, C.CI.upper = C.CI.upper, C.CI.lower = C.CI.lower)
  return(out)
}