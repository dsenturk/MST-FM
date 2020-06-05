MST_FM_inference <- function(FPCAout, # output from function MST_FM_decomposition, including the follows:
                              # Definition: n: number of regions, Ni: number of facilities in region i, T: number of time points
                              # L: number of first-level eigencomponents, M: number of second-level eigencomponents
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
                              
                              MCMCout # output from function MST_FM_MCMC, including the follows:
                              # alpha: posterior samples of spatial variance parameter (matrix of dimension L*2000)
                              # nv: posterior samples of spatial correlation parameter (vector of length 2000)
                              # sigma2: posterior samples of measurement error variance (vector of length 2000)
                              # xi: posterior samples of region-specific PC scores (matrix of dimension n*2000)
                              # zeta: posterior samples of facility-specific PC scores (matrix of dimension sum(Ni)*2000)
                              
){
  
  #############################################################################
  ## Description: Function for obtaining prediction and inference for multilevel trajectories described in Section 2.2
  ## Definition:  n: number of regions, Ni: number of facilities in region i, T: number of time points, 
  ##              L: number of first-level eigencomponents, M: number of second-level eigencomponents
  ## Args: see above
  ## Returns:     list()
  ##              R.traj.Est: Region-specific trajectory predictions (matrix of dimension n*T)
  ##              R.traj.SD:  Pointwise standard errors of the region-specific trajectory predictions (matrix of dimension n*T)
  ##              F.traj.Est: Facility-specific trajectory predictions (matrix of dimension sum(Ni)*T)
  ##              F.traj.SD:  Pointwise standard errors of the facility-specific trajectory predictions (matrix of dimension sum(Ni)*T)
  ############################################################################# 
  
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
  
  # Extract estimated parameters from functions MST_FM_decomposition and MST_FM_MCMC
  muEst <- FPCAout$mu
  nei1 <- dim(FPCAout$psi1)[2]
  nei2 <- dim(FPCAout$psi2)[2]
  psi1.t <- FPCAout$psi1
  psi2.t <- FPCAout$psi2
  xi <- MCMCout$xi
  zeta <- MCMCout$zeta
  
  
  R.traj.Est <- matrix(0,nrow = nregion, ncol = ngrid)
  R.traj.SD <- matrix(0,nrow = nregion, ncol = ngrid)
  for(i in 1:nregion){
    # Predicted region-specific trajectory
    R.traj.Est[i,] <- muEst + psi1.t[,1] * mean(xi[i,]) + psi1.t[,2] * mean(xi[i+nregion,])
    # Pointwise 95% confidence intervals for region-specific trajectory 
    Cov.Est <- cov(t(xi[c(i,i+nregion),]))
    R.traj.SD[i,] <- sqrt(diag(psi1.t %*% Cov.Est %*% t(psi1.t)))
  }
  
  # Facility-specific predicted trajectories 
  F.traj.Est <- matrix(0,nrow = numFac, ncol = ngrid)
  F.traj.SD <- matrix(0,nrow = numFac, ncol = ngrid)
  for(i in 1:numFac){
    # Find region id for the facility
    F.c.id <- (df$rid[df$fid==i])[1]
    # Predicted region-specific trajectory
    F.traj.Est[i,] <- muEst + psi1.t[,1] * mean(xi[F.c.id,]) + psi1.t[,2] * mean(xi[F.c.id+nregion,]) + 
      psi2.t[,1] * mean(zeta[i,]) + psi2.t[,2] * mean(zeta[i+numFac,])
    # Pointwise 95% confidence intervals for facility-specific trajectory 
    Cov.Est <- cov(t(rbind(xi[c(F.c.id,F.c.id+nregion),], zeta[c(i,i+numFac),])))
    F.traj.SD[i,] <- sqrt(diag(cbind(psi1.t,psi2.t) %*% Cov.Est %*% t(cbind(psi1.t,psi2.t))))
  }
  
  out <- list(R.traj.Est = R.traj.Est, R.traj.SD = R.traj.SD, F.traj.Est = F.traj.Est, F.traj.SD = F.traj.SD)
  return(out)
}
  
