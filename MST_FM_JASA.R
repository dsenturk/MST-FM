# Code for generating Table 2 and Table S4 in paper and supplement
source("MST_FM_simulation.R")
source("MST_FM_decomposition.R")
source("MST_FM_MCMC.R")
source("MST_FM_inference.R")
source("MST_FM_correctedinference.R")

numRun <- 50 # Number of simulation runs (200 runs are used for generating results in Table 2 and Table S4)
FacilityPerRegion <- 1 # Number of facilities per region (scalar, input 1 if you want 4-20 facilities per region, input 2 if you want 10-30 facilities per region)

# Functions for calculating MSDE for mean function, predicted trajectories and eigenfunctions
MSDE.fun <- function(x,y){
  err.MSDE <- trapz(gridPoints, (x-y)^2) / trapz(gridPoints, y^2)
  return(err.MSDE)
}

MSDE.eifun <- function(x,y){
  err1 <- trapz(gridPoints, (x-y)^2)
  err2 <- trapz(gridPoints, (x+y)^2)
  return(c(min(err1, err2),which.min(c(err1, err2))))
}

# True values of model parameters
# Define the grid points used for the mean function and eigenfunctions
ngrid <- 24 # 2 year follow up
gridPoints <- seq(0,1,length.out = ngrid)

# True mean function and eigenfunctions
# mean function
mu <- function(x){
  return(.6*(x-1)^2+1.6)
}

# first-level eigenfunctions
psi1.1 <- function(x){
  return(sqrt(2)*sin(2*pi*x))
}
psi1.2 <- function(x){
  return(sqrt(2)*cos(2*pi*x))
}

# second-level eigenfunctions
psi2.1 <- function(x){
  return(sqrt(3)*(2*x-1))
}
psi2.2 <- function(x){
  return(sqrt(5)*(6*x^2-6*x+1))
}

mu.t <- mu(gridPoints)
psi1.1.t <- psi1.1(gridPoints)
psi1.2.t <- psi1.2(gridPoints)
psi2.1.t <- psi2.1(gridPoints)
psi2.2.t <- psi2.2(gridPoints)

alpha.T <- c(1, 0.25)
nv.T <- 0.9
sigma.T <- 0.1

# Storing results from multiple runs
MSDE.mu <- c()
MSDE.psi1 <- c()
MSDE.psi2 <- c()
MSE.alpha <- c()
MSE.nv <- c()
MSE.sigma2 <- c()
MSE.xi1 <- c()
MSE.xi2 <- c()
MSE.zeta1 <- c()
MSE.zeta2 <- c()
MSE.lambda1 <- c()
MSE.lambda2 <- c()
MSE.lambda1.S <- c()
MSE.lambda2.S <- c()
MSE.lambda1.M <- c()
MSE.lambda2.M <- c()
MSE.lambda1.L <- c()
MSE.lambda2.L <- c()
MSDE.fac <- c()
CP.fac <- c()
Len.fac <- c()
MSDE.reg <- c()
MSDE.regS <- c()
MSDE.regM <- c()
MSDE.regL <- c()
CP.reg <- c()
CP.regS <- c()
CP.regM <- c()
CP.regL <- c()
Len.reg <- c()
Len.regS <- c()
Len.regM <- c()
Len.regL <- c()
MSDE.reg.C <- c()
MSDE.regS.C <- c()
MSDE.regM.C <- c()
MSDE.regL.C <- c()
CP.reg.C <- c()
CP.regS.C <- c()
CP.regM.C <- c()
CP.regL.C <- c()
Len.reg.C <- c()
Len.regS.C <- c()
Len.regM.C <- c()
Len.regL.C <- c()
for(r in 1:numRun){
  #############################################################################
  # 1. Simulate hospitalization rate outcome data
  #############################################################################
  
  # Simulate one dataset from the simulation design described in Web Appendix E with 49 regions and 4-20 facilities per region
  data.G <- MST_FM_simulation(numRegion = 2, numFacility = FacilityPerRegion, sigma2 = .1)  # MST_FM_simulation.R
  
  # Data frame used for estimation and inference
  data <- data.G[[1]]
  
  # Adjacency matrix from the map
  Adj.Mat <- data.G[[2]]
  
  # Save the underlying true values of the generated data
  data.T <- data.G[[3]]
  
  #############################################################################
  # 2. Perform MST-FM estimation
  #############################################################################

  FPCAout <- MST_FM_decomposition(data = data)  # MST_FM_decompostion.R
  MCMCout <- MST_FM_MCMC(FPCAout = FPCAout, data = data, AdjMat = Adj.Mat) # MST_FM_MCMC.R
  
  #############################################################################
  # 3. Prediction and inference on multilevel hospitalization trajectories
  #############################################################################
  
  PREDout <- MST_FM_inference(FPCAout = FPCAout, data = data, MCMCout = MCMCout)
  
  #############################################################################
  # 4. Corrected inference on region-specific hospitalization trajectories via bootstrap
  #############################################################################
  
  # NOTE: the corrected inference procedure are only used for cases with 49 regions
  Traj.corrected <- MST_FM_correctedInference(FPCAout, data, MCMCout, nboot = 100, AdjMat = Adj.Mat)
  
  #############################################################################
  # 5. Collect results from the current run and combine results from runs
  #############################################################################
  # MSDE
  # Extract estimates
  nregion <- length(unique(data$rid))
  numFac <- length(unique(data$fid))
  mu.Est <- FPCAout$mu
  psi1.1.Est <- FPCAout$psi1[,1]
  psi1.2.Est <- FPCAout$psi1[,2]
  psi2.1.Est <- FPCAout$psi2[,1]
  psi2.2.Est <- FPCAout$psi2[,2]
  alphaEst <- rowMeans(MCMCout$alpha)
  nvEst <- mean(MCMCout$nv)
  sigmaEst <- mean(MCMCout$sigma2)
  
  MSDEeifun1.2 <- MSDE.eifun(psi1.1.Est, psi1.1.t)
  MSDEeifun2.2 <- MSDE.eifun(psi1.2.Est, psi1.2.t)
  MSDEeifun1.fac <- MSDE.eifun(psi2.1.Est, psi2.1.t)
  MSDEeifun2.fac <- MSDE.eifun(psi2.2.Est, psi2.2.t)
  
  xi1Est <- rowMeans(MCMCout$xi[1:nregion,]) * (MSDEeifun1.2[2] * (-2) + 3)
  xi2Est <- rowMeans(MCMCout$xi[(nregion + 1): (2*nregion),]) * (MSDEeifun2.2[2] * (-2) + 3)
  zeta1Est <- rowMeans(MCMCout$zeta[1:numFac,]) * (MSDEeifun1.fac[2] * (-2) + 3)
  zeta2Est <- rowMeans(MCMCout$zeta[(numFac + 1): (2*numFac),]) * (MSDEeifun2.fac[2] * (-2) + 3)
  
  lambda1Est <- FPCAout$tau2 * FPCAout$lambda[1]
  lambda2Est <- FPCAout$tau2 * FPCAout$lambda[2]
  
  # Mean function 
  MSDE.mu <- c(MSDE.mu, MSDE.fun(mu.Est, mu.t))
  # Eigenfunctions
  MSDE.psi1 <- cbind(MSDE.psi1, c(MSDEeifun1.2[1], MSDEeifun2.2[1]))
  MSDE.psi2 <- cbind(MSDE.psi2, c(MSDEeifun1.fac[1], MSDEeifun2.fac[1]))
  # MSE for time-invariant parameters
  MSE.alpha <- cbind(MSE.alpha, (alphaEst - alpha.T)^2)
  MSE.nv <- c(MSE.nv, (nvEst - nv.T)^2)
  MSE.sigma2 <- c(MSE.sigma2, (sigmaEst - sigma.T)^2)
  
  # True values of PC scores and second-level eigenvalues
  region.index <- (c(1, diff(data$rid))) == 1
  facility.index <- (c(1, diff(data$fid))) == 1
  xi1.T <- data.T$xi1[region.index]
  xi2.T <- data.T$xi2[region.index]
  zeta1.T <- data.T$zeta1[facility.index]
  zeta2.T <- data.T$zeta2[facility.index]
  lambda1.T <- data.T$f.sig1[region.index]
  lambda2.T <- data.T$f.sig2[region.index]
  MSE.xi1 <- c(MSE.xi1, mean((xi1Est - xi1.T)^2))
  MSE.xi2 <- c(MSE.xi2, mean((xi2Est - xi2.T)^2))
  MSE.zeta1 <- c(MSE.zeta1, mean((zeta1Est - zeta1.T)^2))
  MSE.zeta2 <- c(MSE.zeta2, mean((zeta2Est - zeta2.T)^2))
  # Region size
  Size.region <- data.T$r.size[region.index]

  MSE.lambda1 <- c(MSE.lambda1, mean((lambda1Est - lambda1.T)^2))
  MSE.lambda2 <- c(MSE.lambda2, mean((lambda2Est - lambda2.T)^2))
  # Lambda_{im} for small region
  MSE.lambda1.S <- c(MSE.lambda1.S, mean((lambda1Est[Size.region==1] - lambda1.T[Size.region==1])^2))
  MSE.lambda2.S <- c(MSE.lambda2.S, mean((lambda2Est[Size.region==1] - lambda2.T[Size.region==1])^2))
  # Lambda_{im} for medium region
  MSE.lambda1.M <- c(MSE.lambda1.M, mean((lambda1Est[Size.region==2] - lambda1.T[Size.region==2])^2))
  MSE.lambda2.M <- c(MSE.lambda2.M, mean((lambda2Est[Size.region==2] - lambda2.T[Size.region==2])^2))
  # Lambda_{im} for large region
  MSE.lambda1.L <- c(MSE.lambda1.L, mean((lambda1Est[Size.region==3] - lambda1.T[Size.region==3])^2))
  MSE.lambda2.L <- c(MSE.lambda2.L, mean((lambda2Est[Size.region==3] - lambda2.T[Size.region==3])^2))
  
  # MSDE, coverage probability and length of 95% confidence intervals for multilevel predicted trajectories
  # Region-specific predicted trajectories 
  R.traj.cov <- matrix(0,nrow = nregion, ncol = ngrid)
  R.traj.len <- matrix(0,nrow = nregion, ncol = ngrid)
  R.traj.MSDE <- NULL
  for(R.ID in 1:nregion){
    # Underlying true region-specific trajectory
    R.traj.T <- mu.t + data.T$r.eff1[data.T$rid==R.ID][1:ngrid] + data.T$r.eff2[data.T$rid==R.ID][1:ngrid]
    # Predicted region-specific trajectory
    R.traj.Est <- PREDout$R.traj.Est[R.ID,]
    # MSDE
    R.traj.MSDE[R.ID] <- MSDE.fun(R.traj.Est, R.traj.T)
    # Pointwise 95% confidence intervals for region-specific trajectory 
    R.traj.SD <- PREDout$R.traj.SD[R.ID,]
    R.traj.upper <- R.traj.Est + 1.96 * R.traj.SD
    R.traj.lower <- R.traj.Est - 1.96 * R.traj.SD
    # Compute coverage and length of the pointwise 95% confidence intervals
    R.traj.cov[R.ID,] <- (R.traj.T > R.traj.lower) & (R.traj.T < R.traj.upper)
    R.traj.len[R.ID,] <- 2 * 1.96 * R.traj.SD
  }
  
  # MSDE
  MSDE.reg <- c(MSDE.reg, mean(R.traj.MSDE))
  MSDE.regS <- c(MSDE.regS, mean(R.traj.MSDE[Size.region==1]))
  MSDE.regM <- c(MSDE.regM, mean(R.traj.MSDE[Size.region==2]))
  MSDE.regL <- c(MSDE.regL, mean(R.traj.MSDE[Size.region==3]))
  
  # Coverage probability
  # All regions
  CP.reg <- c(CP.reg, mean(R.traj.cov))
  # Small regions
  CP.regS <- c(CP.regS, mean(R.traj.cov[Size.region==1,]))
  # Medium regions
  CP.regM <- c(CP.regM, mean(R.traj.cov[Size.region==2,]))
  # Large regions
  CP.regL <- c(CP.regL, mean(R.traj.cov[Size.region==3,]))
  # Length 
  # All regions
  Len.reg <- c(Len.reg, mean(R.traj.len))
  # Small regions
  Len.regS <- c(Len.regS, mean(R.traj.len[Size.region==1,]))
  # Medium regions
  Len.regM <- c(Len.regM, mean(R.traj.len[Size.region==2,]))
  # Large regions
  Len.regL <- c(Len.regL, mean(R.traj.len[Size.region==3,]))
  
  # MSDE, coverage and length of the corrected region-specific predicted trajectories
  R.traj.cov.C <- matrix(0,nrow = nregion, ncol = ngrid)
  R.traj.len.C <- matrix(0,nrow = nregion, ncol = ngrid)
  R.traj.MSDE.C <- NULL
  for(R.ID in 1:nregion){
    # Underlying true region-specific trajectory
    R.traj.T <- mu.t + data.T$r.eff1[data.T$rid==R.ID][1:ngrid] + data.T$r.eff2[data.T$rid==R.ID][1:ngrid]
    # Predicted region-specific trajectory
    R.traj.Est <- Traj.corrected[[1]][R.ID,]
    # MSDE
    R.traj.MSDE.C[R.ID] <- MSDE.fun(R.traj.Est, R.traj.T)
    # Pointwise 95% confidence intervals for region-specific trajectory
    R.traj.SD <- Traj.corrected[[2]][R.ID,]
    R.traj.upper <- Traj.corrected[[3]][R.ID,]
    R.traj.lower <- Traj.corrected[[4]][R.ID,]
    # Compute coverage and length of the pointwise 95% confidence intervals
    R.traj.cov.C[R.ID,] <- (R.traj.T > R.traj.lower) & (R.traj.T < R.traj.upper)
    R.traj.len.C[R.ID,] <- 2 * 1.96 * R.traj.SD
  }

  # MSDE
  MSDE.reg.C <- c(MSDE.reg.C, mean(R.traj.MSDE.C))
  MSDE.regS.C <- c(MSDE.regS.C, mean(R.traj.MSDE.C[Size.region==1]))
  MSDE.regM.C <- c(MSDE.regM.C, mean(R.traj.MSDE.C[Size.region==2]))
  MSDE.regL.C <- c(MSDE.regL.C, mean(R.traj.MSDE.C[Size.region==3]))

  # Coverage probability
  # All regions
  CP.reg.C <- c(CP.reg.C, mean(R.traj.cov.C))
  # Small regions
  CP.regS.C <- c(CP.regS.C, mean(R.traj.cov.C[Size.region==1,]))
  # Medium regions
  CP.regM.C <- c(CP.regM.C, mean(R.traj.cov.C[Size.region==2,]))
  # Large regions
  CP.regL.C <- c(CP.regL.C, mean(R.traj.cov.C[Size.region==3,]))
  # Length
  # All regions
  Len.reg.C <- c(Len.reg.C, mean(R.traj.len.C))
  # Small regions
  Len.regS.C <- c(Len.regS.C, mean(R.traj.len.C[Size.region==1,]))
  # Medium regions
  Len.regM.C <- c(Len.regM.C, mean(R.traj.len.C[Size.region==2,]))
  # Large regions
  Len.regL.C <- c(Len.regL.C, mean(R.traj.len.C[Size.region==3,]))
  # 
  # Facility-specific predicted trajectories 
  F.traj.cov <- matrix(0,nrow = numFac, ncol = ngrid)
  F.traj.len <- matrix(0,nrow = numFac, ncol = ngrid)
  F.traj.MSDE <- NULL
  for(F.ID in 1:numFac){
    F.traj.T <- mu.t + data.T$r.eff1[data.T$fid==F.ID][1:ngrid] + data.T$r.eff2[data.T$fid==F.ID][1:ngrid] +
      data.T$f.eff1[data.T$fid==F.ID] + data.T$f.eff2[data.T$fid==F.ID]
    # Predicted region-specific trajectory
    F.traj.Est <- PREDout$F.traj.Est[F.ID,]
    # MSDE
    F.traj.MSDE[F.ID] <- MSDE.fun(F.traj.Est, F.traj.T)
    # Pointwise 95% confidence intervals for facility-specific trajectory 
    F.traj.SD <- PREDout$F.traj.SD[F.ID,]
    F.traj.upper <- F.traj.Est + 1.96 * F.traj.SD
    F.traj.lower <- F.traj.Est - 1.96 * F.traj.SD
    # Compute coverage and length of the pointwise 95% confidence intervals
    F.traj.cov[F.ID,] <- (F.traj.T > F.traj.lower) & (F.traj.T < F.traj.upper)
    F.traj.len[F.ID,] <- 1.96 * 2 * F.traj.SD
  }
  
  MSDE.fac <- c(MSDE.fac, mean(F.traj.MSDE))
  CP.fac <- c(CP.fac, mean(F.traj.cov))
  Len.fac <- c(Len.fac, mean(F.traj.len))
}
# Table 2
Tab2 <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(Tab2) <- c("MSDE.mu","MSDE.psi1.1","MSDE.psi1.2","MSDE.psi2.1","MSDE.psi2.2","MSE.alpha1","MSE.alpha2","MSE.nv","MSE.sigma2",
                  "MSE.xi1","MSE.xi2","MSE.zeta1","MSE.zeta2","MSE.lambda1","MSE.lambda1S","MSE.lambda1M","MSE.lambda1L","MSE.lambda2","MSE.lambda2S","MSE.lambda2M","MSE.lambda2L")
Tab2$MSDE.mu <- mean(MSDE.mu)
Tab2$MSDE.psi1.1 <- mean(MSDE.psi1[1,])
Tab2$MSDE.psi1.2 <- mean(MSDE.psi1[2,])
Tab2$MSDE.psi2.1 <- mean(MSDE.psi2[1,])
Tab2$MSDE.psi2.2 <- mean(MSDE.psi2[2,])
Tab2$MSE.alpha1 <- mean(MSE.alpha[1,])
Tab2$MSE.alpha2 <- mean(MSE.alpha[2,])
Tab2$MSE.nv <- mean(MSE.nv)
Tab2$MSE.sigma2 <- mean(MSE.sigma2)
Tab2$MSE.xi1 <- mean(MSE.xi1)
Tab2$MSE.xi2 <- mean(MSE.xi2)
Tab2$MSE.zeta1 <- mean(MSE.zeta1)
Tab2$MSE.zeta2 <- mean(MSE.zeta2)
Tab2$MSE.lambda1 <- mean(MSE.lambda1)
Tab2$MSE.lambda1S <- mean(MSE.lambda1.S)
Tab2$MSE.lambda1M <- mean(MSE.lambda1.M)
Tab2$MSE.lambda1L <- mean(MSE.lambda1.L)
Tab2$MSE.lambda2 <- mean(MSE.lambda2)
Tab2$MSE.lambda2S <- mean(MSE.lambda2.S)
Tab2$MSE.lambda2M <- mean(MSE.lambda2.M)
Tab2$MSE.lambda2L <- mean(MSE.lambda2.L)

# Output table
write.csv(Tab2, file = paste("Table2 NumRegion", nregion, "NumFacility", FacilityPerRegion, "NoiseLevel", sigma.T, ".csv"))

# Table S4 from supplement
TabS4 <- data.frame(matrix(ncol = 15, nrow = 1))
colnames(TabS4) <- c("MSDE.fac","MSDE.reg","MSDE.regS","MSDE.regM","MSDE.regL","CP.fac","CP.reg","CP.regS","CP.regM","CP.regL",
                    "Len.fac","Len.reg","Len.regS","Len.regM","Len.regL")

TabS4$MSDE.fac <- mean(MSDE.fac)
TabS4$MSDE.reg <- mean(MSDE.reg)
TabS4$MSDE.regS <- mean(MSDE.regS)
TabS4$MSDE.regM <- mean(MSDE.regM)
TabS4$MSDE.regL <- mean(MSDE.regL)
TabS4$CP.fac <- mean(CP.fac)
TabS4$CP.reg <- mean(CP.reg)
TabS4$CP.regS <- mean(CP.regS)
TabS4$CP.regM <- mean(CP.regM)
TabS4$CP.regL <- mean(CP.regL)
TabS4$Len.fac <- mean(Len.fac)
TabS4$Len.reg <- mean(Len.reg)
TabS4$Len.regS <- mean(Len.regS)
TabS4$Len.regM <- mean(Len.regM)
TabS4$Len.regL <- mean(Len.regL)

# Output table
write.csv(TabS4, file = paste("TableS4 NumRegion", nregion, "NumFacility", FacilityPerRegion, "NoiseLevel", sigma.T, ".csv"))

# MSDE, coverage probability and length of 95% confidence intervals for corrected region-specific predicted trajectories 
# NOTE: the corrected inference procedure are only used for cases with 49 regions

TabS4.C <- data.frame(matrix(ncol = 12, nrow = 1))
colnames(TabS4.C) <- c("MSDE.reg.C","MSDE.regS.C","MSDE.regM.C","MSDE.regL.C","CP.reg.C","CP.regS.C","CP.regM.C","CP.regL.C",
                     "Len.reg.C","Len.regS.C","Len.regM.C","Len.regL.C")

TabS4.C$MSDE.reg.C <- mean(MSDE.reg.C)
TabS4.C$MSDE.regS.C <- mean(MSDE.regS.C)
TabS4.C$MSDE.regM.C <- mean(MSDE.regM.C)
TabS4.C$MSDE.regL.C <- mean(MSDE.regL.C)
TabS4.C$CP.reg.C <- mean(CP.reg.C)
TabS4.C$CP.regS.C <- mean(CP.regS.C)
TabS4.C$CP.regM.C <- mean(CP.regM.C)
TabS4.C$CP.regL.C <- mean(CP.regL.C)
TabS4.C$Len.reg.C <- mean(Len.reg.C)
TabS4.C$Len.regS.C <- mean(Len.regS.C)
TabS4.C$Len.regM.C <- mean(Len.regM.C)
TabS4.C$Len.regL.C <- mean(Len.regL.C)

# Output table
write.csv(TabS4.C, file = paste("TableS4Corrected NumRegion", nregion, "NumFacility", FacilityPerRegion, "NoiseLevel", sigma.T, ".csv"))
