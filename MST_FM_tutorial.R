## MST_FM_tutorial.R
#############################################################################
## Description: A step-by-step implementation of MST-FM and the associated  
## procedures described in "Multilevel Modeling of Spatially Nested Functional Data:
## Spatiotemporal Patterns of Hospitalization Rates in the U.S. Dialysis Population".  
#############################################################################
## Functions implemented: 
## MST_FM_simulation.R, MST_FM_decomposition.R, MST_FM_MCMC.R, MST_FM_inference.R, MST_FM_correctedInference.R
#############################################################################
## Tutorial Outline:
## 1. Simulate hospitalization rate outcome data (MST_FM_simulation.R)
## 2. Perform MST-FM estimation (MST_FM_decomposition.R, MST_FM_MCMC.R)
## 3. Prediction and inference on multilevel hospitalization trajectories (MST_FM_inference.R)
## 4. Corrected inference on region-specific hospitalization trajectories via bootstrap (MST_FM_correctedInference.R)
## 5. Visualization of MST-FM results
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


#############################################################################
# 1. Simulate hospitalization rate outcome data
#############################################################################

# NOTE: Generating one dataset with 49 regions and 4-20 facilities per region will take approximately three minutes.

# Simulate one dataset from the simulation design described in Web Appendix E with 49 regions and 4-20 facilities per region
data.G <- MST_FM_simulation(numRegion = 2, numFacility = 1, sigma2 = 0.1)  # MST_FM_simulation.R

# Data frame used for estimation and inference
data <- data.G[[1]]

# Adjacency matrix from the map
Adj.Mat <- data.G[[2]]

# Save the underlying true values of the generated data
data.T <- data.G[[3]]

#############################################################################
# 2. Perform MST-FM estimation
#############################################################################

# NOTE: Performing MST-FM estimation steps 1-5 (FPCA decomposition) with 49 regions and 4-20 facilities per region will take approximately two minutes.

FPCAout <- MST_FM_decomposition(data = data)  # MST_FM_decompostion.R

# NOTE: Performing MST-FM estimation step 6 (MCMC estimation) with 49 regions and 4-20 facilities per region will take approximately 15 minutes.

MCMCout <- MST_FM_MCMC(FPCAout = FPCAout, data = data, AdjMat = Adj.Mat) # MST_FM_MCMC.R

#############################################################################
# 3. Prediction and inference on multilevel hospitalization trajectories
#############################################################################

PREDout <- MST_FM_inference(FPCAout = FPCAout, data = data, MCMCout = MCMCout)

#############################################################################
# 4. Corrected inference on region-specific hospitalization trajectories via bootstrap
#############################################################################

# NOTE: The corrected inference procedure via bootstrap may take several hours
#       depending on the size of the data set and processor speed. 

Traj.corrected <- MST_FM_correctedInference(FPCAout, data, MCMCout, nboot = 100, AdjMat = Adj.Mat)

#############################################################################
# 5. Visualization of MST-FM results
#############################################################################  

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

## Align eigenfunctions
# Note: the estimated eigenfunctions may flip signs. Here we match the sign of the
# estimated eigenfunctions to the true eigenfunctions for plotting  
MSDE.eifun <- function(x,y){
  err1 <- trapz(gridPoints, (x-y)^2)
  err2 <- trapz(gridPoints, (x+y)^2)
  return(c(min(err1, err2),which.min(c(err1, err2))))
}

mu.t <- mu(gridPoints)
psi1.1.t <- psi1.1(gridPoints)
psi1.2.t <- psi1.2(gridPoints)
psi2.1.t <- psi2.1(gridPoints)
psi2.2.t <- psi2.2(gridPoints)

mu.Est <- FPCAout$mu
psi1.1.Est <- FPCAout$psi1[,1]
psi1.2.Est <- FPCAout$psi1[,2]
psi2.1.Est <- FPCAout$psi2[,1]
psi2.2.Est <- FPCAout$psi2[,2]


MSDEeifun1.2 <- MSDE.eifun(psi1.1.Est, psi1.1.t)
MSDEeifun2.2 <- MSDE.eifun(psi1.2.Est, psi1.2.t)

MSDEeifun1.fac <- MSDE.eifun(psi2.1.Est, psi2.1.t)
MSDEeifun2.fac <- MSDE.eifun(psi2.2.Est, psi2.2.t)

fun1.sign.2 <- MSDEeifun1.2[2] * (-2) + 3
psi1.1.Est <- psi1.1.Est * fun1.sign.2
fun2.sign.2 <- MSDEeifun2.2[2] * (-2) + 3
psi1.2.Est <- psi1.2.Est * fun2.sign.2

fun1.sign.fac <- MSDEeifun1.fac[2] * (-2) + 3
psi2.1.Est <- psi2.1.Est * fun1.sign.fac
fun2.sign.fac <- MSDEeifun2.fac[2] * (-2) + 3
psi2.2.Est <- psi2.2.Est * fun2.sign.fac

# Plot estimates of mean function and eigenfunctions
par(mfrow = c(3,2))
# Mean function (add labels for true and estimated functions)
plot(gridPoints, mu.t,"l", ylim = c(1,2.8), xaxs = "i", main = "(a)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t : months",ylab=expression(widehat(mu)(t)), line=2, cex.lab=1.6)
lines(gridPoints, mu.Est, col = "grey", lwd = 2, lty = 1)

# First-level eigenfunctions
plot(gridPoints, psi1.1.t,"l", ylim = c(-1.4,1.4), xaxs = "i", main = "(b)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t : months",ylab=expression(paste(widehat(psi)[1]^(1),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi1.1.Est, col = "grey", lwd = 2, lty = 1)

plot(gridPoints, psi1.2.t,"l", ylim = c(-1.4,1.5), xaxs = "i", main = "(c)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t : months",ylab=expression(paste(widehat(psi)[2]^(1),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi1.2.Est, col = "grey", lwd = 2, lty = 1)

# Second-level eigenfunctions
plot(gridPoints, psi2.1.t,"l", ylim = c(-1.9,1.9), xaxs = "i", main = "(d)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t : months",ylab=expression(paste(widehat(psi)[1]^(2),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi2.1.Est, col = "grey", lwd = 2, lty = 1)

plot(gridPoints, psi2.2.t,"l", ylim = c(-1.1,2.3), xaxs = "i", main = "(e)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t : months",ylab=expression(paste(widehat(psi)[2]^(2),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi2.2.Est, col = "grey", lwd = 2, lty = 1)


# Region-specific predicted trajectory for the first region
R.ID <- 1
# True hospitalization rate trajectory for the first region
R.traj.T <- mu.t + data.T$r.eff1[data.T$rid==R.ID][1:ngrid] + data.T$r.eff2[data.T$rid==R.ID][1:ngrid]
# Predicted hospitalization rate trajectory for the first region
R.traj.Est <- PREDout$R.traj.Est[R.ID,]
# Pointwise 95% confidence intervals
R.traj.SD <- PREDout$R.traj.SD[R.ID,]
R.traj.upper <- R.traj.Est + 1.96 * R.traj.SD
R.traj.lower <- R.traj.Est - 1.96 * R.traj.SD
# Visualization of the region-specific trajectory prediction for the first region
ylim1 <- min(R.traj.lower)
ylim2 <- max(R.traj.upper)
plot(gridPoints, R.traj.T,"l", xaxs = "i", ylim = c(ylim1, ylim2), main = "Region-specific trajectory prediction",cex.main=2,xlab = "", ylab = "")
title(xlab = "t : months",ylab="Hospitalization rate", line=2, cex.lab=1.6)
lines(gridPoints, R.traj.Est, col = "grey", lwd = 2, lty = 1)
lines(gridPoints, R.traj.upper, col = "grey", lwd = 2, lty = 2)
lines(gridPoints, R.traj.lower, col = "grey", lwd = 2, lty = 2)

# The corrected region-specific predicted trajectory for the first region
C.traj.Est <- Traj.corrected[[1]][R.ID,]
C.traj.upper <- Traj.corrected[[3]][R.ID,]
C.traj.lower <- Traj.corrected[[4]][R.ID,]
ylim1 <- min(C.traj.lower)
ylim2 <- max(C.traj.upper)
# Visualization of the region-specific trajectory prediction for the first region
plot(gridPoints, R.traj.T,"l", ylim = c(ylim1,ylim2), xaxs = "i", main = "Region-specific trajectory prediction (corrected)",cex.main=2,xlab = "", ylab = "")
title(xlab = "t : months",ylab="Hospitalization rate", line=2, cex.lab=1.6)
lines(gridPoints, C.traj.Est, col = "grey", lwd = 2, lty = 1)
lines(gridPoints, C.traj.upper, col = "grey", lwd = 2, lty = 2)
lines(gridPoints, C.traj.lower, col = "grey", lwd = 2, lty = 2)

# Facility-specific predicted trajectory for the first facility from the first region
F.ID <- 1
# True hospitalization rate trajectory for the first facility from the first region
F.traj.T <- mu.t + data.T$r.eff1[data.T$fid==F.ID][1:ngrid] + data.T$r.eff2[data.T$fid==F.ID][1:ngrid] +
  data.T$f.eff1[data.T$fid==F.ID] + data.T$f.eff2[data.T$fid==F.ID]
# Predicted hospitalization rate trajectory for the first facility from the first region
F.traj.Est <- PREDout$F.traj.Est[F.ID,]
# Pointwise 95% confidence intervals 
F.traj.SD <- PREDout$F.traj.SD[F.ID,]
F.traj.upper <- F.traj.Est + 1.96 * F.traj.SD
F.traj.lower <- F.traj.Est - 1.96 * F.traj.SD
# Visualization of the facility-specific trajectory prediction for the first facility from the first region
ylim1 <- min(F.traj.lower)
ylim2 <- max(F.traj.upper)
plot(gridPoints, F.traj.T,"l", ylim = c(ylim1,ylim2), xaxs = "i", main = "Facility-specific trajectory",cex.main=2,xlab = "", ylab = "")
title(xlab = "t : months",ylab="Hospitalization rate", line=2, cex.lab=1.6)
lines(gridPoints, F.traj.Est, col = "grey", lwd = 2, lty = 1)
lines(gridPoints, F.traj.upper, col = "grey", lwd = 2, lty = 2)
lines(gridPoints, F.traj.lower, col = "grey", lwd = 2, lty = 2)
