setwd("D:/Rt/Gab/")
require(wnl)
dPK35 = read.csv("PK35.csv", skip=1)
colnames(dPK35) = c("TIME", "DV", "Category") ; dPK35

TypVal = c(1.8, 500)
OMEGA = diag(c(1.8^2, 500^2))
SIGMA = 1
invOM = solve(OMEGA)

dTime = seq(0, 456, by=24)
nDose = length(dTime) ; nDose

TIME = dPK35[dPK35$Category==1, "TIME"]
nTIME = length(TIME)

AMT = 200
ETAi = c(0, 0)

fPK35 = function(THETA)
{
  CL = THETA[1]
  V  = THETA[2]
  ETAi = c(CL, V) - TypVal
  
  tRes = matrix(rep(0, nTIME*nDose), nrow=nTIME, ncol=nDose)
  for (i in 1:nDose) {
    cTime = dTime[i]
    ty = AMT/V*exp(-CL/V*(TIME - cTime))
    ty[ty < 0] = 0 # no case
    tRes[,i] = ty 
  }
  Cp = rowSums(tRes)
  Ri = dPK35[dPK35$Category==1, "DV"] - Cp
  return(nTIME*log(SIGMA) + sum(Ri*Ri)/SIGMA + t(ETAi) %*% invOM %*% ETAi)
}

optim(c(2, 350), fPK35, method="L-BFGS-B")

##
fPK35a = function(THETA)
{
  CL = THETA[1]
  V  = THETA[2]
  
  tRes = matrix(rep(0, nTIME*nDose), nrow=nTIME, ncol=nDose)
  for (i in 1:nDose) {
    cTime = dTime[i]
    ty = AMT/V*exp(-CL/V*(TIME - cTime))
    ty[ty < 0] = 0 # no case
    tRes[,i] = ty 
  }
  Cp = rowSums(tRes)
  Ri = dPK35[dPK35$Category==1, "DV"] - Cp
  return(sum(Ri*Ri/Cp) + (CL - 1.8)^2/1.8^2 + (V - 500)^2/500^2)
}

optim(c(2, 350), fPK35a, method="L-BFGS-B")

## Poisson error

fPK35b = function(THETA)
{
  CL = THETA[1]
  V  = THETA[2]
  ETAi = c(CL, V) - TypVal
  
  tRes = matrix(rep(0, nTIME*nDose), nrow=nTIME, ncol=nDose)
  for (i in 1:nDose) {
    cTime = dTime[i]
    ty = AMT/V*exp(-CL/V*(TIME - cTime))
    ty[ty < 0] = 0 # no case
    tRes[,i] = ty 
  }
  Cp = rowSums(tRes)
  Ri = dPK35[dPK35$Category==1, "DV"] - Cp
  return(sum(log(Cp) + Ri*Ri/Cp) + t(ETAi) %*% invOM %*% ETAi)
}
fPK35b(c(1.8, 500))
fPK35b(c(2, 350))
fPK35b(c(5.7, 119.6))
fPK35b(c(3.267, 646.39))
optim(c(2, 350), fPK35b, method="L-BFGS-B")

## compare with parameters from OLS
fPK35c = function(THETA)
{
  CL = THETA[1]
  V  = THETA[2]
  
  tRes = matrix(rep(0, nTIME*nDose), nrow=nTIME, ncol=nDose)
  for (i in 1:nDose) {
    cTime = dTime[i]
    ty = AMT/V*exp(-CL/V*(TIME - cTime))
    ty[ty < 0] = 0 # no case
    tRes[,i] = ty 
  }
  Cp = rowSums(tRes)
  Ri = dPK35[dPK35$Category==1, "DV"] - Cp
  return(sum(Ri*Ri))
}

optim(c(2, 350), fPK35c, method="L-BFGS-B")


## compare with parameters from OLS
fPK35d = function(THETA)
{
  CL = THETA[1]
  V  = THETA[2]
  
  tRes = matrix(rep(0, nTIME*nDose), nrow=nTIME, ncol=nDose)
  for (i in 1:nDose) {
    cTime = dTime[i]
    ty = AMT/V*exp(-CL/V*(TIME - cTime))
    ty[ty < 0] = 0 # no case
    tRes[,i] = ty 
  }
  Cp = rowSums(tRes)
  Ri = dPK35[dPK35$Category==1, "DV"] - Cp
  return(sum(Ri*Ri/Cp))
}

optim(c(2, 350), fPK35d, method="L-BFGS-B")
