setwd("D:/Rt/Gab")
PK02 = read.csv("PK02.csv", skip=1)
colnames(PK02) = c("TIME","DV")
DV2 = c(4.7, 4.48, 4.22, 3.87, 3.57, 2.97, 2.25, 1.74, 1.02, 0.77, 0.61, 0.36, 0.2)
PK02 = cbind(PK02, DV2)
PK02

## Plot
plot(DV2 ~ TIME, data=PK02, type="o")
lines(DV ~ TIME, data=PK02, type="o", col="red")
dev.new()
plot(log(DV2) ~ TIME, data=PK02, type="o")
lines(log(DV) ~ TIME, data=PK02, type="o", col="red")

## NCA
library(NonCompart)
R1 = sNCA(PK02$TIME, PK02$DV2, dose=100, adm="Bolus", doseUnit="ug", timeUnit="min") ; R1
R2 = sNCA(PK02$TIME, PK02$DV, dose=100, adm="Extravascular", doseUnit="ug", timeUnit="min") ; R2
BA = R2["AUCIFO"]/R1["AUCIFO"]; BA * 100 # Absolute Bioavailability (BA)


### Additive error models
source("D:/RP/wnl/wnl_0.1.0/R/CovStep.R")
source("D:/RP/wnl/wnl_0.1.0/R/EstStep.R")
source("D:/RP/wnl/wnl_0.1.0/R/InitStep.R")
source("D:/RP/wnl/wnl_0.1.0/R/MinUtil.R")
e = new.env()

Obj = function (vPara)
{
  if (e$STEP == "EST") vPara = DECN(vPara)
  e$THETA = vPara[1:e$nTheta]
  e$SG = diag(vPara[e$SGindex], nrow = e$nEps)
  FHi = e$PRED(e$THETA, e$Data)
  Ri = e$Data[, "DV"] - FHi[, "F"]
  Hi = FHi[, e$HNames, drop = FALSE]
  Ci = diag(diag(Hi %*% e$SG %*% t(Hi)), nrow = nrow(e$Data))
  return(determinant(Ci, logarithm = TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri)
}

## Model with tlag
nTheta = 4
nEps = 1

THETAinit = c(0.05, 0.1, 30, 20) # Initial estimate
SGinit = matrix(c(0.2), nrow=nEps, ncol=nEps) ; SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

Dose = 100

PRED = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  BA   = BA
  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]
  tlag = THETA[4]

  F  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag)))
  H1 = 1

  return(cbind(F, H1))
}
BA = 1 
colnames(PK02) = c("TIME", "DV", "DV2")
InitStep(PK02, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED)
EstRes = EstStep() ; EstRes ; EstRes$'Final Estimates'[3]/BA
CovRes = CovStep() ; CovRes

## p480   Fig 2.3
dev.new()
plot(PK02$TIME, PK02$DV2)
points(PK02$TIME, PK02$DV, col="red")

xp1 = 0:420
yp1 = PRED(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp1, type="l")

## Model without tlag
nTheta = 3
nEps = 1

THETAinit = c(0.05, 0.1, 30) # Initial estimate
SGinit = matrix(c(0.2), nrow=nEps, ncol=nEps) ; SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

Dose = 100

PRED1 = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  BA   = BA
  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]

  F  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*TIME) - exp(-Ka*TIME))
  H1 = 1

  return(cbind(F, H1))
}

colnames(PK02) = c("TIME", "DV", "DV2")
InitStep(PK02, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED1)
EstRes = EstStep() ; EstRes ; EstRes$'Final Estimates'[3]/BA
CovRes = CovStep() ; CovRes

yp2 = PRED1(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp2, type="l", col="red")

## Obj for simultaneous fitting
nTheta = 5
nEps = 1
THETAinit = c(0.05, 0.1, 30, 10, 0.6) # Initial estimate
SGinit = matrix(c(0.2), nrow=nEps, ncol=nEps) ; SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

Dose = 100

PRED0 = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]
  tlag = THETA[4]
  BA   = THETA[5]

  F  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag)))
  H1 = 1

  return(cbind(F, H1))
}


PRED2 = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  K    = THETA[1]
  V    = THETA[3]

  F  = DOSE/V*exp(-K*TIME)
  H1 = 1

  return(cbind(F, H1))
}

Obj = function (vPara)
{
  if (e$STEP == "EST") vPara = DECN(vPara)
  e$THETA = vPara[1:e$nTheta]
  e$SG = diag(vPara[e$SGindex], nrow = e$nEps)
  FHi = PRED0(e$THETA, e$Data)
  FHi = rbind(FHi, PRED2(e$THETA, e$Data))
  Ri = matrix(c(e$Data[, "DV"], PK02[,"DV2"]) - FHi[, "F"], ncol=1)
  Hi = FHi[, e$HNames, drop = FALSE]
  Ci = diag(diag(Hi %*% e$SG %*% t(Hi)), nrow = 2*nrow(e$Data))
  return(determinant(Ci, logarithm = TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri)
}

InitStep(PK02, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED1)
EstRes = EstStep() ; EstRes
CovRes = CovStep() ; CovRes

yp3 = PRED(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp3, type="l", col="blue")

yp4 = PRED2(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp4, type="l", col="green")



##########################################################################
### Combined error models
## Oral 1 Comp Model fitting

Obj = function (vPara)
{
  if (e$STEP == "EST") vPara = DECN(vPara)
  e$THETA = vPara[1:e$nTheta]
  e$SG = diag(vPara[e$SGindex], nrow = e$nEps)
  FHi = e$PRED(e$THETA, e$Data)
  Ri = e$Data[, "DV"] - FHi[, "F"]
  Hi = FHi[, e$HNames, drop = FALSE]
  Ci = diag(diag(Hi %*% e$SG %*% t(Hi)), nrow = nrow(e$Data))
  return(determinant(Ci, logarithm = TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri)
}

## Model with tlag
nTheta = 4
nEps = 2

THETAinit = c(0.05, 0.1, 20, 10) # Initial estimate
SGinit = matrix(c(0.2, 0, 0, 0.2), nrow=nEps, ncol=nEps) ; SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

Dose = 100

PRED = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  BA   = BA
  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]
  tlag = THETA[4]

  F  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag)))
  H1 = F
  H2 = 1

  return(cbind(F, H1, H2))
}

colnames(PK02) = c("TIME", "DV", "DV2")
InitStep(PK02, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED)
EstRes = EstStep() ; EstRes ; EstRes$'Final Estimates'[3]/BA
CovRes = CovStep() ; CovRes

## p480   Fig 2.3
dev.new()
plot(PK02$TIME, PK02$DV2)
points(PK02$TIME, PK02$DV, col="red")

xp1 = 0:420
yp1 = PRED(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp1, type="l")

## Model without tlag
nTheta = 3
nEps = 2

THETAinit = c(0.05, 0.1, 30) # Initial estimate
SGinit = matrix(c(0.2, 0, 0, 0.2), nrow=nEps, ncol=nEps) ; SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

Dose = 100

PRED1 = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  BA   = BA
  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]

  F  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*TIME) - exp(-Ka*TIME))
  H1 = F
  H2 = 1

  return(cbind(F, H1, H2))
}

colnames(PK02) = c("TIME", "DV", "DV2")
InitStep(PK02, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED1)
EstRes = EstStep() ; EstRes
CovRes = CovStep() ; CovRes

yp2 = PRED1(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp2, type="l", col="red")

## Obj for simultaneous fitting
nTheta = 5
nEps = 2
THETAinit = c(0.05, 0.1, 30, 10, 0.6) # Initial estimate
SGinit = matrix(c(0.2, 0, 0, 0.2), nrow=nEps, ncol=nEps) ; SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

Dose = 100

PRED0 = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]
  tlag = THETA[4]
  BA   = THETA[5]

  F  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag)))
  H1 = F
  H2 = 1

  return(cbind(F, H1, H2))
}

PRED2 = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  K    = THETA[1]
  V    = THETA[3]

  F  = DOSE/V*exp(-K*TIME)
  H1 = F
  H2 = 1

  return(cbind(F, H1, H2))
}

Obj = function (vPara)
{
  if (e$STEP == "EST") vPara = DECN(vPara)
  e$THETA = vPara[1:e$nTheta]
  e$SG = diag(vPara[e$SGindex], nrow = e$nEps)
  FHi = PRED0(e$THETA, e$Data)
  FHi = rbind(FHi, PRED2(e$THETA, e$Data))
  Ri = matrix(c(e$Data[, "DV"], PK02[,"DV2"]) - FHi[, "F"], ncol=1)
  Hi = FHi[, e$HNames, drop = FALSE]
  Ci = diag(diag(Hi %*% e$SG %*% t(Hi)), nrow = 2*nrow(e$Data))
  return(determinant(Ci, logarithm = TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri)
}

InitStep(PK02, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED1)
EstRes = EstStep() ; EstRes
CovRes = CovStep() ; CovRes

yp3 = PRED(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp3, type="l", col="blue")

yp4 = PRED2(EstRes$'Final Estimates', data.frame(TIME=xp1))[,1]
lines(xp1, yp4, type="l", col="green")
