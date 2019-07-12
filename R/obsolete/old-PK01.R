setwd("D:/Rt/Gab")
Data = read.csv("PK01.csv")
colnames(Data) = c("ID", "TIME", "DV")

## Plot
library(lattice)
dev.new()
xyplot(DV ~ TIME | as.factor(ID), data=Data, type="o")
dev.new()
xyplot(log(DV) ~ TIME | as.factor(ID), data=Data, type="o")


## NCA
library(NonCompart)
tblNCA(Data, key="ID", colTime="TIME", colConc="DV", dose=10, adm="Bolus")


## 1 Comp modeling, IV
library(wnl)
nTheta = 2
nEps = 1

THETAinit = c(20, 0.02) # Initial estimate
SGinit = matrix(c(1), nrow=nEps, ncol=nEps)
SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

Dose = 10000

PRED = function(THETA, DATAi) # Prediction function
{
  DOSE = Dose
  TIME = DATAi[,"TIME"]

  V  = THETA[1]
  K  = THETA[2]

  F  = DOSE/V * exp(-K*TIME)
  H1 = 1

  return(cbind(F, H1))
}

grCL = deriv(~V*k, c("V","k","sigsq"), function.arg=c("V", "k", "sigsq"), func=TRUE)
grAUC = deriv(~Dose/V/k, c("V", "k", "sigsq"), function.arg=c("V", "k", "sigsq"), func=TRUE)
grAUMC = deriv(~Dose/V/k/k, c("V", "k", "sigsq"), function.arg=c("V", "k", "sigsq"), func=TRUE)
grThalf = deriv(~log(2)/k, c("V", "k", "sigsq"), function.arg=c("V", "k", "sigsq"), func=TRUE)
grMRT = deriv(~1/k, c("V", "k", "sigsq"), function.arg=c("V", "k", "sigsq"), func=TRUE)

#########
IDs = unique(Data$ID)
nID = length(IDs)

Result = matrix(nrow=nID, ncol=27)
colnames(Result) = c("V_PE", "V_SE", "V_CV",
                     "k_PE", "k_SE", "k_CV",
                     "sgsq_PE", "sgsq_SE", "sgsq_CV",
                     "sg_PE", "sg_SE", "sg_CV",
                     "CL_PE", "CL_SE", "CL_CV",
                     "AUC_PE", "AUC_SE", "AUC_CV",
                     "AUMC_PE", "AUMC_SE", "AUMC_CV",
                     "THALF_PE", "THALF_SE", "THALF_CV",
                     "MRT_PE", "MRT_SE", "MRT_CV")

for (i in 1:nID) {
  Datai = Data[Data$ID == IDs[i], ]
  InitStep(Datai, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED)
  EstRes = EstStep()
  CovRes = CovStep()

  Result[i,1] = EstRes$'Final Estimates'[1]
  Result[i,2] = CovRes$'Standard Error'[1]
  Result[i,3] = Result[i,2] / Result[i,1] * 100

  Result[i,4] = EstRes$'Final Estimates'[2]
  Result[i,5] = CovRes$'Standard Error'[2]
  Result[i,6] = Result[i,5] / Result[i,4] * 100

  Result[i,7] = EstRes$'Final Estimates'[3]
  Result[i,8] = CovRes$'Standard Error'[3]
  Result[i,9] = Result[i,8] / Result[i,7] * 100

  Result[i,10] = sqrt(EstRes$'Final Estimates'[3])
  Result[i,11] = CovRes$'Standard Error'[3]/2/sqrt(EstRes$'Final Estimates'[3]) # See Wackerly p484  gr = (x^0.5)' = 0.5x^(-0.5)   gr^2 = 1/(4*x)
  Result[i,12] = Result[i,11] / Result[i,10] * 100

  gr1 = attr(grCL(Result[i,1], Result[i,4], Result[i,7]), "gradient")
  Result[i,13] = EstRes$'Final Estimates'[1]*EstRes$'Final Estimates'[2]
  Result[i,14] = sqrt(gr1 %*% CovRes$`Covariance Matrix of Estimates` %*% t(gr1))
  Result[i,15] = Result[i,14] / Result[i,13] * 100

  gr2 = attr(grAUC(Result[i,1], Result[i,4], Result[i,7]), "gradient")
  Result[i,16] = Dose/(EstRes$'Final Estimates'[1]*EstRes$'Final Estimates'[2])
  Result[i,17] = sqrt(gr2 %*% CovRes$`Covariance Matrix of Estimates` %*% t(gr2))
  Result[i,18] = Result[i,17] / Result[i,16] * 100

  gr3 = attr(grAUMC(Result[i,1], Result[i,4], Result[i,7]), "gradient")
  Result[i,19] = Dose/(EstRes$'Final Estimates'[1]*EstRes$'Final Estimates'[2]^2)
  Result[i,20] = sqrt(gr3 %*% CovRes$`Covariance Matrix of Estimates` %*% t(gr3))
  Result[i,21] = Result[i,20] / Result[i,19] * 100

  gr4 = attr(grThalf(Result[i,1], Result[i,4], Result[i,7]), "gradient")
  Result[i,22] = log(2)/EstRes$'Final Estimates'[2]
  Result[i,23] = sqrt(gr4 %*% CovRes$`Covariance Matrix of Estimates` %*% t(gr4))
  Result[i,24] = Result[i,23] / Result[i,22] * 100

  gr5 = attr(grMRT(Result[i,1], Result[i,4], Result[i,7]), "gradient")
  Result[i,25] = 1/EstRes$'Final Estimates'[2]
  Result[i,26] = sqrt(gr5 %*% CovRes$`Covariance Matrix of Estimates` %*% t(gr5))
  Result[i,27] = Result[i,26] / Result[i,25] * 100
} ; Result

Result[,c("V_PE", "V_CV", "k_PE", "k_CV", "sgsq_PE", "sgsq_CV", "sg_PE", "sg_CV",
          "CL_PE", "CL_CV", "AUC_PE", "AUC_CV", "AUMC_PE", "AUMC_CV", "THALF_PE", "THALF_CV", "MRT_PE", "MRT_CV")]


## Planar CI
ColNames = c("V_PE", "V_SE", "V_CV", "V_uLL", "V_uUL", "V_pLL", "V_pUL",
             "k_PE", "k_SE", "k_CV", "k_uLL", "k_uUL", "k_pLL", "k_pUL")
Res2 = matrix(nrow=nID, ncol=length(ColNames))
colnames(Res2) = ColNames

ci = 0.95
p = 1 - (1 - ci)/2

for (i in 1:nID) {
  Datai = Data[Data$ID == IDs[i], ]
  InitStep(Datai, THETAinit=THETAinit, SGinit=SGinit, nTheta=nTheta, LB=LB, UB=UB, Pred=PRED)
  EstRes = EstStep()
  CovRes = CovStep()

  Res2[i,"V_PE"] = EstRes$'Final Estimates'[1]
  Res2[i,"V_SE"] = CovRes$'Standard Error'[1]
  Res2[i,"V_CV"] = Res2[i,"V_SE"] / Res2[i,"V_PE"] * 100
  
  df = nrow(Datai) - 2
  Res2[i,"V_uLL"] = Res2[i,"V_PE"] - qt(p, df) * Res2[i,"V_SE"]
  Res2[i,"V_uUL"] = Res2[i,"V_PE"] + qt(p, df) * Res2[i,"V_SE"]
  
  Res2[i,"k_PE"] = EstRes$'Final Estimates'[2]
  Res2[i,"k_SE"] = CovRes$'Standard Error'[2]
  Res2[i,"k_CV"] = Res2[i,"k_SE"] / Res2[i,"k_PE"] * 100
  
  Res2[i,"k_uLL"] = Res2[i,"k_PE"] - qt(p, df) * Res2[i,"k_SE"]
  Res2[i,"k_uUL"] = Res2[i,"k_PE"] + qt(p, df) * Res2[i,"k_SE"]

  mu = EstRes$'Final Estimates'[1:2]
  mCov = CovRes$`Covariance Matrix of Estimates`[1:2, 1:2]
  eg = eigen(mCov) # sorted as descening order
  alpha = atan(eg$vectors[2,1]/eg$vectors[1,1])
  dimR = 2
  radius = sqrt(dimR*qf(ci, dimR, df)*eg$values) # plotting dimension * qf(ci, plotting dimension, npoints - npara)

  x0max = sqrt(radius[1]^2*cos(alpha)^2 + radius[2]^2*sin(alpha)^2) 
  y0max = sqrt(radius[1]^2*sin(alpha)^2 + radius[2]^2*cos(alpha)^2)
  xrange = mu[1] + c(-1,1)*x0max
  yrange = mu[2] + c(-1,1)*y0max
  
  Res2[i,"V_pLL"] = xrange[1]
  Res2[i,"V_pUL"] = xrange[2]
  
  Res2[i,"k_pLL"] = yrange[1]
  Res2[i,"k_pUL"] = yrange[2]
  
} ; Res2

## Simpler Objs

Obj = function(vPara) 
{
  if (e$STEP == "EST") vPara = DECN(vPara)
  e$THETA = vPara[1:e$nTheta]
  e$SG = diag(vPara[e$SGindex], nrow = e$nEps)
  FHi = e$PRED(e$THETA, e$Data)
  Ri = e$Data[, "DV"] - FHi[, "F"]
  Hi = FHi[, e$HNames, drop = FALSE]
  Ci = diag(diag(Hi %*% e$SG %*% t(Hi)), nrow = nrow(e$Data))
  return(t(Ri) %*% solve(Ci) %*% Ri)
}


Obj = function(vPara) 
{
  if (e$STEP == "EST") vPara = DECN(vPara)
  e$THETA = vPara[1:e$nTheta]
  e$SG = diag(vPara[e$SGindex], nrow = e$nEps)
  FHi = e$PRED(e$THETA, e$Data)
  Ri = e$Data[, "DV"] - FHi[, "F"]
  return(sum(Ri*Ri)/e$SG)
}

