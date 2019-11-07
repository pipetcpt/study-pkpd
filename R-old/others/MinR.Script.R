setwd("D:/Rt/Gab")
DataAll = read.csv("PK01.csv")
colnames(DataAll) = c("ID", "TIME", "DV")

IDs = unique(DataAll[,"ID"])
nID = length(IDs)

AddSD = min(DataAll[,"DV"])/4
PropSD = 0.2

i = 1

e = new.env()

e$Data = DataAll[DataAll$ID == IDs[i],]
e$nRec = nrow(e$Data)
e$pNames = c("V", "k")
e$IE = c(20, 0.2)
e$THETA = e$IE
e$nTheta = length(e$IE)
e$Error = "A" # "A" for Addtive error model
              # "P" for proportional error model
              # "C" for combined error model
e$AddSD = AddSD
e$PropSD = PropSD

if(e$Error == "A") {
  e$nEps = 1
  e$IE = c(e$IE, AddSD)
} else if (e$Error == "P") {
  e$nEps = 1
  e$IE = c(e$IE, PropSD)
} else if (e$Error == "C") {
  e$nEps = 2
  e$IE = c(e$IE, AddSD, PropSD)
}
e$nPara = e$nTheta + e$nEps
e$LB = rep(0, e$nPara)
e$UB = rep(10^6, e$nPara)

e$alpha  = 0.1 - log(e$IE/(e$UB - e$LB)/(1 - e$IE/(e$UB - e$LB)))
e$SGindex = (e$nTheta + 1):(e$nTheta + e$nEps)

PRED = function(THETA, DATAi) # Prediction function
{
  DOSE = 10000
  TIME = DATAi[,"TIME"]
  V  = THETA[1]
  K  = THETA[2]

  F  = DOSE/V * exp(-K*TIME)
  return(F)
}

OBJ = function(vPara)
{
  if (e$STEP=="EST") {
    vPara = exp(vPara - e$alpha)/(exp(vPara - e$alpha) + 1)*(e$UB - e$LB)
  }

  e$THETA = vPara[1:e$nTheta]

  Fi = e$Pred(e$THETA, e$Data)
  Ri = e$Data[,"DV"] - Fi

  if (e$Error == "A") {
    Hi = rep(vPara[e$SGindex], e$nRec)
  } else if (e$Error == "P") {
    Hi = vPara[e$SGindex]*Fi
  } else if (e$Error == "C") {
    Hi = rep(vPara[e$SGindex[1]]^2, e$nRec) + vPara[e$SGindex[2]]^2*Fi^2
  }

  Ci    = diag(diag(Hi %*% t(Hi)), nrow=e$nRec)
  return(determinant(Ci, logarithm=TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri)
}

e$Pred = PRED
e$Obj = OBJ

source("D:/RP/wnl/wnl_0.2.0/R/Hessian.R")

  e$STEP = "EST"
  Res = optim(rep(0.1, e$nPara), e$Obj, method="L-BFGS-B")
  e$FinalEst = exp(Res$par - e$alpha)/(exp(Res$par - e$alpha) + 1)*(e$UB - e$LB)
  EstRes = list(Res, e$FinalEst)
  names(EstRes) = list("Optim", "Final Estimates") ; EstRes

  e$STEP = "COV"
  InvCov = Hessian(e$Obj, e$FinalEst)/2  # FinalEst from EstStep()
  Cov = solve(InvCov)
  SE = sqrt(diag(Cov))
  Correl = cov2cor(Cov)
  EigenVal = eigen(Correl)$values

  CovRes = list(SE, Cov, Correl, InvCov, EigenVal)
  names(CovRes) = list("Standard Error", "Covariance Matrix of Estimates", "Correlation Matrix of Estimates", "Inverse Covariance Matrix of Estimates", "Eigen Values")
  CovRes


## invalid test end

#Data = system.file("tests", "5th_ed_data.xlsx", package="wnl")

require(xlsx)
Data = read.xlsx("D:/RP/wnl/wnl_0.2.0/tests/5th_ed_data.xlsx", "PK1")
Data


require(RODBC)
Conn = odbcConnectExcel2007("D:/RP/wnl/wnl_0.2.0/tests/5th_ed_data.xlsx")
Data = sqlFetch(Conn, "PK1") ; Data


odbcDriverConnect(paste("DRIVER=Microsoft Excel Driver (*.xls)", "DBQ=D:\bdr\hills.xlsx",
## invalid test end
######################################################################################

source("D:/RP/wnl/wnl_0.2.0/R/BasicUtil.R")
source("D:/RP/wnl/wnl_0.2.0/R/nlr.R")
source("D:/RP/wnl/wnl_0.2.0/R/Objs.R")
#source("D:/RP/wnl/wnl_0.2.0/R/PredPK01.R")
e = new.env()

## PK01
dPK01 = read.csv("D:/Rt/Gab/PK01.csv")
colnames(dPK01) = c("ID", "TIME", "DV")

fPK01 = function(THETA) # Prediction function
{
  DOSE = Dose
  TIME = e$DATA[,"TIME"]
  V  = THETA[1]
  K  = THETA[2]

  P  = DOSE/V * exp(-K*TIME)
  return(P)
}

Dose = 100000
IDs = unique(dPK01[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  Data = dPK01[dPK01$ID == IDs[i],]
  Res = nlr(fPK01, Data, pNames=c("V", "k"), IE=c(20, 0.2))
#  Res = nlr(fPK01, Data, pNames=c("V", "k"), IE=c(20, 0.2), Error="P")
#  Res = nlr(fPK01, Data, pNames=c("V", "k"), IE=c(20, 0.2), Error="C")
  print(Res)
}

## PK02
dPK02 = read.csv("D:/Rt/Gab/PK02.csv", skip=1)
colnames(dPK02) = c("TIME", "DV")

fPK02 = function(THETA) # Prediction function
{
  DOSE = Dose
  TIME = e$DATA[,"TIME"]

  BA   = BA
  Ka   = THETA[1]
  V    = THETA[2]
  K    = THETA[3]
  tlag = THETA[4]

  P  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag)))
  return(P)
}

Dose = 100
BA = 0.6100976 # from NCA   AUCinf.po / AUCinf.iv
nlr(fPK02, dPK02, pNames=c("ka", "V", "k", "tlag"), IE=c(0.1, 30, 0.05, 20))
#nlr(fPK02, dPK02, pNames=c("k", "ka", "V", "tlag"), IE=c(0.05, 0.1, 30, 20), Error="C?)



## Theoph
tData = Theoph
colnames(tData) = c("ID", "BWT", "DOSE", "TIME", "DV")

fPK = function(THETA) # Prediction function
{
  DOSE = Dose
  TIME = e$DATA[,"TIME"]

  K    = THETA[1]
  Ka   = THETA[2]
  V    = THETA[3]

  P  = DOSE/V*Ka/(Ka - K) * (exp(-K*TIME) - exp(-Ka*TIME))
  return(P)
}
Dose = 320000

IDs = unique(tData[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  Data = tData[tData$ID == IDs[i],]
  Res = nlr(fPK, Data, pNames=c("k", "ka", "V"), IE=c(0.1, 3, 500), Error="A")
  print(Res)
}




## Indometh

## no use
fPK02(c(0.05, 0.1, 30, 20))

Fx = fPK02
Data = dPK02
pNames = c("k", "ka", "V", "tlag")
IE = c(0.05, 0.1, 30, 20)
Error = "A"
Obj = OBJ


e$Fx = fPK02
e$DATA = dPK02
e$Y = dPK02[,"DV"]
vPara = c(0.05, 0.1, 30, 20, 0.1)
e$Obj(vPara)

