setwd("D:/Rt/Gab/")
require(wnl)
dPK52 = read.csv("PK52.csv", as.is=TRUE)
colnames(dPK52) = c("TIME", "DV", "ID")
dPK52 = dPK52[tolower(dPK52[,"DV"]) != "missing",]
dPK52[,"DV"] = as.numeric(dPK52[,"DV"]) ; dPK52

IDs = unique(dPK52[,"ID"])
nID = length(IDs)

AMT = 20000
TAU = 15/60
infRate = AMT/TAU # 80000

require(deSolve)
PKde = function(t, y, p)
{
  if (t < TAU) {
    RateIn = infRate
  } else {
    RateIn = 0
  }
  dy1dt = (RateIn - p["CL"]*y[1] - p["CLd"]*y[1] + p["CLd"]*y[2])/p["Vc"]
  dy2dt = (p["CLd"]*y[1] - p["CLd"]*y[2])/p["Vt"]  
  return(list(c(dy1dt, dy2dt)))
}

fPK52 = function(THETA)
{
  Vc =  THETA[1]
  Vt =  THETA[2]
  CL =  THETA[3]
  CLd = THETA[4]
  cTIME = c(0, e$DATA[,"TIME"])
  iTime = 2:length(cTIME)
  y = lsoda(y=c(0, 0), times=cTIME, func=PKde, parms=c(Vc=Vc, Vt=Vt, CL=CL, CLd=CLd))
  return(y[iTime,"1"])
}

Result = cbind(ID=IDs, Vc=0, Vt=0, CL=0, CLd=0, Vc.RSE=0, Vt.RSE=0, CL.RSE=0, CLd.RSE=0)
for (i in 1:nID) {
  cID = IDs[i]
  cData = dPK52[dPK52$ID == cID,]
  r = nlr(fPK52, cData, pNames=c("Vc", "Vt", "CL", "CLd"), IE=c(17, 7, 0.21, 0.1))
  Result[i, 2:5] = r$Est[1,1:e$nTheta]
  Result[i, 6:9] = r$Est[3,1:e$nTheta]
} ; Result

Result.Pois = cbind(ID=IDs, Vc=0, Vt=0, CL=0, CLd=0, Vc.RSE=0, Vt.RSE=0, CL.RSE=0, CLd.RSE=0)
for (i in 1:nID) {
  cID = IDs[i]
  cData = dPK52[dPK52$ID == cID,]
  r = nlr(fPK52, cData, pNames=c("Vc", "Vt", "CL", "CLd"), IE=c(17, 7, 0.21, 0.1), Error="POIS")
  Result.Pois[i, 2:5] = r$Est[1,1:e$nTheta]
  Result.Pois[i, 6:9] = r$Est[3,1:e$nTheta]
} ; Result.Pois

Result.P = cbind(ID=IDs, Vc=0, Vt=0, CL=0, CLd=0, Vc.RSE=0, Vt.RSE=0, CL.RSE=0, CLd.RSE=0)
for (i in 1:nID) {
  cID = IDs[i]
  cData = dPK52[dPK52$ID == cID,]
  r = nlr(fPK52, cData, pNames=c("Vc", "Vt", "CL", "CLd"), IE=c(17, 7, 0.21, 0.1), Error="P")
  Result.P[i, 2:5] = r$Est[1,1:e$nTheta]
  Result.P[i, 6:9] = r$Est[3,1:e$nTheta]
} ; Result.P # Table 52.1, p 719


# Figure 52.1, p 716
pRes = Result.P[,1:5]
plot(0, 0.01, type="n", xlim=c(0, 500), ylim=c(0.01, 10000), log="y", xlab="Time (min)", ylab="Concentration (ug/mL)")

y = vector()
for (i in 1:nID) {
  cID = IDs[i]
  cData = dPK52[dPK52$ID == cID,]  
  cTIME = c(0, dPK52[dPK52$ID == cID, "TIME"])
  iTime = 2:length(cTIME)  
  cy = lsoda(y=c(0, 0), times=cTIME, func=PKde, parms=pRes[i,])
  points(dPK52[dPK52$ID == cID,"TIME"], dPK52[dPK52$ID == cID,"DV"], pch=17-i)
  lines(cy[iTime,"time"], cy[iTime,"1"])    
  y = rbind(y, cy)
} ; y
###
