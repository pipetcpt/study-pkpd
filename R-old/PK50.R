setwd("D:/Rt/Gab/")
require(wnl)
dPK50 = read.csv("PK50.csv", skip=1, as.is=TRUE)
colnames(dPK50) = c("TIME", "ID", "Response", "DV") # dPK50
dPK50 = dPK50[dPK50[,"DV"] != "Missing",] # dPK50
dPK50[,"TIME"] = as.numeric(dPK50[,"TIME"]) 
dPK50[,"Response"] = as.numeric(dPK50[,"Response"]) 
dPK50[,"DV"] = as.numeric(dPK50[,"DV"]) ; dPK50

IDs = unique(dPK50[,"ID"])
nID = length(IDs)

AMT = 566
TAU = 5
infRate = AMT/TAU # 113.2

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

## Naive model
fPK50a = function(THETA)
{
  Vc  = THETA[1]
  CL  = THETA[2]
  Vt  = THETA[3]
  CLd = THETA[4]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cTIME = c(0, dPK50[dPK50$ID == cID, "TIME"])
    iTime = 2:length(cTIME)
    cy = lsoda(y=c(0,0), times=cTIME, func=PKde, parms=c(Vc=Vc, CL=CL, Vt=Vt, CLd=CLd))
    y = c(y, cy[iTime, "1"])
  }
  return(y)
}
y = fPK50a(c(19.4481, 15.7538, 11.7851, 3.68409))
length(y)
length(dPK50[,"DV"])

nlr(fPK50a, dPK50, pNames=c("Vc", "CL", "Vt", "CLd"), IE=c(10, 30, 10, 10), Error="POIS") # Note Vt at Table 50.1

## Individual fitting
fPK50b = function(THETA)
{
  Vc  = THETA[1]
  CL  = THETA[2]
  Vt  = THETA[3]
  CLd = THETA[4]
  
  cTIME = c(0, e$DATA[, "TIME"])
  iTime = 2:length(cTIME)
  y = lsoda(y=c(0,0), times=cTIME, func=PKde, parms=c(Vc=Vc, CL=CL, Vt=Vt, CLd=CLd))
  return(y[iTime, "1"])
}

Result = cbind(ID=IDs, Vc=0, CL=0, Vt=0, CLd=0, Vc.RSE=0, CL.RSE=0, Vt.RSE=0, CLd.RSE=0)
for (i in 1:nID) {
  cID = IDs[i]
  cData = dPK50[dPK50$ID == cID,]
  Res = nlr(fPK50b, cData, pNames=c("Vc", "CL", "Vt", "CLd"), IE=c(10, 30, 10, 10), Error="POIS")
  Result[i,2:5] = Res$Est[1,1:e$nTheta]
  Result[i,6:9] = Res$Est[3,1:e$nTheta]
} ; Result # Table 50.3 (p710) is quite different from the WinNonlin result

Result.Add = cbind(ID=IDs, Vc=0, CL=0, Vt=0, CLd=0, Vc.RSE=0, CL.RSE=0, Vt.RSE=0, CLd.RSE=0)
for (i in 1:nID) {
  cID = IDs[i]
  cData = dPK50[dPK50$ID == cID,]
  Res = nlr(fPK50b, cData, pNames=c("Vc", "CL", "Vt", "CLd"), IE=c(10, 30, 10, 10))
  Result.Add[i,2:5] = Res$Est[1,1:e$nTheta]
  Result.Add[i,6:9] = Res$Est[3,1:e$nTheta]
} ; Result.Add # Table 50.3 (p710) is quite different from the WinNonlin result

## Emax model
dPK50b = dPK50[!is.na(dPK50[,"Response"]),] ; dPK50b
colnames(dPK50b) = c("TIME", "ID", "DV", "Cp") ; dPK50b

fPK50c = function(THETA)
{
  EC50 = THETA[1]
  Emax = THETA[2]
  n    = THETA[3]
  y = Emax*dPK50b[,"Cp"]^n/(EC50^n + dPK50b[,"Cp"]^n)
  return(y)
}

nlr(fPK50c, dPK50b, pNames=c("EC50", "Emax", "n"), IE=c(2, 3, 3)) # Table 50.2 has different result

## Figure 50.1 upper only (There is no unbound concentration data available.)
plot(0, 0.1, type="n", log="y", xlim=c(0, 25), ylim=c(0.1, 100), xlab="Time (h)", ylab="Concentration (ug/L)")
for (i in 1:nID) {
  cID = IDs[i]
  x = dPK50[dPK50$ID == cID, "TIME"]
  y = dPK50[dPK50$ID == cID, "DV"]
  lines(x, y, col="red")
}
###

## Figure 50.2 left only (There is no unbound concentration data available.)
plot(0.1, 0, type="n", log="x", xlim=c(0.1, 20), ylim=c(0, 4), xlab="Concentration (ug/L", ylab="Response")
for (i in 1:nID) {
  cID = IDs[i]
  x = dPK50b[dPK50b$ID == cID, "Cp"]
  y = dPK50b[dPK50b$ID == cID, "DV"]
  lines(x, y, col="red")
}
###

