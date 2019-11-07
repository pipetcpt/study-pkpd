setwd("D:/Rt/Gab/")
require(wnl)
dPK45 = read.csv("PK45.csv", as.is=TRUE)
colnames(dPK45) = c("TIME", "MOL", "DV", "ID") # dPK45
dPK45 = dPK45[dPK45[,"DV"] != "Missing",] # dPK45
dPK45[,"DV"] = as.numeric(dPK45[,"DV"]) # dPK45
dPK45 = dPK45[order(dPK45$ID, dPK45$MOL, dPK45$TIME),] ; dPK45

## Infusion rate different between the book and WinNonlin data
# Book                                  <-> WinNonlin 
#   5 umol/kg, 15 sec = 1200 umol/kg/hr <->   5 umol/kg, 1190.476 umol/kg/hr (15.2 seconds)
# 4.3 umol/kg, 15 sec = 1032 umol/kg/hr <-> 4.3 umol/kg, 1023.81 umol/kg/hr (15.2 seconds)
# We use WinNonlin value  for the comparison
infRate = cbind(ID=c(1, 2), RATE1=c(0, 1023.81), RATE2=c(1190.476, 0)) ; infRate

require(deSolve)
PKde = function(t, y, p)
{
  if (t < 15.2/3600) {
    Rate1 = infRate[infRate[,"ID"] == cID, "RATE1"]
    Rate2 = infRate[infRate[,"ID"] == cID, "RATE2"]
  } else {
    Rate1 = 0
    Rate2 = 0 
  }

  dy1dt = (Rate1 - p["CL1"]*y[1] - p["CLd1"]*y[1] + p["CLd1"]*y[2] - p["CL12"]*y[1] + p["CL21"]*y[3])/p["Vc1"]
  dy2dt = (p["CLd1"]*y[1] - p["CLd1"]*y[2])/p["Vt1"]
  dy3dt = (Rate2 - p["CL2"]*y[3] - p["CLd2"]*y[3] + p["CLd2"]*y[4] + p["CL12"]*y[1] - p["CL21"]*y[3])/p["Vc2"]
  dy4dt = (p["CLd2"]*y[3] - p["CLd2"]*y[4])/p["Vt2"]
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt)))
}

## Figure 45.1, p 678
plot(0, 0001, type="n", xlim=c(0, 35), ylim=c(0.001, 10), log="y", xlab="Time (h)", ylab="Concentration (uM)")

IDs = unique(dPK45[,"ID"])
nID = length(IDs)

TIME = sort(unique(c(0, dPK45[,"TIME"]))) ; TIME

y = vector()
for (i in 1:nID) {
  cID <<- IDs[i] # referencing within PKde
  cy = lsoda(y=c(0, 0, 0, 0), times=TIME, func=PKde, parms=c(Vc1=0.6738, CL1=0.35, Vt1=0.291, CLd1=0.0618435, CL12=0.0146138, CL21=0.04441, Vc2=0.258623, CL2=0.0686, Vt2=0.722415, CLd2=16.8))
  iTime1 = TIME %in% dPK45[dPK45$ID==cID & dPK45$MOL == "Cp", "TIME"]
  iTime2 = TIME %in% dPK45[dPK45$ID==cID & dPK45$MOL == "Cm", "TIME"]
  points(dPK45[dPK45$ID==cID & dPK45$MOL == "Cp", "TIME"], dPK45[dPK45$ID==cID & dPK45$MOL == "Cp", "DV"], pch=19, col=i)
  points(dPK45[dPK45$ID==cID & dPK45$MOL == "Cm", "TIME"], dPK45[dPK45$ID==cID & dPK45$MOL == "Cm", "DV"], pch=15, col=i)
  lines(TIME, cy[,"1"], col=i)
  lines(TIME, cy[,"3"], lty=2, col=i)
  y = c(y, cy[iTime2,"3"], cy[iTime1,"1"]) # This is the order of arranged dPK45
} ; y
###

fPK45 = function(THETA)
{
  Vc1  = THETA[1]
  CL1  = THETA[2]
  Vt1  = THETA[3]
  CLd1 = THETA[4]
  CL12 = THETA[5]
  CL21 = THETA[6]
  Vc2  = THETA[7]
  CL2  = THETA[8]
  Vt2  = THETA[9]
  CLd2 = THETA[10]

  y = vector()
  for (i in 1:nID) {
    cID <<- IDs[i] # referencing within PKde
    cy = lsoda(y=c(0, 0, 0, 0), times=TIME, func=PKde, parms=c(Vc1=Vc1, CL1=CL1, Vt1=Vt1, CLd1=CLd1, CL12=CL12, CL21=CL21, Vc2=Vc2, CL2=CL2, Vt2=Vt2, CLd2=CLd2))
    iTime1 = TIME %in% dPK45[dPK45$ID==cID & dPK45$MOL == "Cp", "TIME"]
    iTime2 = TIME %in% dPK45[dPK45$ID==cID & dPK45$MOL == "Cm", "TIME"]
    y = c(y, cy[iTime2,"3"], cy[iTime1,"1"]) # This is the order of arranged dPK45
  }
  
  return(y)
}
fPK45(c(0.6738, 0.35, 0.291, 0.0618435, 0.0146138, 0.04441, 0.258623, 0.0686, 0.722415, 16.8))

nlr(fPK45, dPK45, pNames=c("Vc1", "CL1", "Vt1", "CLd1", "CL12", "CL21", "Vc2", "CL2", "Vt2", "CLd2"), IE=c(1, 0.5, 1, 0.5, 1, 0.01, 1, 0.01, 1, 10), Error="P")



## Reduced model
PKde2 = function(t, y, p)
{
  if (t < 15.2/3600) {
    Rate1 = infRate[infRate[,"ID"] == cID, "RATE1"]
    Rate2 = infRate[infRate[,"ID"] == cID, "RATE2"]
  } else {
    Rate1 = 0
    Rate2 = 0 
  }

  dy1dt = (Rate1 - p["CL1"]*y[1] - p["CLd1"]*y[1] + p["CLd1"]*y[2] - p["CL12"]*y[1] + p["CL21"]*y[3])/p["Vc1"]
  dy2dt = (p["CLd1"]*y[1] - p["CLd1"]*y[2])/p["Vt1"]
  dy3dt = (Rate2 - p["CL2"]*y[3] + p["CL12"]*y[1] - p["CL21"]*y[3])/p["Vc2"]
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

fPK45b = function(THETA)
{
  Vc1  = THETA[1]
  CL1  = THETA[2]
  Vt1  = THETA[3]
  CLd1 = THETA[4]
  CL12 = THETA[5]
  CL21 = THETA[6]
  Vc2  = THETA[7]
  CL2  = THETA[8]

  y = vector()
  for (i in 1:nID) {
    cID <<- IDs[i] # referencing within PKde
    cy = lsoda(y=c(0, 0, 0), times=TIME, func=PKde2, parms=c(Vc1=Vc1, CL1=CL1, Vt1=Vt1, CLd1=CLd1, CL12=CL12, CL21=CL21, Vc2=Vc2, CL2=CL2))
    iTime1 = TIME %in% dPK45[dPK45$ID==cID & dPK45$MOL == "Cp", "TIME"]
    iTime2 = TIME %in% dPK45[dPK45$ID==cID & dPK45$MOL == "Cm", "TIME"]
    y = c(y, cy[iTime2,"3"], cy[iTime1,"1"]) # This is the order of arranged dPK45
  }
  
  return(y)
}

nlr(fPK45b, dPK45, pNames=c("Vc1", "CL1", "Vt1", "CLd1", "CL12", "CL21", "Vc2", "CL2"), IE=c(1, 0.5, 1, 0.5, 1, 0.01, 1, 0.01), Error="P")
