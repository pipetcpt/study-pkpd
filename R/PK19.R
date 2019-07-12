setwd("D:/Rt/Gab/")
require(wnl)
dPK19 = read.csv("PK19.csv")
colnames(dPK19) = c("TIME","DV","ID","AMT","DV2") ; dPK19

require(deSolve)
PKde = function(t, y, p)
{
  if (t < TAU) {
    RateIn = cInf
  } else {
    RateIn = 0
  }
  dy1dt = (RateIn - p["Vmax"]*y[1]/(p["Km"] + y[1]) - p["Cld"]*y[1] + p["Cld"]*y[2])/p["Vc"] # eq 19:1
  dy2dt = (p["Cld"]*y[1] - p["Cld"]*y[2])/p["Vt"]                   # eq 19:2, parent peripheral
  dy3dt = p["Vmax"]*y[1]/(p["Km"] + y[1])/p["Vm"] - p["Kme"]*y[3] # eq 19:5, metabolite
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

TAU = 5
plot(0, 0.1, type="n", xlim=c(0, 300), ylim=c(0.1, 1000), xlab="Time (min)", ylab="Concentration (uM)", log="y")

IDs = unique(dPK19[,"ID"])
nID = length(IDs)
y = vector()
for (i in 1:nID) {
  cID = IDs[i]
  DATAi = dPK19[dPK19$ID == cID,]
  cAMT = DATAi[1, "AMT"]
  cInf = cAMT / TAU
  Times = DATAi[,"TIME"]
  y = rbind(y, lsoda(y=c(0, 0, 0), times=Times, func=PKde, parms=c(Vc=1.064, Vmax=1.64429, Km=54.794, Vt=2, Cld=0.1288, Kme=0.14516, Vm=0.29)))
  points(DATAi[,"TIME"], DATAi[,"DV"])
  points(DATAi[,"TIME"], DATAi[,"DV2"], pch=2, col=2)
  lines(y[,"time"], y[,"1"])
  lines(y[,"time"], y[,"3"], lty=2)
} ; y # Figure 19.1, p564

iTime = !is.na(dPK19[,"DV"])

fPK19 = function(THETA)
{
  Vc   = THETA[1]
  Vmax = THETA[2]
  Km   = THETA[3]
  Vt   = THETA[4]
  Cld  = THETA[5]
  Kme  = THETA[6]
  Vm   = THETA[7] # Vd of metabolite

  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    DATAi = dPK19[dPK19$ID == cID,]
    cAMT = DATAi[1, "AMT"]
    cInf <<- cAMT / TAU  # referencing within PKde
    Times = DATAi[,"TIME"]
    y = rbind(y, lsoda(y=c(0, 0, 0), times=Times, func=PKde, parms=c(Vc=Vc, Vmax=Vmax, Km=Km, Vt=Vt, Cld=Cld, Kme=Kme, Vm=Vm)))
  }

  return(c(y[iTime,"1"], y[iTime,"3"]))
}
fPK19(c(1.064, 1.64429, 54.794, 2, 0.1288, 0.14516, 0.29))

dPK19b = rbind(dPK19[iTime, c("ID","TIME","DV")], cbind(dPK19[iTime, c("ID","TIME")], DV=dPK19[iTime, "DV2"])) ; dPK19b[,"DV"]

nlr(fPK19, dPK19b, pNames=c("Vc", "Vmax", "Km", "Vt", "Cld", "Kme", "Vm"),
    IE=c(1, 1.5, 50, 2, 0.15, 0.15, 0.3),
    LB=c(0.5, 1, 25, 1, 0.1, 0.1, 0.1),
    UB=c(2, 2, 100, 3, 0.2, 0.2, 0.5))

nlr(fPK19, dPK19b, pNames=c("Vc", "Vmax", "Km", "Vt", "Cld", "Kme", "Vm"),
    IE=c(1, 1.5, 50, 2, 0.15, 0.15, 0.3),
    LB=c(0.5, 1, 25, 1, 0.1, 0.1, 0.1),
    UB=c(2, 2, 100, 3, 0.2, 0.2, 0.5), Error="POIS")

nlr(fPK19, dPK19b, pNames=c("Vc", "Vmax", "Km", "Vt", "Cld", "Kme", "Vm"),
    IE=c(1, 1.5, 50, 2, 0.15, 0.15, 0.3),
    LB=c(0.5, 1, 25, 1, 0.1, 0.1, 0.1),
    UB=c(2, 2, 100, 3, 0.2, 0.2, 0.5), Error="C")




## Parent alone (reduced model, Result at Table 19.2 p 569
PKde2 = function(t, y, p)
{
  if (t < TAU) {
    RateIn = cInf
  } else {
    RateIn = 0
  }
  dy1dt = (RateIn - p["Vmax"]*y[1]/(p["Km"] + y[1]) - p["Cld"]*y[1] + p["Cld"]*y[2])/p["Vc"] # eq 19:1
  dy2dt = (p["Cld"]*y[1] - p["Cld"]*y[2])/p["Vt"]                   # eq 19:2, parent peripheral
  return(list(c(dy1dt, dy2dt)))
}

fPK19b = function(THETA)
{
  Vc   = THETA[1]
  Vmax = THETA[2]
  Km   = THETA[3]
  Vt   = THETA[4]
  Cld  = THETA[5]

  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    DATAi = dPK19[dPK19$ID == cID,]
    cAMT = DATAi[1, "AMT"]
    cInf <<- cAMT / TAU
    Times = DATAi[,"TIME"]
    y = rbind(y, lsoda(y=c(0, 0), times=Times, func=PKde2, parms=c(Vc=Vc, Vmax=Vmax, Km=Km, Vt=Vt, Cld=Cld)))
  }

  return(y[iTime,"1"])
}

dPK19c = dPK19[iTime, c("ID","TIME","DV")]

nlr(fPK19b, dPK19c, pNames=c("Vc", "Vmax", "Km", "Vt", "Cld"), IE=c(1, 1.5, 50, 2, 0.15), Error="POIS")

