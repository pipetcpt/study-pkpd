setwd("D:/Rt/Gab/")
require(wnl)
dPK41 = read.csv("PK41.csv", skip=1, as.is=TRUE)
colnames(dPK41) = c("TIME", "DV", "ID") ; dPK41
dPK41[,"DV"] = as.numeric(dPK41[,"DV"]) ; dPK41

AMTs = c(310, 520, 780)

require(deSolve)
PKde = function(t, y, p)
{
  if (t < 5) {
    RateIn = AMTs[cID]/5
  } else {
    RateIn = 0
  }
  dy1dt = (RateIn - p["Vmax"]*y[1]/(p["Km"] + y[1]))/p["V"]
  return(list(c(dy1dt))) 
}

# Figure 41.1, p 662
plot(0, 0.3, type="n", xlim=c(0, 10), ylim=c(0.3, 300), log="y", xlab="Time (h)", ylab="Concentration (ug/L)")
Bullets = c(19, 15, 17)

IDs = unique(dPK41[,"ID"])
nID = length(IDs)
y = vector()
for (i in 1:nID) {
  cID <<- IDs[i] # referencing within PKde
  cTIME = c(0, dPK41[dPK41$ID==cID,"TIME"])
  iTime = 2:length(cTIME)
  cy = lsoda(y=c(0), times=cTIME, func=PKde, parms=c(V=1.8, Km=79.8382, Vmax=180.311))
  points(dPK41[dPK41$ID==cID,"TIME"], dPK41[dPK41$ID==cID,"DV"], pch=Bullets[i])
  lines(cTIME[iTime], cy[iTime,"1"])
  y = c(y, cy[iTime,"1"])
} ; y
###

dPK41 = dPK41[!is.na(dPK41$DV),] # removing NA points for fitting
fPK41 = function(THETA)
{
  V    = THETA[1]
  Km   = THETA[2]
  Vmax = THETA[3]

  y = vector()
  for (i in 1:nID) {
    cID <<- IDs[i] # referencing within PKde
    cTIME = c(0, dPK41[dPK41$ID==cID,"TIME"])
    iTime = 2:length(cTIME)
    cy = lsoda(y=c(0), times=cTIME, func=PKde, parms=c(V=V, Km=Km, Vmax=Vmax))
    y = c(y, cy[iTime,"1"])
  }
  return(y)
}
fPK41(c(1.8, 79.8382, 180.311))

nlr(fPK41, dPK41, pNames=c("V", "Km", "Vmax"), IE=c(2, 200, 100), Error="P")
