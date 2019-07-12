setwd("D:/Rt/Gab/")
require(wnl)
dPK46 = read.csv("PK46.csv", skip=1)
colnames(dPK46) = c("TIME", "DV") ; dPK46

AMT = 500000
TAU = 2016  # 2015.99 in WinNonlin
infRate = AMT/TAU ; infRate # 248.0159, while 248.017 in WinNonlin

require(deSolve)
PKde = function(t, y, p)
{
  if (t < 2016) {
    RateIn = infRate
  } else {
    RateIn = 0
  }
  
  dy1dt = (RateIn - p["CL"]*y[1])/p["V"]
  return(list(c(dy1dt)))
}

TIME = c(0, dPK46[,"TIME"])
iTime = 2:length(TIME)
y = lsoda(y=c(0), times=TIME, func=PKde, parms=c(V=35.5385, CL=61.2039))

## Figure 46.1, p 687
oldpar = par(no.readonly=TRUE)
dev.new()
par(mfrow=c(1,2))

plot(dPK46[,"TIME"], dPK46[,"DV"], xlim=c(0, 2500), ylim=c(0.1, 10), log="y", xlab="Time (h)", ylab="Concentration (nM)", pch=19)
lines(TIME[iTime], y[iTime, "1"])

plot(dPK46[,"TIME"], dPK46[,"DV"], xlim=c(2015, 2019), ylim=c(0.1, 10), log="y", xlab="Time (h)", ylab="Concentration (nM)", pch=19)
lines(TIME[iTime], y[iTime, "1"])

par(oldpar)
###


fPK46 = function(THETA)
{
  y = lsoda(y=c(0), times=TIME, func=PKde, parms=c(V=THETA[1], CL=THETA[2]))
  return(y[iTime, "1"])
}
fPK46(c(35.5385, 61.2039))

nlr(fPK46, dPK46, pNames=c("V", "CL"), IE=c(50, 100))

## using analytical equation, much faster
fPK46b = function(THETA)
{
  V   = THETA[1]
  CL  = THETA[2]
  T1 = dPK46[,"TIME"][dPK46[,"TIME"] <= TAU]
  T2 = dPK46[,"TIME"][dPK46[,"TIME"] > TAU]
  y1 = infRate/CL*(1 - exp(-CL/V*T1))
  y2 = infRate/CL*(1 - exp(-CL/V*TAU))*exp(-CL/V*(T2 - TAU))
  return(c(y1, y2))
}
fPK46b(c(35.5385, 61.2039))

nlr(fPK46b, dPK46, pNames=c("V", "CL"), IE=c(50, 100))
