setwd("D:/Rt/Gab/")
require(wnl)
dPK43 = read.csv("PK43.csv", skip=1)
colnames(dPK43) = c("TIME", "DV") ; dPK43

require(deSolve)
PKde = function(t, y, p)
{
  if (t < p["Tlag"]) {
    Ka2 = 0
  } else {
    Ka2 = p["Ka2"]
  }
  dy1dt = -p["Ka1"]*y[1]  # buccal
  dy2dt = -Ka2*y[2]  # mouth
  dy3dt = (p["Ka1"]*y[1] + Ka2*y[2])/p["V"] - p["Ke"]*y[3] # central
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

cTIME = c(0, dPK43$TIME)
iTime = 2:length(cTIME)

# Figure 43.1, p 671
plot(dPK43[,"TIME"], dPK43[,"DV"], xlim=c(0, 25), ylim=c(0, 80), xlab="Time (h)", ylab="Concentration (ug/L)", pch=19)


AMT = 2000
Fr = 0.515023
y = lsoda(y=c(AMT*Fr, AMT*(1 - Fr), 0), times=cTIME, func=PKde, parms=c(V=20.6274, Ke=0.0887, Ka1=7.6237, Ka2=1.0751, Tlag=2.29614))
lines(cTIME, y[,"3"])
###

fPK43 = function(THETA)
{
  V    = THETA[1]
  Ke   = THETA[2]
  Ka1  = THETA[3]
  Ka2  = THETA[4]
  Fr   = THETA[5]
  Tlag = THETA[6]
  
  y = lsoda(y=c(AMT*Fr, AMT*(1 - Fr), 0), times=cTIME, func=PKde, parms=c(V=V, Ke=Ke, Ka1=Ka1, Ka2=Ka2, Tlag=Tlag))
  return(y[iTime,"3"])
}
fPK43(c(20.6274, 0.0887, 7.6237, 1.0751, 0.515023, 2.29614))

nlr(fPK43, dPK43, pNames=c("V", "Ke", "Ka1", "Ka2", "Fract", "Tlag"), IE=c(20, 0.1, 1, 0.5, 0.7, 2))
