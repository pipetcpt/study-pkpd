setwd("D:/Rt/Gab/")
require(wnl)
dPK48 = read.csv("PK48.csv", skip=1)
colnames(dPK48) = c("TIME", "CMT", "DV", "Analyte") ; dPK48

Bolus = 500 # umol

require(deSolve)
PKde =function(t, y, p)
{
  dy1dt = (-p["Vmax"]*y[1]/(p["Km"] + y[1]) - p["CLr"]*y[1])/p["V"] # eq 48:1
  dy2dt = p["CLr"]*y[1]                                             # eq 48:2
  dy3dt = p["Vmax"]*y[1]/(p["Km"] + y[1])                           # eq 48:3
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

TIME = sort(unique(c(0, dPK48[,"TIME"])))
y = lsoda(c(Bolus/24.5279, 0, 0), times=TIME, func=PKde, parms=c(V=24.528, Km=5.3, Vmax=51.4, CLr=2.46764))
y

# Figure 48.1, p 695
plot(0, 0.1, type="n", xlim=c(0, 16), ylim=c(0.1, 1000), log="y", xlab="Time (h)", ylab="Conc (uM) & Amount (umol)")
points(dPK48[dPK48$CMT==1,"TIME"], dPK48[dPK48$CMT==1,"DV"], pch=19)
points(dPK48[dPK48$CMT==2,"TIME"], dPK48[dPK48$CMT==2,"DV"], pch=17)
points(dPK48[dPK48$CMT==3,"TIME"], dPK48[dPK48$CMT==3,"DV"], pch=15)
lines(TIME, y[,"1"])
lines(TIME, y[,"2"])
lines(TIME, y[,"3"])
###


iTime1 = TIME %in% dPK48[dPK48$CMT==1,"TIME"]
iTime2 = TIME %in% dPK48[dPK48$CMT==2,"TIME"]
iTime3 = TIME %in% dPK48[dPK48$CMT==3,"TIME"]

fPK48 = function(THETA)
{
  V    = THETA[1]
  Km   = THETA[2]
  Vmax = THETA[3]
  CLr  = THETA[4]
  
  y = lsoda(c(Bolus/V, 0, 0), times=TIME, func=PKde, parms=c(V=V, Km=Km, Vmax=Vmax, CLr=CLr))
  return(c(y[iTime1, "1"], y[iTime2, "2"], y[iTime3, "3"]))
}

nlr(fPK48, dPK48, pNames=c("V", "Km", "Vmax", "CLr"), IE=c(25, 2, 12, 2.4), Error="P")

