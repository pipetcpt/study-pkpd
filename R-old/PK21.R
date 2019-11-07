setwd("D:/Rt/Gab/")
require(wnl)
dPK21 = read.csv("PK21.csv", skip=2)
colnames(dPK21) = c("TIME", "AMT", "ADDL", "DV") ; dPK21

DKOR = 3796.94
TBP = 216
TBP2 = 512

require(deSolve)
PKde = function(t, y, p)
{
  RKPRE = p["CLpre"]/p["V"]
  RKSS = p["CLss"]/p["V"]

  Ke = RKPRE
  if (t >= TBP & t < TBP2) {
    Ke = RKSS - (RKSS - RKPRE)*exp(-p["RKD"]*(t - TBP))
  } else if (t > TBP2) {
    AS = RKSS - (RKSS - RKPRE)*exp(-p["RKD"]*(TBP2 - TBP))
    Ke = RKPRE - (RKPRE - AS)*exp(-p["RKD"]*(t - TBP2))  # erratum a -> A
  }

  dy1dt = -p["Ka"]*y[1]
  dy2dt = DKOR*p["Ka"]*y[1]/p["V"] - Ke*y[2]

  return(list(c(dy1dt, dy2dt)))
}

Tlag = 0.814121
dTIME = seq(0, 696, by=8) ; dTIME
nDose = length(dTIME) ; nDose
EventDat = data.frame(var = rep("y1", nDose),
                      time = (dTIME + Tlag),
                      value = rep(10, nDose),
                      method = rep("add", nDose))

gTimes = seq(0, max(dPK21[,"TIME"]), by=0.1)
y = lsoda(y=c(y1=0, y2=0), times=gTimes, func=PKde, events=list(data=EventDat),
          parms=c(Ka=1.8406, CLss=114.344, CLpre=46.296, RKD=0.00547243, V=1679.4))
y

# Figure 21.1 p 576
plot(dPK21[,"TIME"], dPK21[,"DV"], xlim=c(100, 800), ylim=c(0, 120), xlab="Time (h)", ylab="Concentration (ug/L)")
lines(y[,"time"], y[,"y2"]) # Figure 18.2 p558


Times = c(0, dPK21[, "TIME"])
iTime = 2:length(Times)

fPK21 = function(THETA)
{
  Ka    = THETA[1]
  CLss  = THETA[2]
  CLpre = THETA[3]
  RKD   = THETA[4]
  V     = THETA[5]
  Tlag  = THETA[6]

  Times = c(0, dPK21[, "TIME"])

  EventDat = data.frame(var = rep("y1", nDose),
                        time = (dTIME + Tlag),
                        value = rep(10, nDose),
                        method = rep("add", nDose))

  y = lsoda(y=c(y1=0, y2=0), times=Times, func=PKde, events=list(data=EventDat),
            parms=c(Ka=Ka, CLss=CLss, CLpre=CLpre, RKD=RKD, V=V), hmax=8)

  return(y[iTime,"y2"])
}
fPK21(c(1.8406, 114.344, 46.296, 0.00547243, 1679.4, 0.814121))

nlr(fPK21, dPK21, pNames=c("Ka", "CLss", "CLpre", "RKD", "V", "Tlag"),
    IE=c(3, 140, 40, 0.005, 1260, 0.7),
    LB=c(0.1, 20, 4, 0, 100, 0),
    UB=c(5, 500, 200, 1, 3000, 3)) # fitting failure -> See NONMEM output

