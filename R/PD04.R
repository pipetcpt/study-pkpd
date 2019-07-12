# PD 4
require(wnl)
setwd("D:/Rt/PD")

dPD4 = read.csv("PD4.csv")
colnames(dPD4) = c("Time", "DV")

# Model 1
require(deSolve)
fPD4ade = function(t, y, p)
{
  Cw   = 1.0511*exp(-0.0228*t)
  Inh  = 1/(1 + (Cw/p["IC50"])^p["n"])
  dy1dt = p["Kout"]*(p["P0"]*Inh - y[1])
  return(list(c(dy1dt)))
}
TIMEs = dPD4[,"Time"]
lsoda(y=c(120.951), times=TIMEs, func=fPD4ade, parms=c(IC50=0.26, n=3, Kout=0.03, P0=120.951))

fPD4a = function(THETA)
{
 F1 =lsoda(y=THETA[4], times=TIMEs, func=fPD4ade, parms=c(IC50=THETA[1], n=THETA[2], Kout=THETA[3], P0=THETA[4]))
 return(F1[,"1"])
}
r1 = nlr(fPD4a, dPD4, pNames=c("IC50", "n", "Kout", "P0"), IE=c(0.35, 3.5, 0.3, 130)) ; r1

# Model 2
fPD4bde = function(t, y, p)
{
  if (t > p["Tlag"]) {
    Cw = 1.0511*exp(-0.0228*(t - p["Tlag"]))
  } else {
    Cw = 0
  }  
  Inh  = 1/(1 + (Cw/p["IC50"])^p["n"])
  dy1dt = p["Kout"]*(p["P0"]*Inh - y[1])
  return(list(c(dy1dt)))
}

fPD4b = function(THETA)
{
  F1 = lsoda(y=THETA[4], times=TIMEs, func=fPD4bde, parms=c(IC50=THETA[1], n=THETA[2], Kout=THETA[3], P0=THETA[4], Tlag=THETA[5]))
  return(F1[,"1"])
}
r2 = nlr(fPD4b, dPD4, pNames=c("IC50", "n", "Kout", "P0", "Tlag"), IE=c(0.35, 3.5, 0.3, 130, 1)) ; r2

# Figure 4.2
plot(dPD4[,"Time"], dPD4[,"DV"], xlim=c(0, 160), ylim=c(0, 140), xlab="Time (h)", ylab="Response (sec)", pch=16)
abline(h=seq(20, 140, 20), lty=3)
F2 = lsoda(y=r2$Est["PE","P0"], times=0:150, func=fPD4bde, parms=r2$Est["PE",c("IC50", "n", "Kout", "P0", "Tlag")])
lines(F2[,"time"], F2[,"1"])

# Figure 4.5
F3 = lsoda(y=c(130), times=0:150, func=fPD4ade, parms=c(IC50=0.35, n=3.5, Kout=0.3, P0=130))
F4 = lsoda(y=c(130), times=0:150, func=fPD4ade, parms=c(IC50=2*0.35, n=3.5, Kout=0.3, P0=130))
F5 = lsoda(y=c(130), times=0:150, func=fPD4ade, parms=c(IC50=0.5*0.35, n=3.5, Kout=0.3, P0=130))
dev.new()
plot(F3[,1], F3[,2], type="l", xlim=c(0, 160), ylim=c(0, 140), xlab="Time (h)", ylab="Response (sec)")
abline(h=seq(20, 140, 20), lty=3)
lines(F4[,1], F4[,2], lty=2)
lines(F5[,1], F5[,2], lty=2)

# Figure 4.6
F6 = lsoda(y=c(2*130), times=0:150, func=fPD4ade, parms=c(IC50=0.35, n=3.5, Kout=0.5*0.3, P0=2*130))
F7 = lsoda(y=c(130/2), times=0:150, func=fPD4ade, parms=c(IC50=0.35, n=3.5, Kout=2*0.3, P0=130/2))
dev.new()
plot(F3[,1], F3[,2], log="y", type="l", xlim=c(0, 160), ylim=c(1, 1000), xlab="Time (h)", ylab="Response (sec)")
abline(h=c(1:9, seq(10, 90, 10), seq(100, 1000, 100)), lty=3)
lines(F6[,1], F6[,2], lty=2)
lines(F7[,1], F7[,2], lty=2)

# Figure 4.7
F8 = lsoda(y=c(2*130), times=0:150, func=fPD4ade, parms=c(IC50=0.35, n=3.5, Kout=0.3, P0=2*130))
F9 = lsoda(y=c(130/2), times=0:150, func=fPD4ade, parms=c(IC50=0.35, n=3.5, Kout=0.3, P0=130/2))
dev.new()
plot(F3[,1], F3[,2], type="l", xlim=c(0, 160), ylim=c(0, 300), xlab="Time (h)", ylab="Response (sec)")
abline(h=seq(50, 300, 50), lty=3)
lines(F8[,1], F8[,2], lty=2)
lines(F9[,1], F9[,2], lty=2)

# Model 4
fPD4cde = function(t, y, p)
{
  Cw = 1.0511*exp(-0.0228*t)
  Inh  = 1/(1 + (Cw/p["IC50"])^p["n"])
  dy1dt = p["Kout"]*(p["P0"]*Inh - y[1])
  dy2dt = p["Kout"]*y[1] - p["Kout"]*y[2]
  return(list(c(dy1dt, dy2dt)))
}
lsoda(y=c(122.693, 122.693), times=TIMEs, func=fPD4cde, parms=c(IC50=0.27, n=1.476, Kout=0.09, P0=122.693))

fPD4c = function(THETA)
{
  F1 = lsoda(y=c(THETA[4], THETA[4]), times=TIMEs, func=fPD4cde, parms=c(IC50=THETA[1], n=THETA[2], Kout=THETA[3], P0=THETA[4]))
  return(F1[,"2"])
}
r3 = nlr(fPD4c, dPD4, pNames=c("IC50", "n", "Kout", "P0"), IE=c(0.35, 3.5, 0.3, 130)) ; r3
r4 = nlr(fPD4c, dPD4, pNames=c("IC50", "n", "Kout", "P0"), IE=c(0.35, 3.5, 0.3, 130), Error="POIS") ; r4
