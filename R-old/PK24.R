setwd("D:/Rt/Gab/")
require(wnl)
dPK24 = read.csv("PK24.csv", skip=1)
colnames(dPK24) = c("TIME", "DV") ; dPK24

require(deSolve)
PKde1 = function(t, y, p)
{
  RateIn = 0
  if (t < 2) RateIn = 5000

  dy1dt = (RateIn - p["CL"]*y[1] - p["CLd1"]*y[1] + p["CLd1"]*y[2] - p["CLd2"]*y[1] + p["CLd2"]*y[3])/p["Vc"] # eq 24:1
  dy2dt = (p["CLd1"]*y[1] - p["CLd1"]*y[2])/p["Vt1"] # eq 24:2
  dy3dt = (p["CLd2"]*y[1] - p["CLd2"]*y[3])/p["Vt2"] # eq 24:3

  return(list(c(dy1dt, dy2dt, dy3dt)))
}

PKde2 = function(t, y, p)
{
  RateIn = 0
  if (t < 2) RateIn = 5000

  CL = p["CLo"] - p["a"]*y[1] # eq 24:4
  dy1dt = (RateIn - CL*y[1] - p["CLd1"]*y[1] + p["CLd1"]*y[2] - p["CLd2"]*y[1] + p["CLd2"]*y[3])/p["Vc"] # eq 24:1
  dy2dt = (p["CLd1"]*y[1] - p["CLd1"]*y[2])/p["Vt1"] # eq 24:2
  dy3dt = (p["CLd2"]*y[1] - p["CLd2"]*y[3])/p["Vt2"] # eq 24:3

  return(list(c(dy1dt, dy2dt, dy3dt)))
}

# Figure 24.1, p 591
gTIME = seq(0, 7, by=0.1)
y1 = lsoda(y=c(0, 0, 0), times=gTIME, func=PKde1, parms=c(Vc=0.762136, Vt1=1.525, Vt2=2.38, CL=5.1, CLd1=6.638, CLd2=0.725))
y2 = lsoda(y=c(0, 0, 0), times=gTIME, func=PKde2, parms=c(Vc=0.684231, Vt1=1.772, Vt2=3.2, CLo=6.624, a=0.0025, CLd1=5.9, CLd2=0.93))
plot(dPK24[,"TIME"], dPK24[,"DV"], ylim=c(0, 1000), xlab="Time (h)", ylab="Concentration (ug/L)", pch=19)
lines(y1[,"time"], y1[,"1"], lty=2)
lines(y2[,"time"], y2[,"1"])
##

Times = c(0, dPK24[, "TIME"])
iTime = 2:length(Times)

fPK24a = function(THETA)
{
  Vc   = THETA[1]
  Vt1  = THETA[2]
  Vt2  = THETA[3]
  CL   = THETA[4]
  CLd1 = THETA[5]
  CLd2 = THETA[6]

  y = lsoda(y=c(0, 0, 0), times=Times, func=PKde1, parms=c(Vc=Vc, Vt1=Vt1, Vt2=Vt2, CL=CL, CLd1=CLd1, CLd2=CLd2))
  return(y[iTime, "1"])
}
fPK24a(c(0.762136, 1.525, 2.38, 5.1, 6.638, 0.725))

nlr(fPK24a, dPK24, pNames=c("Vc", "Vt1", "Vt2", "CL", "CLd1", "CLd2"), IE=c(1, 2, 4, 7, 6, 1), Error="P")


fPK24b = function(THETA)
{
  Vc   = THETA[1]
  Vt1  = THETA[2]
  Vt2  = THETA[3]
  CLo  = THETA[4]
  CLd1 = THETA[5]
  CLd2 = THETA[6]
  a    = THETA[7]

  y = lsoda(y=c(0, 0, 0), times=Times, func=PKde2, parms=c(Vc=Vc, Vt1=Vt1, Vt2=Vt2, CLo=CLo, CLd1=CLd1, CLd2=CLd2, a=a))
  return(y[iTime, "1"])
}
fPK24b(c(0.684231, 1.772, 3.2, 6.624, 5.9, 0.93, 0.0025))

nlr(fPK24b, dPK24, pNames=c("Vc", "Vt1", "Vt2", "CLo", "CLd1", "CLd2", "a"), IE=c(1, 2, 4, 7, 6, 1, 0.003), Error="P")




