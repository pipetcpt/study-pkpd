setwd("D:/Rt/Gab/")
require(wnl)
dPK33 = read.csv("PK33.csv", skip=1)
colnames(dPK33) = c("TIME", "DV") ; dPK33

Dose = 15890

require(deSolve)
PKde = function(t, y, p)
{
  V     = p["V"]
  CL    = p["CL"]
  Ftau  = p["Ftau"]
  Stau  = p["Stau"]
  SDose = p["SDose"]
  
  if (t <= Ftau) {
    Frate = (Dose - SDose)/Ftau  # eq 33:2
  } else {
    Frate = 0
  }
   
  if (t <= Stau) {
    Srate = SDose/Stau
  } else {
    Srate = 0
  }
  
  dy1dt = (Frate + Srate - CL*y[1])/V # eq 33:1
  
  return(list(c(dy1dt)))
}

TIME = dPK33[,"TIME"]
y = lsoda(y=c(2), times=TIME, func=PKde, parms=c(V=253.7, CL=79.3431, Ftau=7.35265, Stau=19.6, SDose=10950.5))
y # note initial value and value at time 0 in figure 33.1 

# Figure 33.1, p 635
plot(dPK33[,"TIME"], dPK33[,"DV"], xlim=c(0, 25), ylim=c(2, 16), xlab="Time (h)", ylab="Concentration (ug/L)", pch=19, col="red")
lines(y[,"time"], y[,"1"])
###

# Frank version
y = lsoda(y=c(0), times=TIME, func=PKde, parms=c(V=253.7, CL=79.3431, Ftau=7.35265, Stau=19.6, SDose=10950.5))
dev.new()
plot(dPK33[,"TIME"], dPK33[,"DV"], xlim=c(0, 25), ylim=c(0, 16), xlab="Time (h)", ylab="Concentration (ug/L)", pch=19, col="red")
lines(y[,"time"], y[,"1"])
###


fPK33 = function(THETA)
{
  V     = THETA[1]
  CL    = THETA[2]
  Ftau  = THETA[3]
  Stau  = THETA[4]
  SDose = THETA[5]
  y = lsoda(y=c(0), times=TIME, func=PKde, parms=c(V=V, CL=CL, Ftau=Ftau, Stau=Stau, SDose=SDose))

  return(y[,"1"])
}
y = fPK33(c(253.7, 79.3431, 7.35265, 19.6,10950.5)) ; y
length(y)
length(dPK33[,"DV"])


nlr(fPK33, dPK33, pNames=c("V", "CL", "Ftau", "Stau", "SDose"), IE=c(140, 78, 6, 17, 10000)) # note unstable V

wnl5(fPK33, dPK33, pNames=c("V", "CL", "Ftau", "Stau", "SDose"), IE=c(140, 78, 6, 17, 10000), Error="POIS")

fPK33 = function(THETA)
{
  V     = THETA[1]
  CL    = THETA[2]
  Ftau  = THETA[3]
  Stau  = THETA[4]
  SDose = THETA[5]
  y = lsoda(y=c(2), times=TIME, func=PKde, parms=c(V=V, CL=CL, Ftau=Ftau, Stau=Stau, SDose=SDose))

  return(y[,"1"])
}

nlr(fPK33, dPK33, pNames=c("V", "CL", "Ftau", "Stau", "SDose"), IE=c(140, 78, 6, 17, 10000), Error="POIS")

wnl5(fPK33, dPK33, pNames=c("V", "CL", "Ftau", "Stau", "SDose"), IE=c(140, 78, 6, 17, 10000), Error="POIS")
