setwd("D:/Rt/Gab/")
require(wnl)
dPK53 = read.csv("PK53_data.csv", skip=1)
colnames(dPK53) = c("TIME", "DV") ; dPK53

#dPK53dose = read.csv("PK53_dose.csv", skip=1)
#colnames(dPK53dose) = c("AMT", "Start", "End", "Rate") ; dPK53dose
#R1 = cbind(TIME=dPK53dose[,"Start"], RATE=dPK53dose[,"Rate"])
#R2 = cbind(TIME=dPK53dose[,"End"], RATE=0)
#infRate = rbind(R1, R2)
#infRate = infRate[order(infRate[,"TIME"]),] ; infRate

require(deSolve)
PKde = function(t, y, p)
{
  RateIn = 0
  if (t <= 0.41667)            RateIn = 1.848
  if (t > 72.17  & t < 72.67)  RateIn = 15.4
  if (t > 144.17 & t < 144.67) RateIn = 154
  if (t > 216.6  & t < 217.0)  RateIn = 642.5
  if (t > 288.52 & t < 289.02) RateIn = 1542.0

  dy1dt = (RateIn - p["CL"]*y[1] - p["CLd"]*y[1] + p["CLd"]*y[2])/p["Vc"]
  dy2dt = (p["CLd"]*y[1] - p["CLd"]*y[2])/p["Vt"]
  return(list(c(dy1dt, dy2dt)))
}

TIME = sort(unique(c(0:2000, dPK53[,"TIME"])))
iTime = TIME %in% dPK53[,"TIME"]
y = rk4(y=c(0, 0), times=TIME, func=PKde, parms=c(Vc=2.14, Vt=1.5858, CL=0.00541264, CLd=0.0164))
y[iTime,]
nrow(y)
nrow(y[iTime,])
length(dPK53[,"DV"])

# Fiture 53.1, p 722
plot(dPK53[,"TIME"], dPK53[,"DV"], xlim=c(0,2000), ylim=c(0.1, 1000), log="y", xlab="Time (h)", ylab="Concentration (uM)", pch=19)
lines(y[-1,"time"],y[-1,"1"])
###

#cTIME = c(0, dPK53[,"TIME"])
#iTime = 2:length(cTIME)

fPK53 = function(THETA)
{
  Vc  = THETA[1]
  Vt  = THETA[2]
  CL  = THETA[3]
  CLd = THETA[4]
  y = rk4(y=c(0, 0), times=TIME, func=PKde, parms=c(Vc=Vc, Vt=Vt, CL=CL, CLd=CLd))
  return(y[iTime, "1"])
}
y2 = fPK53(c(2.02889, 1.75136, 0.00542873, 0.023733))
length(y2)
length(dPK53[,"DV"])
y2
dPK53[,"DV"]

nlr(fPK53, dPK53, pNames=c("Vc", "Vt", "CL", "CLd"), IE=c(3, 2, 0.01, 0.1)) # fitting failure, see NONMEM, different, note large sigma
nlr(fPK53, dPK53, pNames=c("Vc", "Vt", "CL", "CLd"), IE=c(3, 2, 0.01, 0.1), Error="POIS") # fitting failure, see NONMEM, same
nlr(fPK53, dPK53, pNames=c("Vc", "Vt", "CL", "CLd"), IE=c(3, 2, 0.01, 0.1), Error="P") # fitting failure, see NONMEM, different
