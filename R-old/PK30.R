setwd("D:/Rt/Gab/")
require(wnl)
dPK30 = read.csv("PK30.csv")
colnames(dPK30) = c("TIME", "DV") ; dPK30

DOSE = 40
TIME = dPK30[,"TIME"]

require(deSolve)
PKde = function(t, y, p)
{
  dy1dt = -p["Ka"]*y[1]
  dy2dt = (p["Syn"] + p["Ka"]*y[1] - p["CL"]*y[2])/p["V"]
  return(list(c(dy1dt, dy2dt)))
}
lsoda(y=c(DOSE, 0.78/0.028), times=TIME, func=PKde, parms=c(Ka=0.54, V=0.1, CL=0.028, Syn=0.78))

## No need of the above differential equation.


fPK30 = function(THETA)
{
  Ka =  THETA[1]
  V  =  THETA[2]
  CL =  THETA[3]
  Syn = THETA[4]
  Cp = Syn/CL + 40*Ka/V/(Ka - CL/V)*(exp(-CL/V*TIME) - exp(-Ka*TIME))
  return(Cp)  
}
fPK30(c(0.54, 0.1, 0.028, 0.78))

nlr(fPK30, dPK30, pNames=c("Ka", "V", "CL", "Syn"), IE=c(0.5, 0.2, 0.1, 1))

# compare with the following, flip-flop
nlr(fPK30, dPK30, pNames=c("Ka", "V", "CL", "Syn"), IE=c(1, 0.5, 0.1, 1))
e$PE[3]/e$PE[2]
