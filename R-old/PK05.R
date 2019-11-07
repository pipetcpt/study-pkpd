require(wnl)

dPK05 = read.csv("data/PK05.csv", skip=1)
colnames(dPK05) = c("TIME", "DV", "CMT") ; dPK05

AMT = 250 # mg

## Closed form model
fPK05a = function(THETA)
{
  V  = THETA[1]
  Cl = THETA[2]
  fe = THETA[3]
  
  Div = AMT
  
  T1 = e$DATA[e$DATA[,"CMT"]==1,"TIME"]
  Civ = Div / V * exp(-Cl/V*T1)       # Eq 5:1
  
  T2 = e$DATA[e$DATA[,"CMT"]==2,"TIME"]
  Xu = fe * Div * (1 - exp(-Cl/V*T2)) # Eq 5:2
  
  return(c(Civ, Xu))  
}

nlr(fPK05a, dPK05, pNames=c("V", "Cl", "fe"), IE=c(10.5, 1.15, 0.3), Error="P")
wnl5(fPK05a, dPK05, pNames=c("V", "Cl", "fe"), IE=c(10.5, 1.15, 0.3), Error="PROP")

## Differential Equation Model
require(deSolve)

ivPK1c = function(t, y, p)
{
  dy1dt = -p["Cl"]/p["V"]*y[1]     # Eq 5:3
  dy2dt = p["fe"]*p["Cl"]*y[1]     # Eq 5:4
  return(list(c(dy1dt, dy2dt)))
}

Times = c(0, dPK05[dPK05[,"CMT"]==1,"TIME"])
lsoda(y=c(AMT/10.72, 0), times=Times, func=ivPK1c, parms=c(V=10.72, Cl=1.23, fe=0.3539))

iTime = 2:length(Times)

fPK05b = function(THETA)
{
  Fs = lsoda(y=c(AMT/THETA[1],0), times=Times, func=ivPK1c, parms=c(V=THETA[1], Cl=THETA[2], fe=THETA[3]))
  return(c(Fs[iTime,2], Fs[iTime,3]))
}

nlr(fPK05b, dPK05, pNames=c("V", "Cl", "fe"), IE=c(10.5, 1.15, 0.3), Error="P")
wnl5(fPK05b, dPK05, pNames=c("V", "Cl", "fe"), IE=c(10.5, 1.15, 0.3), Error="PROP")
