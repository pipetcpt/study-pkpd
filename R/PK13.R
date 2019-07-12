setwd("D:/Rt/Gab/")
require(wnl)
dPK13 = read.csv("PK13.csv", skip=1)
colnames(dPK13) = c("TIME", "DV") ; dPK13

Dbolus = 400
Dinf = 800
TAU = 26
Ro = Dinf/TAU  # 30.76923 # 800 mg for 26 min

require(deSolve)
PKde = function(t, y, p)
{
  if (t < 26) {
    RateIn = Ro # 30.76923 = 800/26
  } else {
    RateIn = 0
  }

  dy1dt = (RateIn - p["Cl"]*y[1] - p["Cld"]*y[1] + p["Cld"]*y[2])/p["Vc"] # central, plasma
  dy2dt = (p["Cld"]*y[1] - p["Cld"]*y[2])/p["Vt"] # peripheral, tissue
  return(list(c(dy1dt, dy2dt)))
}

Times = c(0, dPK13[,"TIME"]) # simultaneous dosing <- not for this case
iTime = 2:length(Times)
y = lsoda(y=c(Dbolus/2.895, 0), times=Times, func=PKde, parms=c(Vc=2.895, Vt=2.18368, Cl=0.3447, Cld=0.1784)) ; y
plot(dPK13[,"TIME"], dPK13[,"DV"])
lines(y[,1], y[,2])

fPK13 = function(THETA)
{
  Vc   = THETA[1]
  Vt   = THETA[2]
  Cl   = THETA[3]
  Cld  = THETA[4]

  Fs = lsoda(y=c(Dbolus/Vc, 0), times=Times, func=PKde, parms=c(Vc=Vc, Vt=Vt, Cl=Cl, Cld=Cld))
  return(Fs[iTime,"1"])
}

## Additive error model
nlr(fPK13, dPK13, pNames=c("Vc", "Vt", "Cl", "Cld"),
    IE=c(3, 2, 1, 1),
    LB=c(2, 1, 0.1, 0.1),
    UB=c(4, 3, 10, 10))

wnl5(fPK13, dPK13, pNames=c("Vc", "Vt", "Cl", "Cld"),
    IE=c(3, 2, 1, 1),
    LB=c(2, 1, 0.1, 0.1),
    UB=c(4, 3, 10, 10))

## Proportional (multiplicative) error model
nlr(fPK13, dPK13, pNames=c("Vc", "Vt", "Cl", "Cld"),
    IE=c(3, 3, 0.5, 0.2),
    LB=c(1, 1, 0.1, 0.1),
    UB=c(5, 5, 1, 1), Error="P")

wnl5(fPK13, dPK13, pNames=c("Vc", "Vt", "Cl", "Cld"),
    IE=c(3, 3, 0.5, 0.2),
    LB=c(1, 1, 0.1, 0.1),
    UB=c(5, 5, 1, 1), Error="P") # different result

## Combined (mixed) error model
nlr(fPK13, dPK13, pNames=c("Vc", "Vt", "Cl", "Cld"), IE=c(2, 1, 1, 5), Error="C")



### analytical solution
fPK13b = function(THETA)
{
  Vc  = THETA[1]
  k21 = THETA[2]
  a   = THETA[3] # alpha
  b   = THETA[4] # beta
  TIME = e$DATA[,"TIME"]

  A = (k21 - a)/(b - a)
  B = (k21 - b)/(a - b)

  Cbolus = Dbolus/Vc*(A*exp(-a*TIME) + B*exp(-b*TIME))

  T1 = TIME[TIME < TAU]
  T2 = TIME[TIME >= TAU]
  Cinf1 = Ro/Vc*(A/a*(1 - exp(-a*T1)) + B/b*(1 - exp(-b*T1)))
  CT1 = Ro/Vc*A/a*(1 - exp(-a*TAU))
  CT2 = Ro/Vc*B/b*(1 - exp(-b*TAU))
  Cinf2 = CT1*exp(-a*(T2 - TAU)) + CT2*exp(-b*(T2 - TAU))

  return(Cbolus + c(Cinf1, Cinf2))
}

#fPK13b(c(2, 0.1, 0.1, 0.01))

nlr(fPK13b, dPK13, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(2, 0.1, 0.2, 0.05))
wnl5(fPK13b, dPK13, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(2, 0.1, 0.2, 0.05))

nlr(fPK13b, dPK13, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(2, 0.1, 0.2, 0.05), Error="P")
wnl5(fPK13b, dPK13, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(2, 0.1, 0.2, 0.05), Error="P")

nlr(fPK13b, dPK13, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(2, 0.1, 0.2, 0.05), Error="C")

nlr(fPK13b, dPK13, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(2, 0.1, 0.1, 0.01), Error="POIS")
wnl5(fPK13b, dPK13, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(2, 0.1, 0.1, 0.01), Error="POIS")




