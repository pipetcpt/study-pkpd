setwd("D:/Rt/Gab/")
require(wnl)
dPK20 = read.csv("PK20.csv", skip=1)
colnames(dPK20) = c("TIME", "DV", "ID") ; dPK20

# iv bolus 1 comp, MM elimination

require(deSolve)
PKde = function(t, y, p)
{
  if (cID == 1) {
    Km = p["Km1"]
  } else if (cID == 2) {
    Km = p["Km2"]
  }
  dy1dt = -p["Vmax"]*y[1]/(Km + y[1])/p["V"] # eq 20:1
  return(list(dy1dt))
}

# Figure 20.1, p 571
plot(0, 0.1, type="n", xlim=c(0, 8), ylim=c(10, 10000), xlab="Time (h)", ylab="Concentration (ug/L)", log="y")

IDs = unique(dPK20[,"ID"])
nID = length(IDs)
DOSE = c(25000, 100000) # ug
y = vector()
for (i in 1:nID) {
  cID <<- IDs[i]  # referencing wihtin PKde
  DATAi = dPK20[dPK20$ID == cID,]
  cAMT = DOSE[i]
  Times = DATAi[,"TIME"]
  ty = lsoda(y=c(cAMT/48.6253), times=Times, func=PKde, parms=c(V=48.6253, Vmax=37609.5, Km1=279.638, Km2=804.325))
  points(DATAi[,"TIME"], DATAi[,"DV"])
  lines(ty[,"time"], ty[,"1"])
  y = rbind(y, ty)
} ; y


fPK20 = function(THETA)
{
  V    = THETA[1]
  Vmax = THETA[2]
  Km1  = THETA[3]
  Km2  = THETA[4]

  y = vector()
  for (i in 1:nID) {
    cID <<- IDs[i] # referencing wihtin PKde
    DATAi = dPK20[dPK20$ID == cID,]
    cAMT = DOSE[i]
    Times = DATAi[,"TIME"]
    y = rbind(y, lsoda(y=c(cAMT/V), times=Times, func=PKde, parms=c(V=V, Vmax=Vmax, Km1=Km1, Km2=Km2)))
  }

  return(y[,"1"])
}
fPK20(c(48.6253, 37609.5, 279.638, 804.325))

nlr(fPK20, dPK20, pNames=c("V", "Vmax", "Km1", "Km2"), IE=c(45, 45000, 350, 1000), Error="P")

## Km corrected with fu
dPK20b = read.csv("PK20b.csv", skip=1)
colnames(dPK20b) = c("TIME", "DV", "ID") ; dPK20b


PKde2 = function(t, y, p)
{
  cfu = fu[cID]
  dy1dt = -p["Vmax"]*cfu*y[1]/(p["Km"] + cfu*y[1])/p["V"] # eq 20:1
  return(list(dy1dt))
}

## Figure 20.3, p 573
plot(0, 0.1, type="n", xlim=c(0, 20), ylim=c(1, 1000000), xlab="Time (h)", ylab="Concentration (ug/L)", log="y")

IDs = unique(dPK20b[,"ID"])
nID = length(IDs)
DOSE = c(25000, 100000, 600000) # ug
fu = c(0.03, 0.01, 0.02)

y = vector()
for (i in 1:nID) {
  cID <<- IDs[i]  # referencing wihtin PKde
  DATAi = dPK20b[dPK20b$ID == cID,]
  cAMT = DOSE[i]
  Times = DATAi[,"TIME"]
  ty = lsoda(y=c(cAMT/48.6), times=Times, func=PKde2, parms=c(V=48.6, Vmax=39108.63, Km=8.8912))
  points(DATAi[,"TIME"], DATAi[,"DV"])
  lines(ty[,"time"], ty[,"1"])
  y = rbind(y, ty)
} ; y

fPK20b = function(THETA)
{
  V    = THETA[1]
  Vmax = THETA[2]
  Km   = THETA[3]

  y = vector()
  for (i in 1:nID) {
    cID <<- IDs[i] # referencing wihtin PKde
    DATAi = dPK20b[dPK20b$ID == cID,]
    cAMT = DOSE[i]
    Times = DATAi[,"TIME"]
    y = rbind(y, lsoda(y=c(cAMT/V), times=Times, func=PKde2, parms=c(V=V, Vmax=Vmax, Km=Km)))
  }

  return(y[,"1"])
}
fPK20b(c(48.6, 39108.63, 8.8912))

nlr(fPK20b, dPK20b, pNames=c("V", "Vmax", "Km"), IE=c(45, 45000, 10), Error="P")

