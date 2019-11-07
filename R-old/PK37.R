setwd("D:/Rt/Gab/")
require(wnl)
dPK37 = read.csv("PK37.csv", skip=1)
colnames(dPK37) = c("Cp", "DV") ; dPK37

fPK37 = function(THETA)
{
  Vmax1 = THETA[1]
  Km1   = THETA[2]
  Vmax2 = THETA[3]
  Km2   = THETA[4]
  Rate = Vmax1*e$DATA[,"Cp"]/(Km1 + e$DATA[,"Cp"]) +  Vmax2*e$DATA[,"Cp"]/(Km2 + e$DATA[,"Cp"])
  return(Rate)
}

nlr(fPK37, dPK37, pNames=c("Vmax1", "Km1", "Vmax2", "Km2"), IE=c(250, 3, 2, 0.5))

nlr(fPK37, dPK37, pNames=c("Vmax1", "Km1", "Vmax2", "Km2"), IE=c(250, 3, 2, 0.5), Error="POIS") # flip-flop

nlr(fPK37, dPK37, pNames=c("Vmax1", "Km1", "Vmax2", "Km2"), IE=c(250, 3, 2, 0.5), Error="P")

# Figure 37.1, p 647
  Vmax1 = e$PE[1]
  Km1   = e$PE[2]
  Vmax2 = e$PE[3]
  Km2   = e$PE[4]
  sCp   = seq(0, 2000, by=1)
  Rate = Vmax1*sCp/(Km1 + sCp) +  Vmax2*sCp/(Km2 + sCp)

oldpar = par(no.readonly=TRUE)
dev.new()
par(mfrow = c(1,2))

plot(dPK37[,"Cp"], dPK37[,"DV"], xlim=c(0, 2000), ylim=c(0, 2.5), xlab="Concentration (uM)", ylab="Rate (umol/g/min)", pch=19, col="red")
lines(sCp, Rate)

plot(dPK37[,"Cp"], dPK37[,"DV"], xlim=c(0.01, 2000), ylim=c(0, 2.5), log="x", xlab="Concentration (uM)", ylab="Rate (umol/g/min)", pch=19, col="red")
lines(sCp, Rate)
###

