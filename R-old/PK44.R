setwd("D:/Rt/Gab/")
require(wnl)
dPK44 = read.csv("PK44.csv", skip=1)
colnames(dPK44) = c("Cp", "DV", "ID", "Exposure") ; dPK44

# Competitive model, better fitting
fPK44 = function(THETA)
{
  Ki = THETA[1]
  Km = THETA[2]
  Vm = THETA[3]
  Cp = dPK44[,"Cp"]
  Ex = dPK44[,"Exposure"]
  Rate = Vm*Cp/(Km*(1 + Ex/Ki) + Cp)
  return(Rate)
}
fPK44(c(5.506, 11.3192, 99.8))

# Figure 44.1, p 675
oldpar = par(no.readonly=TRUE)
dev.new()
par(mfrow=c(1,2))
IDs = unique(dPK44[,"ID"])
nID = length(IDs)

Bullets = c(1, 0, 18, 17, 6, 16)

# Left figure
plot(0, 0, type="n", xlim=c(0, 1000), ylim=c(0, 100), xlab="Concentration (uM)", ylab="Rate (mole/min/mg protein)")
for (i in 1:nID) {
  cID = IDs[i]
  tx = dPK44[dPK44$ID == cID, "Cp"]
  ty = dPK44[dPK44$ID == cID, "DV"]
  points(tx, ty, pch=Bullets[i])

  Ex = unique(dPK44[dPK44$ID == cID, "Exposure"])
  Cp = seq(0, 1000, by=1)
  y2 = 99.8*Cp/(11.3192*(1 + Ex/5.506) + Cp)
  lines(Cp, y2)
}

# Right figure
plot(0.1, 0, type="n", xlim=c(1, 1000), ylim=c(0, 100), log="x", xlab="Concentration (uM)", ylab="Rate (mole/min/mg protein)")
for (i in 1:nID) {
  cID = IDs[i]
  tx = dPK44[dPK44$ID == cID, "Cp"]
  ty = dPK44[dPK44$ID == cID, "DV"]
  points(tx, ty, pch=Bullets[i])

  Ex = unique(dPK44[dPK44$ID == cID, "Exposure"])
  Cp = exp(seq(log(1), log(1000), by=0.1))
  y2 = 99.8*Cp/(11.3192*(1 + Ex/5.506) + Cp)
  lines(Cp, y2)
}
par(oldpar)
### end of plotting

nlr(fPK44, dPK44, pNames=c("Ki", "Km", "Vm"), IE=c(10, 10, 100))


# Noncompetitive model, worse fitting
fPK44b = function(THETA)
{
  Ki = THETA[1]
  Km = THETA[2]
  Vm = THETA[3]
  Cp = dPK44[,"Cp"]
  Ex = dPK44[,"Exposure"]
  Rate = Vm*Cp/(Km + Cp)*Ki/(Ki + Ex)
  return(Rate)
}

nlr(fPK44b, dPK44, pNames=c("Ki", "Km", "Vm"), IE=c(10, 10, 100))

# Figure 44.2, p 676
oldpar = par(no.readonly=TRUE)
dev.new()
par(mfrow=c(1,2))
IDs = unique(dPK44[,"ID"])
nID = length(IDs)

Bullets = c(1, 0, 18, 17, 6, 16)

# Left figure
plot(0, 0, type="n", xlim=c(0, 1000), ylim=c(0, 100), xlab="Concentration (uM)", ylab="Rate (mole/min/mg protein)")
for (i in 1:nID) {
  cID = IDs[i]
  tx = dPK44[dPK44$ID == cID, "Cp"]
  ty = dPK44[dPK44$ID == cID, "DV"]
  points(tx, ty, pch=Bullets[i])

  Ex = unique(dPK44[dPK44$ID == cID, "Exposure"])
  Cp = seq(0, 1000, by=1)
  y2 = 99.8*Cp/(42.7 + Cp)*98.8/(98.8 + Ex) # Different part with Figure 44.1
  lines(Cp, y2)
}

# Right figure
plot(0.1, 0, type="n", xlim=c(1, 1000), ylim=c(0, 100), log="x", xlab="Concentration (uM)", ylab="Rate (mole/min/mg protein)")
for (i in 1:nID) {
  cID = IDs[i]
  tx = dPK44[dPK44$ID == cID, "Cp"]
  ty = dPK44[dPK44$ID == cID, "DV"]
  points(tx, ty, pch=Bullets[i])

  Ex = unique(dPK44[dPK44$ID == cID, "Exposure"])
  Cp = exp(seq(log(1), log(1000), by=0.1))
  y2 = 99.8*Cp/(42.7 + Cp)*98.8/(98.8 + Ex) # Different part with Figure 44.1
  lines(Cp, y2)
}
par(oldpar)
### end of plotting
