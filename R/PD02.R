# PD 2
require(wnl)
setwd("D:/Rt/PD")
dPD2 = read.csv("PD2.csv")
dPD2[,"DV"] = dPD2[,"Binding"] ; dPD2

IDs = unique(dPD2[,"Dataset"]) ; IDs
nID = length(IDs) ; nID

# One-site binding model
fPD2a = function(THETA)
{
  IC51 = THETA[5]

  Res = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    AC = 10^dPD2[dPD2$Dataset == cID, "Ligand"]
    Inh = AC/(IC51 + AC)
    Res = c(Res, THETA[i]*(1 - Inh))
  }
  return(Res)
}
y1 = fPD2a(c(7039.76, 5685.25, 3097.89, 1392.52, 0.549916)) ; y1

# Figure 2.2
dev.new()
plot(dPD2$Ligand, dPD2$Binding, xlab="Log(Ligand concentration)", ylab="Binding", pch=16)
for (i in 1:nID) {
  lines(dPD2[dPD2$Dataset == IDs[i], "Ligand"], Y[dPD2$Dataset == IDs[i]])
}

require(ggplot2)
dev.new()
ggplot(dPD2, aes(x=Ligand, y=Binding, group=Dataset)) +
geom_point() +
geom_line(aes(y=y1, group=Dataset))
###

r1 = nlr(fPD2a, dPD2, c("B1", "B2", "B3", "B4", "IC51"), c(7000, 5500, 3000, 1500, 1))
r1
w1 = wnl5(fPD2a, dPD2, c("B1", "B2", "B3", "B4", "IC51"), c(7000, 5500, 3000, 1500, 1))
w1
# Two-site binding model
fPD2b = function(THETA)
{
  F1   = THETA[5]
  IC51 = THETA[6]
  IC52 = THETA[7]

  Res = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    AC = 10^dPD2[dPD2$Dataset == cID, "Ligand"]
    Inh = F1*AC/(IC51 + AC) + (1 - F1)*AC/(IC52 + AC)
    Res = c(Res, THETA[i]*(1 - Inh))
  }
  return(Res)
}
y2 = fPD2b(c(8050.94, 6512.13, 3542.05, 1594.11, 0.613739, 0.0250628, 37.0361)) ; y2

# Figure 2.1
dev.new()
plot(dPD2$Ligand, dPD2$Binding, xlab="Log(Ligand concentration)", ylab="Binding", pch=16)
for (i in 1:nID) {
  lines(dPD2[dPD2$Dataset == IDs[i],"Ligand"], Y[dPD2$Dataset == IDs[i]])
}
abline(v=log10(0.250628), lty=2)
abline(v=log10(37.0361), lty=2)

require(ggplot2)
dev.new()
ggplot(dPD2, aes(x=Ligand, y=Binding, group=Dataset)) +
geom_point() +
geom_line(aes(y=y2, group=Dataset))
###

r2 = nlr(fPD2b, dPD2, c("B1", "B2", "B3", "B4", "F1", "IC51", "IC52"), c(7000, 5500, 3000, 1500, 0.5, 0.1, 10))
r2

w2 = wnl5(fPD2b, dPD2, c("B1", "B2", "B3", "B4", "F1", "IC51", "IC52"), c(7000, 5500, 3000, 1500, 0.5, 0.1, 10))
w2

cmpChi(r1, r2)

