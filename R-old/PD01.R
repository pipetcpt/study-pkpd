# PD 1
#install.packages("wnl") # stable version
#install.packages("wnl", repos="http://pmx.amc.seoul.kr") # latest version
require(wnl)
setwd("D:/Rt/PD")

dPD1 = read.csv("PD1.csv", header=FALSE, skip=2)
colnames(dPD1) = c("ID", "Ligand", "DV", "RT", "Ka", "D")
dPD1[,"Ligand"] = 1e9*dPD1[,"Ligand"] # M -> nM
dPD1[,"DV"] = 1e9*dPD1[,"DV"] # M -> nM

Ka = 31.25 # nM See p726
D =  2.01  # nM See p726
RT = 0.396 # nM See p726

# one-site binding model
fPD1a = function(THETA)
{
  KI = THETA[1]
  Ligand = e$DATA[,"Ligand"]
  AC = Ligand/(1 + Ka*D)
  ki = KI*1e1
  DR = RT*(1 - ki*AC/(1 + ki*AC))
  return(DR)
}
e$DATA = dPD1
fPD1a(7)

r1 = nlr(fPD1a, dPD1, pNames=c("KI"), IE=c(7), Error="POIS") ; r1

# two-site binding model
fPD1b = function(THETA)
{
  F1  = THETA[1]
  K1  = THETA[2]
  K2  = THETA[3]
  Ligand = e$DATA[,"Ligand"]
  AC  = Ligand/(1 + Ka*D)
  k1  = K1*1e-3
  k2  = K2*1e-2
  Inh = F1*k1*AC/(1 + k1*AC) + (1 - F1)*k2*AC/(1 + k2*AC)
  DR  = RT*(1 - Inh)
  return(DR)
}
fPD1b(c(0.8, 9.5, 1.2))

r2 = nlr(fPD1b, dPD1, pNames=c("F1", "K1", "K2"), IE=c(0.8, 50, 20), Error="POIS") ; r2

cmpChi(r1, r2)

# Figure 1.1
plot(dPD1$Ligand, dPD1$DV, xlim=c(10,1e5), ylim=c(1e-3, 1), log="xy", pch=16, 
     xlab="Ligand Concentration (nM)", ylab="Complex Concentration (M)")
lines(dPD1$Ligand, fPD1a(r1$Est["PE","KI"]), col="red")
lines(dPD1$Ligand, fPD1b(r2$Est["PE",c("F1", "K1", "K2")]))


library(ggplot2)
y1 = fPD1a(r1$Est["PE","KI"])
y2 = fPD1b(r2$Est["PE",c("F1", "K1", "K2")])
ggplot(dPD1, aes(x=Ligand, y=DV)) +
scale_x_log10() +
scale_y_log10() +
coord_cartesian(xlim=c(10, 1e5), ylim=c(1e-3,1)) +
xlab("Ligand Concentration (nM)") + 
ylab("Complex Concentration (M)") +
geom_point() +
geom_line(aes(y=y1, color="red"), show.legend=FALSE) +
geom_line(aes(y=y2))










