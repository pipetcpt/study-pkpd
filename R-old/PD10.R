# PD 10
require(wnl)
setwd("D:/Rt/PD")

dPD10 = read.csv("PD10.csv", header=FALSE, skip=2)
colnames(dPD10) = c("TIME", "DV")

V   = 0.7633
V2  = 1.72876
V3  = 3.43857
Cl  = 6.2417
Cl2 = 5.4595
Cl3 = 0.85806

require(deSolve)
fPD10ade = function(t, y, p)
{
  RateIn = 0
  if (t <= 0.3) RateIn = 21333.33
  if (t > 1.051  & t < (1.051 + 0.2667)) RateIn = 22122.23
  dy1dt = (RateIn - Cl*y[1] - Cl2*(y[1] - y[2]) - Cl3*(y[1] - y[3]))/V
  dy2dt = Cl2*(y[1] - y[2])/V2
  dy3dt = Cl3*(y[1] - y[3])/V3
  dy4dt = p["Kin"]*(1 - p["Imax"]*y[1]^p["n"]/(y[1]^p["n"] + p["IC50"]^p["n"])) - p["Kout"]*y[4] 
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt)))
}

Kin  = 143.313
Kout = 12.5527
Imax = 0.830642
IC50 = 244.203
n    = 7.06292

TIMEs = dPD10[,"TIME"]
y =  lsoda(c(0, 0, 0, Kin/Kout), TIMEs, fPD10ade, c(Kin=Kin, Kout=Kout, Imax=Imax, IC50=IC50, n=n))
y
 
# Figure 10.1
plot(y[-1,"time"], y[-1,"1"], log="y", xlim=c(0, 2.5), ylim=c(100, 3000), xlab="Time (h)", ylab="Concentration ug/L", type="l")
abline(h=c(seq(100, 900, 100), seq(1000, 3000, 1000)), lty=3)
lines(c(0, 0.3), c(100, 100), lwd=10)
lines(c(1.051, 1.051 + 0.2667), c(100, 100), lwd=10)

# Figure 10.2 Right
dev.new()
plot(y[,"1"], dPD10[,"DV"], type="o", xlim=c(0, 2500), ylim=c(0, 14), xlab="Concentration (ug/L)", ylab="Response", pch=16)
abline(h=seq(2, 14, 2), lty=3)

# Figure 10.3
dev.new()
plot(dPD10[,"TIME"], dPD10[,"DV"], xlim=c(0, 2.5), ylim=c(0, 12), xlab="Time (h)", ylab="Response", pch=16)
lines(y[,"time"], y[,"4"])
abline(h=seq(2, 12, 2), lty=3)
###

fPD10a = function(THETA)
{
  Kin  = THETA[1]
  Kout = THETA[2]
  Imax = THETA[3]
  IC50 = THETA[4]
  n    = THETA[5]
  y =  lsoda(c(0, 0, 0, Kin/Kout), TIMEs, fPD10ade, c(Kin=Kin, Kout=Kout, Imax=Imax, IC50=IC50, n=n))
  return(y[,"4"])
}
r1 = nlr(fPD10a, dPD10, c("Kin", "Kout", "Imax", "IC50", "n"), c(140, 13, 0.82, 300, 6.6))
r1
r1b = wnl5(fPD10a, dPD10, c("Kin", "Kout", "Imax", "IC50", "n"), c(140, 13, 0.82, 300, 6.6))
r1b

fPD10bde = function(t, y, p)
{
  RateIn = 0
  if (t <= 0.3) RateIn = 21333.33
  if (t > 1.051  & t < (1.051 + 0.2667)) RateIn = 22122.23
  dy1dt = (RateIn - Cl*y[1] - Cl2*(y[1] - y[2]) - Cl3*(y[1] - y[3]))/V
  dy2dt = Cl2*(y[1] - y[2])/V2
  dy3dt = Cl3*(y[1] - y[3])/V3
  dy4dt = p["Kin"] - p["Kout"]*(1 + p["Emax"]*y[1]^p["n"]/(y[1]^p["n"] + p["EC50"]^p["n"]))*y[4] 
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt)))
}

fPD10b = function(THETA)
{
  Kin  = THETA[1]
  Kout = THETA[2]
  Emax = THETA[3]
  EC50 = THETA[4]
  n    = THETA[5]
  y =  lsoda(c(0, 0, 0, Kin/Kout), TIMEs, fPD10bde, c(Kin=Kin, Kout=Kout, Emax=Emax, EC50=EC50, n=n))
  return(y[,"4"])
}
r2 = nlr(fPD10b, dPD10, c("Kin", "Kout", "Emax", "EC50", "n"), c(40, 3, 6, 350, 13))
r2
r2b = wnl5(fPD10b, dPD10, c("Kin", "Kout", "Emax", "EC50", "n"), c(40, 3, 6, 350, 13))
r2b