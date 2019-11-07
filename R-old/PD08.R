# PD 8
require(wnl)
setwd("D:/Rt/PD")

dPD8 = read.csv("PD8c.csv", header=FALSE, skip=2) # Use data in Phoenix WinNonlin
colnames(dPD8) = c("CONC", "DV")
dim(dPD8)
dPD8

# Imax model
fPD8a = function(THETA)
{
  IC50 = THETA[1]
  n    = THETA[2]
  E0   = THETA[3]
  Imax = THETA[4]
  Conc = e$DATA[,"CONC"]
  Gluc = E0 - Imax*Conc^n/(IC50^n + Conc^n)
  return(Gluc)
}

r1 = nlr(fPD8a, dPD8, c("IC50", "n", "E0", "Imax"), c(2, 1, 6, 4))
r1

r1b = wnl5(fPD8a, dPD8, c("IC50", "n", "E0", "Imax"), c(2, 1, 6, 4)) # See Table 8.2
r1b

# Log-linear model
fPD8b = function(THETA)
{
  A    = THETA[1]
  E0   = THETA[2]
  m    = THETA[3]
  Conc = e$DATA[,"CONC"]
  Gluc = E0 - m*log(A*Conc + 1)
  return(Gluc)
}

r2 = nlr(fPD8b, dPD8, c("A", "E0", "m"), c(10, 6, 1))
r2

r2b = wnl5(fPD8b, dPD8, c("A", "E0", "m"), c(10, 6, 1)) # See Table 8.2
r2b

# Figure 8.5 upper
plot(dPD8[-2:-1,"CONC"], dPD8[-2:-1,"DV"], log="x", xlim=c(0.1, 100), ylim=c(1, 7), xlab="Concentration (uM)", ylab="Response (mM)", pch=15)
abline(h=2:7, lty=3)
lines(dPD8[-2:-1,"CONC"], fPD8a(r1b$PE)[-2:-1], col="red", lty=2)
lines(dPD8[-2:-1,"CONC"], fPD8b(r2b$PE)[-2:-1])


