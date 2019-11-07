# PD 11
require(wnl)
setwd("D:/Rt/PD")

dPD11 = read.csv("PD11.csv")
colnames(dPD11) = c("Distance", "DV")

X = dPD11[,"Distance"]

# Gompertz
fPD11a = function(THETA)
{
  Alpha = THETA[1]
  Beta  = THETA[2]
  Gamma = THETA[3]
  Y = Alpha*exp(-exp(Beta - Gamma*X))
  return(Y)
}

# Logistic
fPD11b = function(THETA)
{
  Alpha = THETA[1]
  Beta  = THETA[2]
  Gamma = THETA[3]
  Delta = THETA[4]
  Y = Alpha/(1 + exp(Beta - Gamma*X))
  return(Y)
}

# Weibull
fPD11c = function(THETA)
{
  Alpha = THETA[1]
  Beta  = THETA[2]
  Gamma = THETA[3]
  Delta = THETA[4]
  Y = Alpha - Beta*exp(-Gamma*X^Delta)
  return(Y)
}

# Richards
fPD11d = function(THETA)
{
  Alpha = THETA[1]
  Beta  = THETA[2]
  Gamma = THETA[3]
  Delta = THETA[4]
  Y = Alpha/(1 + exp(Beta - Gamma*X))^(1/Delta)
  return(Y)
}

# MMF
fPD11e = function(THETA)
{
  Alpha = THETA[1]
  Beta  = THETA[2]
  Gamma = THETA[3]
  Delta = THETA[4]
  Y = (Beta*Gamma + Alpha*X^Delta)/(Gamma + X^Delta)
  return(Y)
}

# Hill
fPD11f = function(THETA)
{
  Alpha = THETA[1]
  Beta  = THETA[2]
  Gamma = THETA[3]
  Delta = THETA[4]
  Y = Alpha + Beta*X^Delta/(Gamma^Delta + X^Delta)
  return(Y)
}

r1 = nlr(fPD11a, dPD11, c("Alpha", "Beta" ,"Gamma"), c(50, 1, 1)) ; r1
r2 = nlr(fPD11b, dPD11, c("Alpha", "Beta" ,"Gamma"), c(50, 1, 1)) ; r2
r3 = nlr(fPD11c, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(50, 20, 0.002, 1)) ; r3
r4 = nlr(fPD11d, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(10, 3, 1, 2)) ; r4
r5 = nlr(fPD11e, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(20, 1, 1000, 5)) ; r5
r6 = nlr(fPD11f, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(50, 1, 1, 1)) ; r6

w1 = wnl5(fPD11a, dPD11, c("Alpha", "Beta" ,"Gamma"), c(50, 1, 1)) ; w1
w2 = wnl5(fPD11b, dPD11, c("Alpha", "Beta" ,"Gamma"), c(50, 1, 1)) ; w2 # Errata of Gamma (0.0621 -> 0.622) in the Book
w3 = wnl5(fPD11c, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(50, 20, 0.002, 1)) ; w3 # Errata Delta 3.17 -> 3.18
w4 = wnl5(fPD11d, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(10, 3, 1, 2)) ; w4 # Errata Beta 5.71 -> 5.69
w5 = wnl5(fPD11e, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(20, 1, 1000, 5)) ; w5 # Errata Gamma 5300 -> 5590
w6 = wnl5(fPD11f, dPD11, c("Alpha", "Beta" ,"Gamma", "Delta"), c(1, 50, 1, 1)) ; w6

# Figure 11.3
plot(dPD11[,"Distance"], dPD11[,"DV"], log="x", xlim=c(1, 20), ylim=c(0, 25), xlab="Distance", ylab="Water content", pch=16)
abline(h=seq(5, 25, 5), lty=3)
lines(X, fPD11c(w3$PE))
lines(X, fPD11f(w6$PE), col="red", lty=2)
legend(2.2, 19, legend=c("Hill equation", "Weibull equation"), col=c("red", "black"), lty=c(2,1), bty="n")

# Figure 11.5
Resid3 = dPD11[,"DV"] - fPD11c(w3$PE)
Resid6 = dPD11[,"DV"] - fPD11f(w6$PE)
plot(X, Resid3, type="l", xlim=c(0, 16), ylim=c(-1, 2), xlab="Distance", ylab="Residual")
lines(X, Resid6, lty=2, col="red", lwd=2)
abline(h=seq(-0.5, 2, 0.5), lty=3)
abline(h=0)
legend(11, 1.9, legend=c("Hill equation", "Weibull equation"), col=c("red", "black"), lty=c(2,1), bty="n")
###


