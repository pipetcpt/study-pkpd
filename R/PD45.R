# PD 45
require(wnl)
setwd("D:/Rt/PD")

dPD45 = read.csv("PD45.csv", header=FALSE, skip=2)
colnames(dPD45) = c("Cp", "DV")

Cp = dPD45[,"Cp"]
fPD45 = function(THETA)
{
  Kd   = THETA[1]
  n    = THETA[2]
  Bmax = THETA[3]
  B = Bmax*Cp^n/(Kd^n + Cp^n)
  return(B) 
}

r1 = nlr(fPD45, dPD45, c("Kd", "n", "Bmax"), c(3, 1, 100),
         SecNames=c("T1/2offA", "T1/2offB"), 
         SecForms=c(~log(2)/Kd/0.7, ~log(2)/Kd/0.027)) ; r1

# Figure 45.1
plot(dPD45[,"Cp"], dPD45[,"DV"], log="x", xlim=c(0.1, 200), ylim=c(0, 100), xlab="Concentration", ylab="Binding (%)", pch=16)
abline(h=1:10*10, lty=3)
abline(v=as.vector(outer(1:9, 10^(-1:3))), lty=3)
lines(dPD45[,"Cp"], r1$Prediction)
