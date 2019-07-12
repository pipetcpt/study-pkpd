require(wnl)
dPK02 = read.csv("data/PK02.csv", skip=1)
colnames(dPK02) = c("TIME", "DV") ; dPK02

## NCA
library(NonCompart)
sNCA(dPK02[,"TIME"], dPK02[,"DV"], dose=100, doseUnit="ug", timeUnit="min")
###

DOSE = 100

#
fPK02a = function(THETA) # Prediction function
{
  Ka = THETA[1]
  V  = THETA[2]
  K  = THETA[3]
  Cp = DOSE/V*Ka/(Ka - K)*(exp(-K*TIME) - exp(-Ka*TIME)) # eq 2:1
  return(Cp)
}
TIME = dPK02[,"TIME"]
r1 = nlr(fPK02a, dPK02, pNames=c("ka", "V", "k"), IE=c(0.1, 30, 0.05)) ; r1

#
fPK02b = function(THETA) # Prediction function
{
  Ka   = THETA[1]
  V    = THETA[2]
  K    = THETA[3]
  tlag = THETA[4]

  Cp  = DOSE/V*Ka/(Ka - K)*(exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag))) # eq 2:2
  return(Cp)
}
TIME = dPK02[,"TIME"]
r2 = nlr(fPK02b, dPK02, pNames=c("ka", "V", "k", "tlag"), IE=c(0.1, 30, 0.05, 20)) ; r2 # Error in WinNonlin and textbook

# Figure 2.3, p 480
plot(dPK02[,"TIME"], dPK02[,"DV"], xlim=c(0, 400), ylim=c(0, 2.5), xlab="Time (min)", ylab="Concentration (ug/L)", pch=16)
TIME = 0:400
lines(TIME, fPK02a(r1$Est["PE", 1:3]), lty=2)
lines(TIME, fPK02b(r2$Est["PE", 1:4]))
###


