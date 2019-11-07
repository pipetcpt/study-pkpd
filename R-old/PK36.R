setwd("D:/Rt/Gab/")
require(wnl)
dPK36 = read.csv("PK36.csv", skip=1)
colnames(dPK36) = c("TIME", "DV") ; dPK36

Rin = 100 / 1000
Cycle = 60
fPK36 = function(THETA)
{
  amplt = THETA[1]
  V     = THETA[2]
  CL    = THETA[3]
  Time  = dPK36[,"TIME"]
	input = Rin + amplt * Rin * cos((2*3.14/Cycle) * (Time - 60))
	Cp    = input/CL*(1 - exp((-CL/V)*Time))
	return(Cp)
}

nlr(fPK36, dPK36, pNames=c("amplt", "V", "CL"), IE=c(0.4, 2.5, 0.1))

## Figure 36.1, p 644
plot(dPK36[,"TIME"], dPK36[,"DV"], xlim=c(0,350), ylim=c(0,3), xlab="Time (min)", ylab="Concentration (uM)", pch=19, col="red")

  amplt = e$PE[1]
  V     = e$PE[2]
  CL    = e$PE[3]
  Time  = seq(0, 350, by=0.1)
	input = Rin + amplt * Rin * cos((2*3.14/Cycle) * (Time - 60))
	Cp    = input/CL*(1 - exp((-CL/V)*Time))

lines(Time, Cp)
abline(h=Rin/CL, lty=2)
####

