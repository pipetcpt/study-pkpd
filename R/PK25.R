setwd("D:/Rt/Gab/")
require(wnl)
dPK25 = read.csv("PK25.csv", skip=1, as.is=TRUE)
colnames(dPK25) = c("TIME", "DV", "CMT") ; dPK25
dPK25 = dPK25[dPK25$DV != "Missing",]
dPK25[,"DV"] = as.numeric(dPK25[,"DV"]) ; dPK25

T1 = dPK25[dPK25$CMT==1, "TIME"]
T2 = dPK25[dPK25$CMT==2, "TIME"]

fPK25 = function(THETA)
{
  A     = THETA[1]
  alpha = THETA[2]
  B     = THETA[3]
  beta  = THETA[4]
  CLr   = THETA[5] 
  
  Cp = A*exp(-alpha*T1) + B*exp(-beta*T1)
  Xu = CLr * (A/alpha*(1 - exp(-alpha*T2)) + B/beta*(1 - exp(-beta*T2)))
  return(c(Cp, Xu))
}
fPK25(c(391, 0.9, 13, 0.1, 5.76))

# Figure 25.1 Left, p 597
plot(0, 0.1, type="n", xlim=c(0,25), ylim=c(1, 10000), log="y", xlab="Time (h)", ylab="Conc (ug/L), Amount (ug)")
i1 = dPK25$CMT==1
i2 = dPK25$CMT==2
points(dPK25[i1,"TIME"], dPK25[i1,"DV"], pch=19)
points(dPK25[i2,"TIME"], dPK25[i2,"DV"], pch=15)

T1 = seq(0, 25, by=0.1)
T2 = T1[2:length(T1)]
y = fPK25(c(391, 0.9, 13, 0.1, 5.76))
lines(T1, y[1:length(T1)])
lines(T2, y[(length(T1) + 1):(length(T1)+length(T2))])
###
 
T1 = dPK25[dPK25$CMT==1, "TIME"]
T2 = dPK25[dPK25$CMT==2, "TIME"]

nlr(fPK25, dPK25, pNames=c("A", "alpha", "B", "beta", "CLr"), IE=c(500, 1, 20, 0.1, 1), Error="POIS")

