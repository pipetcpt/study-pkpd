dPK02 = read.csv("data/hands-on1.csv", skip=1)
colnames(dPK02) = c("TIME", "DV") ; dPK02

## NCA
library(NonCompart)
sNCA(dPK02[,"TIME"], dPK02[,"DV"], dose=100, doseUnit="ug", timeUnit="min")
