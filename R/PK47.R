setwd("D:/Rt/Gab/")
require(wnl)
dPK47a = read.csv("PK47_cmpd1.csv")
dPK47b = read.csv("PK47_cmpd2.csv")
colnames(dPK47a) = c("Cu", "ID", "DV", "Protein") ; dPK47a
colnames(dPK47b) = c("Cu", "ID", "DV", "Protein") ; dPK47b

fPK47 = function(THETA)
{
  n  = THETA[1]
  Ka = THETA[2]
  Cu = e$DATA[,"Cu"]  
  Pt = e$DATA[,"Protein"]
  fu = (1 - 1/(1 +  Cu/n/Pt + 1/Ka/n/Pt))*100 # eq 47:1
  return(fu)
}

nlr(fPK47, dPK47a, pNames=c("n", "Ka"), IE=c(3, 5))
nlr(fPK47, dPK47b, pNames=c("n", "Ka"), IE=c(3, 5))

# Figure 47.1
oldpar = par(no.readonly=TRUE)
dev.new()
par(mfrow=c(2,1))

plot(0.01, 0, type="n", xlim=c(0.01, 100), ylim=c(0, 100), log="x", xlab="Unbound Concentration", ylab="Free Fractio (%)")
n  = 2.833
Ka = 6.097
IDs = unique(dPK47a[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  cID = IDs[i]
  cPt = unique(dPK47a[dPK47a$ID == cID,"Protein"])
  points(dPK47a[dPK47a$ID == cID,"Cu"], dPK47a[dPK47a$ID == cID,"DV"], pch=2*i+15)
  Cu = exp(seq(log(0.01), log(100), by=0.1))  
  fu = (1 - 1/(1 +  Cu/n/cPt + 1/Ka/n/cPt))*100 # eq 47:1
  lines(Cu, fu)
}

plot(0.01, 0, type="n", xlim=c(0.01, 100), ylim=c(0, 100), log="x", xlab="Unbound Concentration", ylab="Free Fractio (%)")
n  = 1.93723
Ka = 10.2353
IDs = unique(dPK47b[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  cID = IDs[i]
  cPt = unique(dPK47b[dPK47b$ID == cID,"Protein"])
  points(dPK47b[dPK47b$ID == cID,"Cu"], dPK47b[dPK47b$ID == cID,"DV"], pch=2*i+15)
  Cu = exp(seq(log(0.01), log(100), by=0.1))  
  fu = (1 - 1/(1 +  Cu/n/cPt + 1/Ka/n/cPt))*100 # eq 47:1
  lines(Cu, fu)
}

par(oldpar)
###