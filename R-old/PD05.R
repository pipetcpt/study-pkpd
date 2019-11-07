# PD 5
require(wnl)
setwd("D:/Rt/PD")

dPD5 = read.csv("PD5.csv", header=FALSE, skip=2)
colnames(dPD5) = c("TIME", "ID", "DV", "AMT")
dPD5

require(deSolve)

AMTs = c(4000, 16000, 80000)
Tinf = 6 # hr
Rinf = AMTs/Tinf

V  = 40 # L
Ke = 0.9 # /h
 
fPD5de = function(t, y, p)
{
  if (t < Tinf) {
    RateIn = Rinf[cID] # cID: current ID, external
  } else {
    RateIn = 0
  }
  
  dy1dt = RateIn/V - Ke*y[1]
  Inh = 1 - p["Imax"]*y[1]/(p["IC50"] + y[1])
  dy2dt = p["Kin"] -  Inh*p["Kout"]*y[2]
  return(list(c(dy1dt, dy2dt)))
}

# Figure 5.1
plot(0, 0, type="n", xlim=c(0, 12), ylim=c(30, 90), xlab="Time (h)", ylab="Response")
abline(h=10*(4:9), lty=3)
abline(v=6, col="red", lty=3)
IDs = unique(dPD5[,"ID"])
nID = length(IDs)
y = matrix(ncol=3)
for (i in 1:nID) {
  cID <<- IDs[i] # referencing within fPD5de
  cTIME = dPD5[dPD5$ID == cID, "TIME"]
  cy = lsoda(y=c(0, 19.0971/0.4318), times=cTIME, func=fPD5de, parms=c(Kin=19.0971, Kout=0.4318, Imax=0.6526, IC50=94.1096))
  points(dPD5[dPD5$ID == cID, "TIME"], dPD5[dPD5$ID == cID, "DV"], pch=16)
  lines(cTIME, cy[, "2"])
  dPD5[dPD5$ID==cID,"CONC"] = cy[,"1"]
  dPD5[dPD5$ID==cID,"PRED"] = cy[,"2"]
} ; dPD5
###


dPD5 = dPD5[!is.na(dPD5$DV),] # removing NA points for fitting

# Figure 5.2
dev.new()
plot(0.1, 0, type="n", xlim=c(0.1, 1000), ylim=c(40, 80), log="x", xlab="Concentration", ylab="Response")
abline(h=seq(45, 80, 5), lty=3)
for (i in 1:nID) {
  cID = IDs[i]
  lines(dPD5[dPD5$ID == cID, "CONC"], dPD5[dPD5$ID == cID, "DV"], type="o", pch=15+i)
}
###

fPD5 = function(THETA)
{
  Kin  = THETA[1]
  Kout = THETA[2]
  Imax = THETA[3]
  IC50 = THETA[4]
  y = vector()
  for (i in 1:nID) {
    cID <<- IDs[i] # referencing within PKde
    cTIME = dPD5[dPD5$ID==cID,"TIME"]
    cy = lsoda(y=c(0, Kin/Kout), times=cTIME, func=fPD5de, parms=c(Kin=Kin, Kout=Kout, Imax=Imax, IC50=IC50))
    y = c(y, cy[,"2"])
  }
  return(y)
}
fPD5(c(19.0971, 0.4318, 0.6526, 94.1096))

r1 = nlr(fPD5, dPD5, pNames=c("Kin", "Kout", "Imax", "IC50"), IE=c(7, 0.2, 0.25, 30))
r1


