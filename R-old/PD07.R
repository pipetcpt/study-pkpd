# PD 7
require(wnl)
setwd("D:/Rt/PD")

dPD7 = read.csv("PD7.csv", header=FALSE, skip=2)
colnames(dPD7) = c("TIME", "DV", "ID", "AMT")
dPD7 = dPD7[!is.na(dPD7$DV),] # DV should not have NAs.
dim(dPD7)

IDs = unique(dPD7[,"ID"])
nID = length(IDs)

V  = 28.6 # L
Ke = 2.8 # /hr

AMTs = c(6400, 32000, 160000)
Tinf = 4 # hr
Rinf = AMTs / Tinf

require(deSolve)
fPD7de = function(t, y, p)
{
  if (t < Tinf) {
    RateIn = Rinf[cID] # cID: current ID, external
  } else {
    RateIn = 0
  }
  dy1dt = RateIn/V - Ke*y[1]
  Stim  = 1 + p["Emax"]*y[1]/(p["EC50"] + y[1])
  dy2dt = p["Kin"] -  Stim*p["Kout"]*y[2]
  return(list(c(dy1dt, dy2dt)))
}

# Figure 7.1
plot(0, 0, type="n", xlim=c(0, 12), ylim=c(0, 35), xlab="Time (h)", ylab="Response")
abline(h=seq(5, 35, 5), lty=3)
Kin  = 27.4742
Kout = 0.923285
Emax = 4.48501
EC50 = 49.8007
for (i in 1:nID) {
  cID <<- IDs[i] # Referencing within fPD7de
  cTIME = c(0, dPD7[dPD7$ID == cID, "TIME"]) # TIME should include 0.
  cy = lsoda(c(0, Kin/Kout), cTIME, fPD7de, c(Kin=Kin, Kout=Kout, Emax=Emax, EC50=EC50))
  points(cTIME[-1], dPD7[dPD7$ID == cID, "DV"], pch=16)
  lines(cTIME, cy[, "2"]) 
  dPD7[dPD7$ID==cID,"CONC"] = cy[-1,"1"]
  dPD7[dPD7$ID==cID,"PRED"] = cy[-1,"2"]
} ; dPD7
###

# Figure 7.2
dev.new()
plot(0, 1, type="n", log="y", xlim=c(0, 8), ylim=c(1, 1000), xlab="Time (h)", ylab="Simulated Concentration")
abline(h=c(1:9, seq(10, 90, 10), seq(100, 1000, 100)), lty=3)
for (i in 1:nID) {
  cID = IDs[i]
  lines(dPD7[dPD7$ID == cID, "TIME"], dPD7[dPD7$ID == cID, "CONC"])
}
###

# Figure 7.3
dev.new()
plot(1, 1, type="n", log="x", xlim=c(1, 1000), ylim=c(5, 35), xlab="Concentration", ylab="Response")
abline(h=seq(5,35,5), lty=3)
abline(v=c(2:9, seq(10, 90, 10), seq(100, 1000, 100)), lty=3)
for (i in 1:nID) {
  cID = IDs[i]
  lines(dPD7[dPD7$ID == cID, "CONC"], dPD7[dPD7$ID == cID, "DV"])  
}
###

fPD7 = function(THETA)
{
  Kin  = THETA[1]
  Kout = THETA[2]
  Emax = THETA[3]
  EC50 = THETA[4]

  y = vector()
  for (i in 1:nID) {
    cID <<- IDs[i] # Referencing within fPD7de
    cTIME = c(0, dPD7[dPD7$ID == cID, "TIME"]) # TIME should include 0.
    cy = lsoda(c(0, Kin/Kout), cTIME, fPD7de, c(Kin=Kin, Kout=Kout, Emax=Emax, EC50=EC50))
    y = c(y, cy[-1, "2"])
  }  
  return(y)
}
y = fPD7(c(Kin, Kout, Emax, EC50)) ; y
length(y)
length(dPD7[,"DV"])

r1 = nlr(fPD7, dPD7, c("Kin", "Kout", "Emax", "EC50"), IE=c(30, 1, 4, 42))
r1
e$r$value + e$nRec*log(2*pi) # -2LL, compare with Phoenix results in Overall table

wnl5(fPD7, dPD7, c("Kin", "Kout", "Emax", "EC50"), IE=c(30, 1, 4, 42)) # to see AIC, WRSS in Table 7.1. Compare with nlr results


