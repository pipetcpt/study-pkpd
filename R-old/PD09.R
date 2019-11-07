# PD 9
require(wnl)
setwd("D:/Rt/PD")

dPD9 = read.csv("PD9.csv", header=FALSE, skip=2)
colnames(dPD9) = c("TIME", "DV", "ID")
dim(dPD9)
dPD9

IDs = unique(dPD9[,"ID"])
nID = length(IDs)
AMTs = c(5, 25)
nDose = 4 
II = 24 # hr, Interdose Interval

Ka = 1.1 # /hr
V  = 5
Ke = 0.128 # /hr

require(deSolve)
fPD8de = function(t, y, p)
{
  dy1dt = -Ka*y[1]
  dy2dt = Ka*y[1] - Ke*y[2]
  Conc  = y[2]/V
  Inh   = 1 - Conc^p["n"]/(p["IC50"]^p["n"] + Conc^p["n"]) 
  dy3dt = Inh*p["Kin"] - p["Kout"]*y[3]
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

# Figure 9.1
plot(0, 0, type="n", xlim=c(0, 100), ylim=c(0, 80), xlab="Time (h)", ylab="Response")
abline(h=seq(10, 80, 10), lty=3)

Kin  = 8.80231
Kout = 0.105859
IC50 = 0.244842
n    = 1.38354

for (i in 1:nID) {
  cID = IDs[i]
  cTIME = sort(unique(c(dPD9[dPD9$ID == cID, "TIME"], seq(0, 72, 0.2))))
  iTime = cTIME %in% dPD9[dPD9$ID == cID, "TIME"]
  EventDat = data.frame(var = rep("y1", nDose),
                        time = seq(0, 72, II),
                        value = rep(AMTs[i], nDose),
                        method = rep("add", nDose))
  cy = lsoda(y=c(y1=0, y2=0, y3=Kin/Kout), times=cTIME, func=fPD8de, events=list(data=EventDat),
             parms=c(Kin=Kin, Kout=Kout, IC50=IC50, n=n))
  dPD9[dPD9$ID == cID, "CONC"] = cy[iTime, "y2"]
  dPD9[dPD9$ID == cID, "PRED"] = cy[iTime, "y3"]
  points(dPD9[dPD9$ID == cID, "TIME"], dPD9[dPD9$ID == cID, "DV"], pch=14 + i)
  lines(cTIME, cy[,"y3"]) 
} ; dPD9
###

fPD9 = function(THETA)
{
  Kin  = THETA[1]
  Kout = THETA[2]
  IC50 = THETA[3]
  n    = THETA[4]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cTIME = sort(unique(c(dPD9[dPD9$ID == cID, "TIME"], seq(0, 72, II))))
    iTime = cTIME %in% dPD9[dPD9$ID == cID, "TIME"]
    EventDat = data.frame(var = rep("y1", nDose),
                          time = seq(0, 72, II),
                          value = rep(AMTs[i], nDose),
                          method = rep("add", nDose))
    cy = lsoda(y=c(y1=0, y2=0, y3=Kin/Kout), times=cTIME, func=fPD8de, events=list(data=EventDat),
               parms=c(Kin=Kin, Kout=Kout, IC50=IC50, n=n))
    y = c(y, cy[iTime, "y3"])
  }
  return(y)
}
y = fPD9(c(Kin, Kout, IC50, n))
y
length(y)

r1 = nlr(fPD9, dPD9, c("Kin", "Kout", "IC50", "n"), c(4.8, 0.06, 0.25, 2))
r1

e$nRec*log(2*pi) + e$r$value # -2LL, compare with Phoenix result
e$nRec*log(2*pi) + e$r$value + e$nPara*log(e$nRec) # BIC, compare with Phoenix result
