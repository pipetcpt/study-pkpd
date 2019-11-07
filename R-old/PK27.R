setwd("D:/Rt/Gab/")
require(wnl)
dPK27c = read.csv("PK27_data3.csv", as.is=TRUE)
colnames(dPK27c) = c("TIME", "AMT", "DV", "ID") ; dPK27c
dPK27c = dPK27c[dPK27c$DV != "Missing",]
dPK27c[,"DV"] = as.numeric(dPK27c[,"DV"]) ; dPK27c

dPK27 = dPK27c[dPK27c$ID %in% 1:4,] ; dPK27
AMTs = dPK27[dPK27$TIME==0,"AMT"] ; AMTs

require(deSolve)
PKde = function(t, y, p)
{
  dy1dt = (-p["CL"]*y[1] - p["CLd"]*(y[1] - y[2]))/0.05 - p["Kon"]*y[1]*y[3] + p["Koff"]*y[4] # Ligand, Erratum in eq 27:1
  dy2dt = p["CLd"]*(y[1] - y[2])/p["Vt"]                                  # eq 27:1, peripheral for Ligand
  dy3dt = p["R0"]*p["Kout"] - p["Kout"]*y[3] - p["Kon"]*y[1]*y[3] + p["Koff"]*y[4] # eq 27:2, Receptor 
  dy4dt = p["Kon"]*y[1]*y[3] - p["Koff"]*y[4] - p["KeRL"]*y[4]            # eq 27:3, R + L complex
  
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt)))
}

# Figure 27.5, p606 -> Figure erratum: Subject 1 last point missing at data file
plot(0, 0.001, type="n", xlim=c(0, 600), ylim=c(0.001, 1000), xlab="Time (h)", ylab="Concentration (mg/L)", log="y")

IDs = unique(dPK27[,"ID"])
nID = length(IDs)

TIME = unique(dPK27[,"TIME"]) ; TIME
nTIME = length(TIME)

y = vector()
for (i in 1:nID) {
  cID = IDs[i]
  cAMT = AMTs[cID]
  cy = lsoda(y=c(cAMT/0.05, 0, 12, 0), times=TIME, func=PKde, parms=c(Vt=0.1, CL=0.001, CLd=0.003, R0=12, Kon=0.102, Koff=0.0009, Kout=0.00887782, KeRL=0.00187727))
  points(dPK27[dPK27$ID == cID, "TIME"], dPK27[dPK27$ID == cID, "DV"], pch=19)
  points(dPK27c[dPK27c$ID == cID+4, "TIME"], dPK27c[dPK27c$ID == cID+4, "DV"], pch=15, col="blue")
  points(dPK27c[dPK27c$ID == cID+8, "TIME"], dPK27c[dPK27c$ID == cID+8, "DV"], pch=17, col="red")
  lines(TIME, cy[,"1"])
  lines(TIME, cy[,"3"], lty=2)
  lines(TIME, cy[,"4"], lty=3)
  y = rbind(y, cy)
} ; y

y[-10,c("time","1")]
y[c(-13,-22,-23,-24,-32,-33,-34,-35,-36,-37),c("time","3")]
y[-38,c("time","4")]
###

fPK27a = function(THETA)
{
  Vt   = THETA[1]
  CL   = THETA[2]
  CLd  = THETA[3]
  R0   = THETA[4]
  Kon  = THETA[5]
  Koff = THETA[6]
  Kout = THETA[7]
  KeRL = THETA[8]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = AMTs[cID]
    cy = lsoda(y=c(cAMT/0.05, 0, R0, 0), times=TIME, func=PKde, parms=c(Vt=Vt, CL=CL, CLd=CLd, R0=R0, Kon=Kon, Koff=Koff, Kout=Kout, KeRL=KeRL))
    y = rbind(y, cy)
  }
  
  return(y[-10,"1"])
}
fPK27a(c(0.1, 0.001, 0.003, 12, 0.102, 0.0009,0.00887782, 0.00187727))
fPK27a(c(0.1, 0.001, 0.003, 12, 0.091, 0.011, 0.0089, 0.003))

nlr(fPK27a, dPK27, pNames=c("Vt", "CL", "CLd", "R0", "Kon", "Koff", "Kout", "KeRL"), 
    IE=c(0.1, 0.001, 0.003, 12, 0.1, 0.001, 0.01, 0.002), 
    LB=c(0.05, 0.0005, 0.001, 6, 0.05, 0.0005, 0.005, 0.001),
    UB=c(0.2, 0.002, 0.005, 18, 0.2, 0.002, 0.02, 0.003), Error="P")


## Dataset 2
dPK27b = dPK27c[dPK27c$ID < 9,] ; dPK27b
fPK27b = function(THETA)
{
  Vt   = THETA[1]
  CL   = THETA[2]
  CLd  = THETA[3]
  R0   = THETA[4]
  Kon  = THETA[5]
  Koff = THETA[6]
  Kout = THETA[7]
  KeRL = THETA[8]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = AMTs[cID]
    cy = lsoda(y=c(cAMT/0.05, 0, R0, 0), times=TIME, func=PKde, parms=c(Vt=Vt, CL=CL, CLd=CLd, R0=R0, Kon=Kon, Koff=Koff, Kout=Kout, KeRL=KeRL))
    y = rbind(y, cy)
  }
  
  return(c(y[-10,"1"], y[c(-13,-22,-23,-24,-32,-33,-34,-35,-36,-37),"3"]))
}
y = fPK27b(c(0.1, 0.001, 0.003, 12, 0.102, 0.0009,0.00887782, 0.00187727))
length(y)
length(dPK27b[,"DV"])

nlr(fPK27b, dPK27b, pNames=c("Vt", "CL", "CLd", "R0", "Kon", "Koff", "Kout", "KeRL"), 
    IE=c(0.1, 0.001, 0.003, 12, 0.1, 0.001, 0.01, 0.002), 
    LB=c(0.05, 0.0005, 0.001, 6, 0.05, 0.0005, 0.005, 0.001),
    UB=c(0.2, 0.002, 0.005, 18, 0.2, 0.002, 0.02, 0.003), Error="P")

## Dataset 3
fPK27c = function(THETA)
{
  Vt   = THETA[1]
  CL   = THETA[2]
  CLd  = THETA[3]
  R0   = THETA[4]
  Kon  = THETA[5]
  Koff = THETA[6]
  Kout = THETA[7]
  KeRL = THETA[8]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = AMTs[cID]
    cy = lsoda(y=c(cAMT/0.05, 0, R0, 0), times=TIME, func=PKde, parms=c(Vt=Vt, CL=CL, CLd=CLd, R0=R0, Kon=Kon, Koff=Koff, Kout=Kout, KeRL=KeRL))
    y = rbind(y, cy)
  }
  
  return(c(y[-10,"1"], y[c(-13,-22,-23,-24,-32,-33,-34,-35,-36,-37),"3"], y[-38,"4"]))
}
y = fPK27c(c(0.1, 0.001, 0.003, 12, 0.102, 0.0009,0.00887782, 0.00187727))
length(y)
length(dPK27c[,"DV"])

nlr(fPK27c, dPK27c, pNames=c("Vt", "CL", "CLd", "R0", "Kon", "Koff", "Kout", "KeRL"), 
    IE=c(0.1, 0.001, 0.003, 12, 0.1, 0.001, 0.01, 0.002), 
    LB=c(0.05, 0.0005, 0.001, 6, 0.05, 0.0005, 0.005, 0.001),
    UB=c(0.2, 0.002, 0.005, 18, 0.2, 0.002, 0.02, 0.003), Error="P")

