# PD 6
require(wnl)
setwd("D:/Rt/PD")

dPD6 = read.csv("PD6.csv", header=FALSE, skip=2)
colnames(dPD6) = c("TIME", "DV", "ID", "AMT")
dim(dPD6)
head(dPD6)

IDs = unique(dPD6[,"ID"])
nID = length(IDs)

V   = 5.205
Ke  = 0.456
AMTs = c(10750, 43000, 172000)

require(deSolve)
# Turn over model (Indirect response model)
fPD6ade = function(t, y, p)
{
  dy1dt = -Ke*y[1]
  dy2dt = p["Kin"]*(1 + y[1]/p["EC50"]) - p["Kout"]*y[2]
  return(list(c(dy1dt, dy2dt)))
}

# Figure 6.1
plot(0, 0, type="n", xlim=c(0, 8), ylim=c(0, 800), xlab="Time (h)", ylab="Response")
abline(h=100*(1:8), lty=3)
Kin   = 233.018
Kout  = 5.56636
EC50  = 1627.99
Tdose = 0.691629
y = vector()
for (i in 1:nID) {
  cID = IDs[i]
  Times = sort(c(dPD6[dPD6$ID == cID, "TIME"], Tdose))   # buggy
  iTdose = sum(Tdose > dPD6[dPD6$ID == cID, "TIME"]) + 1 # buggy
  EventDat = data.frame(var = "y1",
                        time = Tdose,
                        value = AMTs[i]/V,
                        method = "add")
  cy = lsoda(y=c(y1=0, y2=Kin/Kout), times=Times, func=fPD6ade, events=list(data=EventDat),
             parms=c(Kin=Kin, Kout=Kout, EC50=EC50))
  points(dPD6[dPD6$ID == cID, "TIME"], dPD6[dPD6$ID == cID, "DV"], pch=16)
  lines(Times[-iTdose], cy[-iTdose, "y2"])
  y = c(y, cy[, "y2"])
} ; y
length(y)
###

fPD6a = function(THETA)
{
  Kin   = THETA[1]
  Kout  = THETA[2]
  EC50  = THETA[3]
  Tdose = THETA[4]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    Times = sort(c(dPD6[dPD6$ID == cID, "TIME"], Tdose))
    iTdose = sum(Tdose > dPD6[dPD6$ID == cID, "TIME"]) + 1

    EventDat = data.frame(var = "y1",
                          time = Tdose,
                          value = AMTs[i]/V,
                          method = "add")
    cy = lsoda(y=c(y1=0, y2=Kin/Kout), times=Times, func=fPD6ade, events=list(data=EventDat),
               parms=c(Kin=Kin, Kout=Kout, EC50=EC50))
              
    y = c(y, cy[-iTdose, "y2"])
  }
  return(y)
}
y = fPD6a(c(Kin, Kout, EC50, Tdose)) ; y
length(y)
length(dPD6[,"DV"])

r1 = nlr(fPD6a, dPD6, pNames=c("Kin", "Kout", "EC50", "Tdose"), IE=c(340, 5, 2285, 0.8))
r1

# Link model (effect compartment model)
fPD6bde = function(t, y, p)
{
  dy1dt = -Ke*y[1]
  dy2dt = p["Ke0"]*(y[1] - y[2])
  return(list(c(dy1dt, dy2dt)))
}

fPD6b = function(THETA)
{
  Ke0   = THETA[1]
  A     = THETA[2]
  E0    = THETA[3]
  Tdose = THETA[4]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    Times = sort(c(dPD6[dPD6$ID == cID, "TIME"], Tdose))
    iTdose = sum(Tdose > dPD6[dPD6$ID == cID, "TIME"]) + 1

    EventDat = data.frame(var = "y1",
                          time = Tdose,
                          value = AMTs[i]/V,
                          method = "add")
    cy = lsoda(y=c(y1=0, y2=E0), times=Times, func=fPD6bde, events=list(data=EventDat),
               parms=c(Ke0=Ke0))
              
    y = c(y, E0 + A*cy[-iTdose, "y2"])
  }
  return(y)
}
y = fPD6b(c(5.56635, 0.0257139, 41.8618, 0.691629)) ; y
length(y)
length(dPD6[,"DV"])

r2 = nlr(fPD6b, dPD6, pNames=c("Ke0", "A", "E0", "Tdose"), IE=c(0.8, 0.019, 42.2, 0.8))
r2
