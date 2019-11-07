# PD 52
require(wnl)
setwd("D:/Rt/PD")

dPD52 = read.csv("PD52.csv")
colnames(dPD52) = c("TIME", "DV", "ID", "DOSE")

IDs = unique(dPD52[,"ID"]) ; IDs
nID = length(IDs) ; nID
AMTs = unique(dPD52[,"DOSE"]) ; AMTs

require(deSolve)
# First order biophse model
fPD52ade = function(t, y, p)
{
  Cp    = AMTs[i]*p["Kp"]*t*exp(-p["Kp"]*t)
  Stim  = p["Smax"]*Cp^p["n"]/(p["SD50"]^p["n"] + Cp^p["n"])
  Kut   = p["Kout"]/(p["Km"] + y[1])
  dy1dt = Stim - Kut*y[1]
  return(list(c(dy1dt)))
}

fPD52a = function(THETA)
{
  Kp   = THETA[1]
  n    = THETA[2]
  Km   = THETA[3]/1e3
  Kout = THETA[4]
  Smax = THETA[5]
  SD50 = THETA[6]

  y = vector()
  for (i in 1:nID) {
    i <<- i
    cID   = IDs[i]
    cTIME = unique(c(0, dPD52[dPD52$ID == cID, "TIME"]))
    iTIME = cTIME %in% dPD52[dPD52$ID == cID, "TIME"]
    cy = lsoda(c(0), cTIME, fPD52ade, c(Kp=Kp, n=n, Km=Km, Kout=Kout, Smax=Smax, SD50=SD50))
    y = c(y, cy[iTIME, "1"])
  }
  return(y)
}
y0a = fPD52a(c(5.93, 1.68, 1, 30.5, 244, 0.9842)) ; y0a
length(y0a)

r1 = nlr(fPD52a, dPD52, c("Kp", "n", "Km", "Kout", "Smax", "SD50"), c(10, 2, 0.001, 30, 240, 5)) ; r1

# Bolus approximation model
fPD52bde = function(t, y, p)
{
  Cp    = AMTs[i]*exp(-p["Kp"]*t)
  Stim  = p["Smax"]*Cp^p["n"]/(p["SD50"]^p["n"] + Cp^p["n"])
  Kut   = p["Kout"]/(p["Km"] + y[1])
  dy1dt = Stim - Kut*y[1]
  return(list(c(dy1dt)))
}

fPD52b = function(THETA)
{
  Kp   = THETA[1]
  n    = THETA[2]
  Km   = THETA[3]/1e3
  Kout = THETA[4]
  Smax = THETA[5]
  SD50 = THETA[6]

  y = vector()
  for (i in 1:nID) {
    i <<- i
    cID = IDs[i]
    cTIME = unique(c(0, dPD52[dPD52$ID == cID, "TIME"]))
    iTIME = cTIME %in% dPD52[dPD52$ID == cID, "TIME"]
    cy = lsoda(c(0), cTIME, fPD52bde, c(Kp=Kp, n=n, Km=Km, Kout=Kout, Smax=Smax, SD50=SD50))
    y = c(y, cy[iTIME, "1"])
  }
  return(y)
}
y0b = fPD52b(c(2.49348, 2.02, 1, 32.52, 195, 1.556)) ; y0b
length(y0b)

r2 = nlr(fPD52b, dPD52, c("Kp", "n", "Km", "Kout", "Smax", "SD50"), c(6, 1.7, 0.001, 30, 240, 1)) ; r2

# Figure 52.1
plot(0, 0, type="n", xlim=c(0, 3.5), ylim=c(0, 80), xlab="Time (h)", ylab="Locomotor activity")
abline(h=1:8*10, lty=3)
for (i in 1:nID) {
  cID = IDs[i]
  points(dPD52[dPD52$ID == cID, "TIME"], dPD52[dPD52$ID == cID, "DV"], pch=17-i, col=c("red", "blue")[i])
  lines(dPD52[dPD52$ID == cID, "TIME"], y0a[dPD52$ID == cID])
  lines(dPD52[dPD52$ID == cID, "TIME"], y0b[dPD52$ID == cID], lty=2)
}

# Refined model
fPD52cde = function(t, y, p)
{
  Cp    = AMTs[i]*exp(-p["Kp"]*t)
  Stim  = p["Smax"]*Cp^p["n"]/(p["SD50"]^p["n"] + Cp^p["n"])
  Kut   = p["Kout"]/(1e-8 + y[1])
  dy1dt = Stim - Kut*y[1]
  return(list(c(dy1dt)))
}

fPD52c = function(THETA)
{
  Kp   = THETA[1]
  n    = THETA[2]
  Kout = THETA[3]
  Smax = THETA[4]
  SD50 = THETA[5]

  y = vector()
  for (i in 1:nID) {
    i <<- i
    cID = IDs[i]
    cTIME = unique(c(0, dPD52[dPD52$ID == cID, "TIME"]))
    iTIME = cTIME %in% dPD52[dPD52$ID == cID, "TIME"]
    cy = lsoda(c(1e-8), cTIME, fPD52cde, c(Kp=Kp, n=n, Kout=Kout, Smax=Smax, SD50=SD50))
    y = c(y, cy[iTIME, "1"])
  }
  return(y)
}
y0c = fPD52c(c(2.49348, 2.02, 32.52, 195, 1.556)) ; y0c
length(y0c)

r3 = nlr(fPD52c, dPD52, c("Kp", "n", "Kout", "Smax", "SD50"), c(6, 1.7, 30, 240, 1)) ; r3
dx(r3)

# Alternative model 2
fPD52dde = function(t, y, p)
{
  Cp    = AMTs[i]*exp(-p["Kp"]*t)
  Stim  = 1 + p["Smax"]*Cp^p["n"]/(p["SD50"]^p["n"] + Cp^p["n"])
  dy1dt = p["Kin"]*Stim - p["Kout"]*y[1]
  return(list(c(dy1dt)))
}

fPD52d = function(THETA)
{
  Kp   = THETA[1]
  n    = THETA[2]
  Kin  = THETA[3]
  Kout = THETA[4]
  Smax = THETA[5]
  SD50 = THETA[6]

  y = vector()
  for (i in 1:nID) {
    i <<- i
    cID = IDs[i]
    cTIME = unique(c(0, dPD52[dPD52$ID == cID, "TIME"]))
    iTIME = cTIME %in% dPD52[dPD52$ID == cID, "TIME"]
    cy = lsoda(c(Kin/Kout), cTIME, fPD52dde, c(Kp=Kp, n=n, Kin=Kin, Kout=Kout, Smax=Smax, SD50=SD50))
    y = c(y, cy[iTIME, "1"])
  }
  return(y)
}
y0d = fPD52d(c(2.49348, 2.02, 0.1, 32.52, 195, 1.556)) ; y0d
length(y0d)

r4 = nlr(fPD52d, dPD52, c("Kp", "n", "Kin", "Kout", "Smax", "SD50"), c(6, 1.7, 0.1, 30, 240, 1)) ; r4
dx(r4)

