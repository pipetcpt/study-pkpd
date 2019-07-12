# PD 46
require(wnl)
require(deSolve)
setwd("D:/Rt/PD")


# mRNA model
dPD46a = read.csv("PD46_mRNA.csv", header=FALSE, skip=2, na="Missing")
colnames(dPD46a) = c("TIME", "DV", "ID", "Cp")
dPD46a = dPD46a[!is.na(dPD46a$DV),]
dPD46a
dim(dPD46a)

IDs = unique(dPD46a[,"ID"])
nID = length(IDs)

AMTs = c(175, 3000)
Ka = 13.6 # /hr
V  = 2.703
Cl = 2.068
Ke = Cl/V
n  = 6
R0 = 1

fPD46ade = function(t, y, p)
{
  Cp = AMTs[i]/V*Ka/(Ka - Ke)*(exp(-Ke*t) - exp(-Ka*t))
  Kin = p["Kout"]
  Stim = 1 + p["Smax"]*Cp^n/(p["SC50"]^n + Cp^n)
  dy1dt = Kin*Stim - p["Kout"]*y[1]*y[5]
  dy2dt = p["Kout"]*y[1] - p["Kout"]*y[2]
  dy3dt = p["Kout"]*y[2] - p["Kout"]*y[3]
  dy4dt = p["Ktol"]*y[1] - p["Ktol"]*y[4]
  dy5dt = p["Ktol"]*y[4] - p["Ktol"]*y[5]
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt, dy5dt)))
}

fPD46a = function(THETA)
{
  Ktol = THETA[1]
  Kout = THETA[2]
  Smax = THETA[3]
  SC50 = THETA[4]

  y = vector()
  for (i in 1:nID) {
    i <<- i
    cID = IDs[i]
    cTIME = unique(c(0, dPD46a[dPD46a$ID == cID, "TIME"])) # 0 now added
    iTIME = cTIME %in% dPD46a[dPD46a$ID == cID, "TIME"]
    cy = lsoda(c(R0, R0, R0, R0, R0), cTIME, fPD46ade, c(Ktol=Ktol, Kout=Kout, Smax=Smax, SC50=SC50))
    y = c(y, cy[iTIME, "3"])
  }
  return(y)
}
y1 = fPD46a(c(0.198542, 1.01996, 79.211, 31.5192)) ; y1
length(y1)

r1 = nlr(fPD46a, dPD46a, c("Ktol", "Kout", "Smax", "SC50"), c(0.26, 5, 29, 30)) ; r1

# Figure 46.2 : note the DATA difference
Ktol = 0.179
Kout = 1.01996
Smax = 71.2
SC50 = 30.86
y0 = matrix(nrow=0, ncol=6)
for (i in 1:nID) {
  i <<- i
  cID = IDs[i]
  cTIME <<- seq(0, 42, length=201)
  cy = lsoda(c(R0, R0, R0, R0, R0), cTIME, fPD46ade, c(Ktol=Ktol, Kout=Kout, Smax=Smax, SC50=SC50))
  y0 = rbind(y0, cy)
} ; y0
iTime = cbind(1:201, 202:402)

dev.new(width=7, height=5)
plot(1, 1, type="n", log="y", xlim=c(0, 42), ylim=c(0.1, 100), xlab="Time (h)", ylab="Fold mRNA")
abline(h=as.vector(outer(1:9, 10^(-1:3))), lty=3)
for (i in 1:nID) {
  cID = IDs[i]
  points(dPD46a[dPD46a$ID == cID, "TIME"], dPD46a[dPD46a$ID == cID, "DV"], pch=14+i, col=c("red", "blue")[i])
  lines(cTIME, y0[iTime[,i],"3"], col=c("red", "blue")[i])
}

# Protein model
dPD46b = read.csv("PD46_protein.csv", header=FALSE, skip=2, na="Missing")
colnames(dPD46b) = c("TIME", "DV", "ID", "Cp", "mRNA")
dPD46b = dPD46b[!is.na(dPD46b$DV),]
dPD46b
dim(dPD46b)

n = 10
Kin  = Kout
fPD46bde = function(t, y, p)
{
  Cp = AMTs[i]/V*Ka/(Ka - Ke)*(exp(-Ke*t) - exp(-Ka*t))
  Stim1 = 1 + Smax*Cp^n/(SC50^n + Cp^n)

  dy1dt = Kin*Stim1 - Kout*y[1]*y[5]
  dy2dt = Kout*y[1] - Kout*y[2]
  dy3dt = Kout*y[2] - Kout*y[3]
  dy4dt = Ktol*y[1] - Ktol*y[4]
  dy5dt = Ktol*y[4] - Ktol*y[5]

  Kpro = p["Kdeg"]
  Stim2 = y[3]^p["a"]

  dy6dt = Kpro*Stim2 - p["Kdeg"]*y[6]
  dy7dt = p["Kdeg"]*y[6] - p["Kdeg"]*y[7]
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt, dy5dt, dy6dt, dy7dt)))
}

fPD46b = function(THETA)
{
  a    = THETA[1]
  Kdeg = THETA[2]
  y = vector()
  for (i in 1:nID) {
    i <<- i
    cID = IDs[i]
    cTIME = unique(c(0, dPD46b[dPD46b$ID == cID, "TIME"])) # 0 now added
    iTIME = cTIME %in% dPD46b[dPD46b$ID == cID, "TIME"]
    cy = lsoda(c(R0, R0, R0, R0, R0, R0, R0), cTIME, fPD46bde, c(a=a, Kdeg=Kdeg), hmax=1)
    y = c(y, cy[iTIME, "7"])
  }
  return(y)
}
y2 = fPD46b(c(0.779207, 0.249095)) ; y2
length(y2)

r2 = nlr(fPD46b, dPD46b, c("a", "Kdeg"), c(3, 0.3)) ; r2


# Figure 46.4
a    = 0.779207
Kdeg = 0.249095
y0b = matrix(nrow=0, ncol=8)
for (i in 1:nID) {
  i <<- i
  cID = IDs[i]
  cTIME <<- seq(0, 42, length=201)
  cy = lsoda(c(R0, R0, R0, R0, R0, R0, R0), cTIME, fPD46bde, c(a=a, Kdeg=Kdeg))
  y0b = rbind(y0b, cy)
} ; y0b
iTime = cbind(1:201, 202:402)

dev.new(width=7, height=5)
plot(1, 1, type="n", log="y", xlim=c(0, 42), ylim=c(0.5, 10), xlab="Time (h)", ylab="Fold protein")
abline(h=as.vector(outer(1:9, 10^(-1:2))), lty=3)
for (i in 1:nID) {
  cID = IDs[i]
  points(dPD46b[dPD46b$ID == cID, "TIME"], dPD46b[dPD46b$ID == cID, "DV"], pch=14+i, col=c("red", "blue")[i])
  lines(cTIME, y0b[iTime[,i],"7"], col=c("red", "blue")[i])
}

