# PD 43
require(wnl)
setwd("D:/Rt/PD")

dPD43 = read.csv("PD43.csv", header=FALSE, skip=2, na="Missing")
colnames(dPD43) = c("TIME", "DV", "ID")
dPD43b = dPD43[!is.na(dPD43$DV),]
dim(dPD43b) ; dPD43b

IDs = unique(dPD43b[,"ID"])
nID = length(IDs)

AMTs = c(0, 0.5, 2, 10)

Ka   = 4.0488
V    = 2.561
Ke   = 0.1439
Tlag = 0.8454
Pow  = 2/3

require(deSolve)
fPD43de = function(t, y, p)
{
  if (t > Tlag) {
    Cp = AMTs[i]/V*Ka/(Ka - Ke)*(exp(-Ke*(t - Tlag)) - exp(-Ka*(t - Tlag)))
  } else {
    Cp = 0
  }
  Inh = p["Kmax"]*Cp/(p["IC50"] + Cp)
  E11 = max(y[1], 0)^Pow
  E21 = max(y[2], 0)^Pow
  E31 = max(y[3], 0)^Pow
  E41 = max(y[4], 0)^Pow

  dy1dt = (p["Kgro"] - Inh*p["Kout"])*E11
  dy2dt = Inh*p["Kout"]*E11 - p["Kout"]*E21
  dy3dt = p["Kout"]*E21 - p["Kout"]*E31
  dy4dt = p["Kout"]*E31 - p["Kout"]*E41

  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt)))
}

Kgro = 0.0363923
Kout = 0.171686
Kmax = 14.145
IC50 = 1.60506
R0   = 357.576
y0 = matrix(nrow=0, ncol=5)
for (i in 1:nID) {
  i <<- i
  cID = IDs[i]
  cTIME = dPD43b[dPD43b$ID == cID, "TIME"] # O already included
  cy = lsoda(c(R0, 0, 0, 0), cTIME, fPD43de, c(Kgro=Kgro, Kout=Kout, Kmax=Kmax, IC50=IC50), hmax=1)
  y0 = rbind(y0, cy)
} ; y0

fPD43 = function(THETA)
{
  Kgro = THETA[1]
  Kout = THETA[2]
  Kmax = THETA[3]
  IC50 = THETA[4]
  R0   = THETA[5]

  y = vector()
  for (i in 1:nID) {
    i <<- i
    cID = IDs[i]
    cTIME = dPD43b[dPD43b$ID == cID, "TIME"] # O already included
    cy = lsoda(c(R0, 0, 0, 0), cTIME, fPD43de, c(Kgro=Kgro, Kout=Kout, Kmax=Kmax, IC50=IC50), hmax=1)
    y = c(y, rowSums(cy[,2:5]))
  }
  return(y)
}
y1 = fPD43(c(Kgro, Kout, Kmax, IC50, R0)) ; y1
length(y1)

r1 = nlr(fPD43, dPD43b, c("Kgro", "Kout", "Kmax", "IC50", "R0"), c(0.04, 0.2, 15, 2, 360),
         LB=c(0.01, 0.01, 1, 1, 300),
         UB=c(10, 10, 50, 100, 400)) ; r1
e$PE
e$r

r2 = nlr(fPD43, dPD43b, c("Kgro", "Kout", "Kmax", "IC50", "R0"), c(0.1, 0.1, 5, 100, 350)) ; r2
e$PE

w1 = wnl5(fPD43, dPD43b, c("Kgro", "Kout", "Kmax", "IC50", "R0"), c(0.01, 0.1, 5, 300, 350)) ; w1

# Figure 43.1
plot(0, 0, xlim=c(0, 700), ylim=c(0, 1600), xlab="Time (h)", ylab="Tumor volume")
abline(h=seq(200, 1600, 200), lty=3)
for (i in 1:nID) {
  cID = IDs[i]
  points(dPD43b[dPD43b$ID == cID, "TIME"], dPD43b[dPD43b$ID == cID, "DV"], pch=14+i)
  lines(dPD43b[dPD43b$ID == cID, "TIME"], y1[dPD43b$ID == cID])  
}



