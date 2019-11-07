setwd("D:/Rt/Gab/")
require(wnl)
dPK29 = read.csv("PK29.csv", skip=1)
colnames(dPK29) = c("TIME", "DV", "ID", "Species", "BWT") ; dPK29

IDs = unique(dPK29[,"ID"])
nID = length(IDs)
BWTs = unique(dPK29[,"BWT"])
AMTs = c(10, 125, 200, 6000, 12000)

fPK29 = function(THETA)
{
  a = THETA[1]
  b = THETA[2]
  cc = THETA[3]
  d = THETA[4]
  e = THETA[5]
  g = THETA[6]
  
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = AMTs[i]
    cBWT = BWTs[i]
    cTIME = dPK29[dPK29$ID == cID, "TIME"]
    CL  = a*cBWT^b
    Vc  = cc*cBWT^d
    CLd = g*cBWT^b
    Vt  = e*cBWT^d 
    k12 = CLd/Vc
    k21 = CLd/Vt
    k10 = CL/Vc 
    beta = 0.5*(k12 + k21 + k10 - sqrt((k12 + k21 + k10)^2 - 4*k21*k10))
    alpha = k21*k10/beta   
    Cp = cAMT/Vc*((k21 - alpha)/(beta - alpha)*exp(-alpha*cTIME) + (k21 - beta)/(alpha - beta)*exp(-beta*cTIME))
    y = c(y, Cp)
  }
  return(y)
} 
fPK29(c(0.021572, 0.76765, 0.1, 1.2, 0.5158, 0.0658))

nlr(fPK29, dPK29, pNames=c("a", "b", "c", "d", "e", "g"), IE=c(0.1, 1, 1, 1.2, 0.8, 0.1), Error="POIS")


# Figure 29.3 Lower, p 619
plot(0, 0.0001, type="n", xlim=c(0, 400), ylim=c(0.0001, 10), log="y", xlab="Apolysichrons (T*BW^(b-d))", ylab="Conc*Dose^(-d)")
for (i in 1:nID) {
  cID = IDs[i]
  x = dPK29[dPK29$ID == cID, "TIME"] * BWTs[i]^(0.76765 - 1.2)
  y = dPK29[dPK29$ID == cID, "DV"] * BWTs[i]^1.2 / AMTs[i] # Erratum in Figure 29.3
  points(x, y, pch=14+i)
}
##

