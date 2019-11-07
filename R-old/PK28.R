setwd("D:/Rt/Gab")
require(wnl)
dPK28 = read.csv("PK28.csv", skip=1)
colnames(dPK28) = c("TIME", "ID", "DV", "Func", "Species") ; dPK28

IDs = unique(dPK28[,"ID"])
nID = length(IDs)
BWTs = c(0.023, 0.25, 70)
AMTs = c(25, 500, 100000)

fPK28 = function(THETA)
{
  a = THETA[1]
  b = THETA[2]
  cc = THETA[3]
  d = THETA[4]
  
  Cp = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = AMTs[i]
    cBWT = BWTs[i]
    cCL = a*cBWT^b
    cV  = cc*cBWT^d
    cTIME = dPK28[dPK28$ID == cID, "TIME"]
    Cp = c(Cp, cAMT/cV*exp(-cCL/cV*cTIME))
  }
  return(Cp)
} 
fPK28(c(0.319142, 0.636712, 3.07666, 1.03094))

nlr(fPK28, dPK28, pNames=c("a", "b", "c", "d"), IE=c(0.5, 0.8, 3.5, 1.52), Error="P")

# Figure 28.2, p 613
plot(0, 0.01, type="n", xlim=c(0, 20), ylim=c(0.01, 1.00), log="y", xlab="Kallynochrons (h/BW^(1-b))", ylab="Conc*BW/Dose")
for (i in 1:nID) {
  cID = IDs[i]
  x = dPK28[dPK28$ID == cID, "TIME"] / BWTs[i]^(1 - 0.636712)
  y = dPK28[dPK28$ID == cID, "DV"] * BWTs[i] / AMTs[i]  # Erratum in figure 28.2
  lines(x, y)
}
##

