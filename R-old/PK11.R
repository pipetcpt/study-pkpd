require(wnl)
dPK11pk = read.csv("data-old/PK11pk.csv", skip=1)
colnames(dPK11pk) = c("TIME", "DV") ; dPK11pk

dPK11pd = read.csv("PK11pd.csv", skip=1)
colnames(dPK11pd) = c("TIME", "DV") ; dPK11pd

DosingTimes = c(0, seq(96, 192, by=8)) ; DosingTimes
nDose = length(DosingTimes)
AMT = rep(400, nDose) ; AMT
DoseHist = cbind(TIME=DosingTimes, AMT=AMT) ; DoseHist
nRec = nrow(dPK11pk)

# Errata in eq 11:1 and 11:2   one of (alpha - Ka)s -> (beta - Ka)

## fitting PK model
fPK11pk = function(THETA)
{
  Vc  = THETA[1]
  Ka  = THETA[2]
  k21 = THETA[3]
  a   = THETA[4] # alpha
  b   = THETA[5] # beta
  TL  = THETA[6] # Tlag

  T1 = e$DATA[, "TIME"]

  Res = matrix(rep(0, nRec*nDose), nrow=nRec, ncol=nDose)
  for (i in 1:nDose) {
    dTime = DoseHist[i,"TIME"]
    cAMT = DoseHist[i,"AMT"]
    cTime = T1 > (dTime + TL)
    Res[cTime,i] = Ka*cAMT/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*(T1[cTime]-dTime-TL)) + (k21-b)/(Ka-b)/(a-b)*exp(-b*(T1[cTime]-dTime-TL)) + (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*(T1[cTime]-dTime-TL)))
  }
  return(rowSums(Res))
}

r1 = nlr(fPK11pk, dPK11pk, pNames=c("Vc", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(10, 0.7, 0.1, 0.6, 0.07, 0.4), Error="P") ; r1
plot(e$Residual, type="o") ; abline(h=0)




## fitting PD model after fixed PK parameters
fPK11pd = function(THETA)
{
  EC50 = THETA[1]
  Emax = THETA[2]

  Vc  = 19.253893         
  Ka  = 2.048674
  k21 = 0.063026789
  a   = 0.3164563 # alpha
  b   = 0.053598374 # beta
  TL  = 0.31804026 # Tlag

  T1 = e$DATA[,"TIME"]
  Res = matrix(rep(0, e$nRec*nDose), nrow=e$nRec, ncol=nDose)
  for (i in 1:nDose) {
    dTime = DoseHist[i,"TIME"]
    cAMT = DoseHist[i,"AMT"]
    cTime = T1 > (dTime + TL)
    Res[cTime,i] = Ka*cAMT/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*(T1[cTime]-dTime-TL)) + (k21-b)/(Ka-b)/(a-b)*exp(-b*(T1[cTime]-dTime-TL)) + (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*(T1[cTime]-dTime-TL)))
  }
  Cp = rowSums(Res)
  
  E = Emax*Cp/(EC50 + Cp)
  return(E)
}
r2 = nlr(fPK11pd, dPK11pd, pNames=c("EC50", "Emax"), IE=c(2, 100)) ; r2
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK11pd, dPK11pd, pNames=c("EC50", "Emax"), IE=c(2, 100))
