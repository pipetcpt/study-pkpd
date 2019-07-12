require(wnl)

## PK01
dPK01 = read.csv("D:/Rt/Gab/PK01.csv")
colnames(dPK01) = c("ID", "TIME", "DV")

fPK01 = function(THETA) # Prediction function
{
  DOSE = Dose
  TIME = e$DATA[,"TIME"]
  V  = THETA[1]
  K  = THETA[2]

  P  = DOSE/V * exp(-K*TIME)
  return(P)
}

Dose = 100000
IDs = unique(dPK01[,"ID"])
nID = length(IDs)
for (i in 1:nID) {
  Data = dPK01[dPK01$ID == IDs[i],]
  Res = nlr(fPK01, Data, pNames=c("V", "k"), IE=c(20, 0.2))
#  Res = nlr(fPK01, Data, pNames=c("V", "k"), IE=c(20, 0.2), Error="P")
#  Res = nlr(fPK01, Data, pNames=c("V", "k"), IE=c(20, 0.2), Error="C")
  print(Res)
}

## PK02
dPK02 = read.csv("D:/Rt/Gab/PK02.csv", skip=1)
colnames(dPK02) = c("TIME", "DV")

fPK02 = function(THETA) # Prediction function
{
  DOSE = Dose
  TIME = e$DATA[,"TIME"]

  BA   = BA
  Ka   = THETA[1]
  V    = THETA[2]
  K    = THETA[3]
  tlag = THETA[4]

  P  = BA*DOSE/V*Ka/(Ka - K) * (exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag)))
  return(P)
}

Dose = 100
#BA = 0.6100976 # from NCA   AUCinf.po / AUCinf.iv
BA = 1
nlr(fPK02, dPK02, pNames=c("ka", "V", "k", "tlag"), IE=c(0.1, 30, 0.05, 20))
wnl5(fPK02, dPK02, pNames=c("ka", "V", "k", "tlag"), IE=c(0.1, 30, 0.05, 20))
