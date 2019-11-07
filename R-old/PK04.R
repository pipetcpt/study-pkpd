require(wnl)
dPK04 = read.csv("data/PK04.csv", skip=1)
colnames(dPK04) = c("TIME", "DV") ; dPK04

DHPK04 = read.csv("data/PK04-dose.csv")
colnames(DHPK04) = c("TIME", "AMT") ; DHPK04

fPK04a = function(THETA) # Eq 4:3
{
  V  = THETA[1]
  Ka = THETA[2]
  K  = THETA[3]
  
  TIME = e$DATA[,"TIME"]
  nTime = length(TIME)
  nDose = nrow(DHPK04)
  
  Res = matrix(rep(0, nTime*nDose), nrow=nTime, ncol=nDose)
  
  for (i in 1:nDose) {
    dTime = DHPK04[i,"TIME"]
    AMT = DHPK04[i,"AMT"]
    cTime = TIME > dTime
    Res[cTime,i] = Ka * AMT / V / (Ka - K) * (exp(-K*(TIME[cTime] - dTime)) - exp(-Ka*(TIME[cTime] - dTime)))
  }
  
  return(rowSums(Res))  
}

r1 = nlr(fPK04a, dPK04, pNames=c("V/F", "Ka", "K"), IE=c(70, 0.25, 0.1)) ; r1
r1 = nlr(fPK04a, dPK04, pNames=c("V/F", "Ka", "K"), IE=c(70, 0.25, 0.1), Error="POIS") ; r1
condNum1 = sqrt(max(e$EigenVal)/min(e$EigenVal)) ; condNum1
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK04a, dPK04, pNames=c("V/F", "Ka", "K"), IE=c(70, 0.25, 0.1), Error="POIS")

## 
fPK04b = function(THETA) # Eq 4:4
{
  V  = THETA[1]
  Ka = THETA[2]
  K  = THETA[3]
  tlag = THETA[4]
  
  TIME = e$DATA[,"TIME"]
  nTime = e$nRec
  nDose = nrow(DHPK04)
  
  Res = matrix(rep(0, nTime*nDose), nrow=nTime, ncol=nDose)
  for (i in 1:nDose) {
    dTime = DHPK04[i,"TIME"]
    cAMT = DHPK04[i,"AMT"]
    iTime = TIME > (dTime + tlag)
    Res[iTime,i] = Ka * cAMT / V / (Ka - K) * (exp(-K*(TIME[iTime] - dTime - tlag)) - exp(-Ka*(TIME[iTime] - dTime - tlag)))
  }
  
  return(rowSums(Res))  
}

#r2 = nlr(fPK04b, dPK04, pNames=c("V/F", "Ka", "K", "Tlag"), IE=c(70, 0.25, 0.1, 1)) ; r2 # flip-flop phenomenon
r2 = nlr(fPK04b, dPK04, pNames=c("V/F", "Ka", "K", "Tlag"), IE=c(70, 0.25, 0.1, 1), LB=c(50, 0, 0, 0), Error="POIS") ; r2
r2 = nlr(fPK04b, dPK04, pNames=c("V/F", "Ka", "K", "Tlag"), IE=c(70, 0.25, 0.1, 1), LB=c(50, 0, 0, 0), Error="C") ; r2
condNum2 = sqrt(max(e$EigenVal)/min(e$EigenVal)) ; condNum2
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
#wnl5(fPK04b, dPK04, pNames=c("V/F", "Ka", "K", "Tlag"), IE=c(70, 0.25, 0.1, 1), Error="POIS") # fitting failure
wnl5(fPK04b, dPK04, pNames=c("V/F", "Ka", "K", "Tlag"), IE=c(70, 0.25, 0.1, 1), LB=c(50, 0, 0, 0))

## 
fPK04c = function(THETA) # Eq 4:5
{
  V  = THETA[1]
  Kp = THETA[2]
  
  TIME = e$DATA[,"TIME"]
  nTime = length(TIME)
  nDose = nrow(DHPK04)
  
  Res = matrix(rep(0, nTime*nDose), nrow=nTime, ncol=nDose)
  
  for (i in 1:nDose) {
    dTime = DHPK04[i,"TIME"]
    AMT = DHPK04[i,"AMT"]
    cTime = TIME > dTime
    Res[cTime,i] = Kp * AMT / V * (TIME[cTime] - dTime) * exp(-Kp*(TIME[cTime] - dTime))
  }
  
  return(rowSums(Res))  
}

r3 = nlr(fPK04c, dPK04, pNames=c("V/F", "Kp"), IE=c(70, 0.25), Error="POIS") ; r3
condNum3 = sqrt(max(e$EigenVal)/min(e$EigenVal)) ; condNum3
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK04c, dPK04, pNames=c("V/F", "Kp"), IE=c(70, 0.25), Error="POIS")

## 
fPK04d = function(THETA) # Eq 4:6
{
  V  = THETA[1]
  Kp = THETA[2]
  tlag = THETA[3]
  
  TIME = e$DATA[,"TIME"]
  nTime = length(TIME)
  nDose = nrow(DHPK04)
  
  Res = matrix(rep(0, nTime*nDose), nrow=nTime, ncol=nDose)
  
  for (i in 1:nDose) {
    dTime = DHPK04[i,"TIME"]
    AMT = DHPK04[i,"AMT"]
    cTime = TIME > (dTime + tlag)
    Res[cTime,i] = Kp * AMT / V * (TIME[cTime] - dTime - tlag) * exp(-Kp*(TIME[cTime] - dTime - tlag))
  }
  
  return(rowSums(Res))  
}

r4 = nlr(fPK04d, dPK04, pNames=c("V/F", "Kp", "Tlag"), IE=c(70, 0.25, 1), Error="POIS") ; r4 # fitting failure
r4 = nlr(fPK04d, dPK04, pNames=c("V/F", "Kp", "Tlag"), IE=c(70, 0.25, 1)) ; r4
condNum4 = sqrt(max(e$EigenVal)/min(e$EigenVal)) ; condNum4
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK04d, dPK04, pNames=c("V/F", "Kp", "Tlag"), IE=c(70, 0.25, 1))
