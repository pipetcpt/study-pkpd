#setwd("D:/Rt/Gab/")
require(wnl)
dPK14 = read.csv("data-old/PK14.csv", skip=1)
colnames(dPK14) = c("TIME", "DV") ; dPK14

Dpo = 23158

## without lag
fPK14a = function(THETA)
{
  Vc  = THETA[1]
  Ka  = THETA[2]
  k21 = THETA[3]
  a   = THETA[4] # alpha
  b   = THETA[5] # beta

  T1 = e$DATA[,"TIME"]
  Co = Ka*Dpo/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*T1) + (k21-b)/(Ka-b)/(a-b)*exp(-b*T1) + (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*T1)) # Erratum in eq 14:1

  return(Co)
}

nlr(fPK14a, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta"), IE=c(350, 11, 1, 0.1, 0.01))
wnl5(fPK14a, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta"), IE=c(350, 11, 1, 0.1, 0.01))

nlr(fPK14a, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta"), IE=c(350, 11, 1, 0.1, 0.01), Error="P")
wnl5(fPK14a, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta"),
     IE=c(100, 10, 1, 0.1, 0.01),
     LB=c(50, 1, 0.1, 0.01, 0.001),
     UB=c(500, 20, 2, 1, 0.1), Error="P") # fitting failure

nlr(fPK14a, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta"), IE=c(350, 11, 1, 0.1, 0.01), Error="C") # COV step failure
e$PE

## with lag
fPK14b = function(THETA)
{
  Vc  = THETA[1]
  Ka  = THETA[2]
  k21 = THETA[3]
  a   = THETA[4] # alpha
  b   = THETA[5] # beta
  TL  = THETA[6] # Tlag

  T1 = e$DATA[,"TIME"]
  Co = Ka*Dpo/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*(T1-TL)) + (k21-b)/(Ka-b)/(a-b)*exp(-b*(T1-TL)) + (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*(T1-TL))) # Erratum in eq 14:2
  Co[Co < 0] = 0 # remove negative concentrationss before tlag

  return(Co)
}
# fPK14b(c(150, 5, 0.12, 0.1, 0.01, 0.05))

nlr(fPK14b, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(150, 11, 0.12, 0.1, 0.01, 0.05))
wnl5(fPK14b, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(150, 11, 0.12, 0.1, 0.01, 0.05))

nlr(fPK14b, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(150, 11, 0.12, 0.1, 0.01, 0.05), Error="P") # fitting failure
wnl5(fPK14b, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(150, 11, 1, 0.1, 0.01, 0.05), Error="P")

nlr(fPK14b, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(150, 11, 0.12, 0.1, 0.01, 0.05), Error="C") # COV step failure
e$PE

nlr(fPK14b, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(350, 11, 1, 0.1, 0.01, 0.05), Error="POIS") # fitting failure
wnl5(fPK14b, dPK14, pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), IE=c(150, 11, 1, 0.1, 0.01, 0.05), Error="POIS")


## Matrix exponential approach
require(Matrix)
#M = matrix(c(1,0,0,1),nrow=2)
#expm(M)

fPK14c = function(THETA)
{
  Vc  = THETA[1]
  Ka  = THETA[2]
  k21 = THETA[3]
  K   = THETA[4] # ke, k10
  k12 = THETA[5]
  TL  = THETA[6] # Tlag

  M = matrix(c(-Ka, 0, 0, Ka, -K-k12, k21, 0, k12, -k21), nrow=3, byrow=TRUE)
  TIME = e$DATA[,"TIME"]

  A0 = c(Dpo, 0, 0)
  nTime = length(TIME)
  Res = matrix(nrow=nTime, ncol=3)
  for (i in 1:nTime) {
    cT = TIME[i]
    Res[i,] = as.vector(expm(M*(cT - TL)) %*% A0)
  }
  Cp = Res[,2]/Vc
  Cp[Cp < 0] = 0
  return(Cp)
}
# fPK14c(c(83, 10, 0.66, 0.13, 0.10, 0.078))

nlr(fPK14c, dPK14, pNames=c("Vc/F", "Ka", "k21", "K", "k12", "TL"),
    IE=c(100, 10, 1, 0.2, 0.1, 0.1),
    LB=c(50, 5, 0.1, 0.1, 0.01, 0.01),
    UB=c(200, 100, 2, 1, 1, 1))

wnl5(fPK14c, dPK14, pNames=c("Vc/F", "Ka", "k21", "K", "k12", "TL"),
    IE=c(100, 10, 1, 0.2, 0.1, 0.1),
    LB=c(50, 5, 0.1, 0.1, 0.01, 0.01),
    UB=c(200, 100, 2, 1, 1, 1))

nlr(fPK14c, dPK14, pNames=c("Vc/F", "Ka", "k21", "K", "k12", "TL"),
    IE=c(90, 10, 0.2, 0.7, 0.13, 0.08),
    LB=c(80, 5, 0.1, 0.1, 0.01, 0.01),
    UB=c(100, 20, 1, 1, 1, 1), Error="P") # COV failure
e$PE

wnl5(fPK14c, dPK14, pNames=c("Vc/F", "Ka", "k21", "K", "k12", "TL"),
    IE=c(90, 10, 0.2, 0.7, 0.13, 0.08),
    LB=c(80, 5, 0.1, 0.1, 0.01, 0.01),
    UB=c(100, 20, 1, 1, 1, 1), Error="P")

nlr(fPK14c, dPK14, pNames=c("Vc/F", "Ka", "k21", "K", "k12", "TL"),
    IE=c(90, 10, 0.2, 0.7, 0.13, 0.08),
    LB=c(80, 5, 0.1, 0.1, 0.01, 0.01),
    UB=c(100, 20, 1, 1, 1, 1), Error="C") # COV failure
e$PE
