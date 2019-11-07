#setwd("D:/Rt/Gab/")

dPK10 = read.csv("data-old/PK10.csv", skip=1, as.is=TRUE)
colnames(dPK10) = c("TIME", "DV", "GRP", "ADM") ; dPK10

Div = 100
Dpo = 500

## iv only
fPK10a = function(THETA)
{
  Vc  = THETA[1]
  k21 = THETA[2]
  a   = THETA[3] # alpha
  b   = THETA[4] # beta
  
  TIME = e$DATA[e$DATA[,"ADM"]=="iv", "TIME"]
  Cp = Div/Vc*((k21-a)/(b-a)*exp(-a*TIME) + (k21-b)/(a-b)*exp(-b*TIME))
  return(Cp)
}

nlr(fPK10a, dPK10, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(60, 0.06, 0.06, 0.007))
plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK10a, dPK10, pNames=c("Vc", "k21", "alpha", "beta"), IE=c(60, 0.06, 0.06, 0.007))

## iv and oral without tlag
fPK10b = function(THETA)
{
  Vc  = THETA[1]
  Ka  = THETA[2]
  k21 = THETA[3] 
  a   = THETA[4] # alpha
  b   = THETA[5] # beta
  Fa  = THETA[6] # F, bioavailability
  
  T1 = e$DATA[e$DATA[,"ADM"]=="oral", "TIME"]
  Co = Ka*Fa*Dpo/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*T1) + (k21-b)/(Ka-b)/(a-b)*exp(-b*T1) + (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*T1)) # Erratum in eq 10:2

  T2 = e$DATA[e$DATA[,"ADM"]=="iv", "TIME"]
  Ci = Div/Vc*((k21-a)/(b-a)*exp(-a*T2) + (k21-b)/(a-b)*exp(-b*T2))

  return(c(Co, Ci))
}

nlr(fPK10b, dPK10, pNames=c("Vc", "Ka", "k21", "alpha", "beta", "F"), IE=c(60, 0.04, 0.026, 0.06, 0.007, 0.4))
plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK10b, dPK10, pNames=c("Vc", "Ka", "k21", "alpha", "beta", "F"), IE=c(60, 0.04, 0.026, 0.06, 0.007, 0.4))

## iv and oral with tlag
fPK10c = function(THETA)
{
  Vc  = THETA[1]
  Ka  = THETA[2]
  k21 = THETA[3] 
  a   = THETA[4] # alpha
  b   = THETA[5] # beta
  Fa  = THETA[6] # F, bioavailability
  TL  = THETA[7] # Tlag
  
  T1 = e$DATA[e$DATA[,"ADM"]=="oral", "TIME"]
  Co = Ka*Fa*Dpo/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*(T1-TL)) + (k21-b)/(Ka-b)/(a-b)*exp(-b*(T1-TL)) + (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*(T1-TL))) # Erratum in eq 10:2
  Co[Co < 0] = 0 # remove negative concentrationss before tlag

  T2 = e$DATA[e$DATA[,"ADM"]=="iv", "TIME"]
  Ci = Div/Vc*((k21-a)/(b-a)*exp(-a*T2) + (k21-b)/(a-b)*exp(-b*T2))

  return(c(Co, Ci))
}

nlr(fPK10c, dPK10, pNames=c("Vc", "Ka", "k21", "alpha", "beta", "F", "Tlag"), IE=c(60, 0.04, 0.026, 0.06, 0.007, 0.4, 3))
plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK10c, dPK10, pNames=c("Vc", "Ka", "k21", "alpha", "beta", "F", "Tlag"), IE=c(60, 0.04, 0.026, 0.06, 0.007, 0.4, 3))
