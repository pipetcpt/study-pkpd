require(wnl)
dPK03 = read.csv("data/PK03.csv", skip=1)
colnames(dPK03) = c("TIME", "DV") ; dPK03

DOSE = 20000
TIME = dPK03[,"TIME"]

fPK03a = function(THETA) # Prediction function
{
  Ka   = THETA[1]
  K    = THETA[2]
  tlag = THETA[3]
  V    = THETA[4]

  Cp  = DOSE/V*Ka/(Ka - K)*(exp(-K*(TIME - tlag)) - exp(-Ka*(TIME - tlag))) # eq 3:2
  return(Cp)
}

r1 = nlr(fPK03a, dPK03, pNames=c("Ka", "Ke", "tlag", "V/F"), IE=c(1, 0.5, 0.5, 300)) ; r1
Resid1 = e$Residual

fPK03b = function(THETA) # Prediction function
{
  Tabs = THETA[1]
  K    = THETA[2]
  V    = THETA[3]
  Rin  = DOSE/Tabs

  FT = TIME
  FT[FT >= Tabs] = Tabs
  Cp  = Rin/V/K * (1 - exp(-K*FT)) * exp(-K*(TIME - FT)) # eq 3:1
  return(Cp)
}

r2 = nlr(fPK03b, dPK03, pNames=c("Tabs", "ke", "V"), IE=c(1, 0.5, 300), SecNames=c("CL"), SecForms=c(~V*ke)) ; r2
Resid2 = e$Residual

# Figure 3.1 p 484
plot(dPK03[,"TIME"], dPK03[,"DV"], xlim=c(0, 10), ylim=c(0, 90), xlab="Time (h)", ylab="Concentration (ug/L)", pch=16)
lines(TIME, dPK03[,"DV"] - Resid1)
lines(TIME, dPK03[,"DV"] - Resid2, lty=2)
###

# or #
dev.new()
plot(dPK03[,"TIME"], dPK03[,"DV"], xlim=c(0, 10), ylim=c(0, 90), xlab="Time (h)", ylab="Concentration (ug/L)", pch=16)
lines(TIME, fPK03a(r1$Est["PE", 1:4]))
lines(TIME, fPK03b(r2$Est["PE", 1:3]), lty=2)
###
