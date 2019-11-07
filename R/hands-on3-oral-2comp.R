

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
  Co = Ka*Dpo/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*T1) + 
                    (k21-b)/(Ka-b)/(a-b)*exp(-b*T1) + 
                    (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*T1)) # Erratum in eq 14:1
  
  return(Co)
}

r1 <- nlr(fPK14a, dPK14, 
          pNames=c("Vc/F", "Ka", "k21", "alpha", "beta"), 
          IE=c(350, 11, 1, 0.1, 0.01))

r1$Est

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
  Co = Ka*Dpo/Vc*((k21-a)/(Ka-a)/(b-a)*exp(-a*(T1-TL)) + 
                    (k21-b)/(Ka-b)/(a-b)*exp(-b*(T1-TL)) + 
                    (k21-Ka)/(a-Ka)/(b-Ka)*exp(-Ka*(T1-TL))) 
  Co[Co < 0] = 0 # remove negative concentrationss before tlag
  
  return(Co)
}

r2 <- nlr(fPK14b, dPK14, 
          pNames=c("Vc/F", "Ka", "k21", "alpha", "beta", "Tlag"), 
          IE=c(150, 11, 0.12, 0.1, 0.01, 0.05))
r2$Est

# Figure 2.3, p 480

plot(dPK14[,"TIME"], dPK14[,"DV"], xlim=c(0, 25), ylim=c(0, 250), xlab="Time (min)", ylab="Concentration (ug/L)", pch=16)
TIME = dPK14[,"TIME"]
lines(TIME, fPK14a(r1$Est["PE", 1:5]), lty=2)
lines(TIME, fPK14b(r2$Est["PE", 1:6]))
###
