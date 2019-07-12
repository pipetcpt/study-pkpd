# PD 44
require(wnl)
setwd("D:/Rt/PD")

dPD44 = read.csv("PD44.csv")
colnames(dPD44) = c("Cp", "DV", "ID", "GROUP")
dPD44 # Data is different with the Book Table 44.4
dim(dPD44)

Cp1 = dPD44[dPD44$ID == 1, "Cp"]
Cp2 = dPD44[dPD44$ID == 2, "Cp"]
nCp2 = length(Cp2)
fPD44 = function(THETA)
{
  E0   = THETA[1]
  Shft = THETA[2]
  a    = THETA[3]
  b    = THETA[4]
  E1   = E0 + (Cp1/a)^b

  Cp2a = Cp2[Cp2 > Shft]
  Cp2b = Cp2[Cp2 > Shft]
  E2   = rep(E0, nCp2)
  E2[Cp2 > a*Shft] = E0 + (Cp2[Cp2 > a*Shft]/a - Shft)^b
  return(c(E1, E2))  
}  
y0 = fPD44(c(10.7004, 1.51, 0.25646, 1.55399)) ; y0
length(y0)

r1 = nlr(fPD44, dPD44, c("E0", "Shft", "a", "b"), c(10, 0.5, 0.1, 1), 
         SecNames=c("Ce20", "Cen20"), 
         SecForms=c(~(20 - E0)^(1/b)*a, ~((20 - E0)^(1/b) + Shft)*a)) ; r1 # Different

# Figure 44.1 : note the difference in "Normal" curve
dev.new(width=7, height=5)
plot(dPD44[dPD44$ID==1,"Cp"], dPD44[dPD44$ID==1,"DV"], pch=16, xlim=c(0, 3), ylim=c(0, 60), xlab="Concentration", ylab="Reponse (PT)", col="red")
lines(dPD44[dPD44$ID==1,"Cp"], y0[dPD44$ID==1], col="red")
points(dPD44[dPD44$ID==2,"Cp"], dPD44[dPD44$ID==2,"DV"], pch=16)
lines(dPD44[dPD44$ID==2,"Cp"], y0[dPD44$ID==2])
abline(h=1:6*10, lty=3)

# Figure 44.1 again with r1
dev.new(width=7, height=5)
plot(dPD44[dPD44$ID==1,"Cp"], dPD44[dPD44$ID==1,"DV"], pch=16, xlim=c(0, 3), ylim=c(0, 60), xlab="Concentration", ylab="Reponse (PT)", col="red")
lines(dPD44[dPD44$ID==1,"Cp"], r1$Prediction[dPD44$ID==1], col="red")
points(dPD44[dPD44$ID==2,"Cp"], dPD44[dPD44$ID==2,"DV"], pch=16)
lines(dPD44[dPD44$ID==2,"Cp"], r1$Prediction[dPD44$ID==2])
abline(h=1:6*10, lty=3)

