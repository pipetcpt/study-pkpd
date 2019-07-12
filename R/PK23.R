setwd("D:/Rt/Gab/")
require(wnl)
dPK23 = read.csv("PK23.csv", skip=1)
colnames(dPK23) = c("Ca", "DV", "Qh") ; dPK23

fPK23a = function(THETA) # well stirred model
{
  Ca = e$DATA[,"Ca"]
  Qh = e$DATA[,"Qh"]
  CLint = THETA[1]
  Eh = CLint/(Qh + CLint)  # eq 23:1
  Cav = Ca * Eh
  return(Cav)
}

nlr(fPK23a, dPK23, pNames=c("CLint"), IE=c(5))


## Figure 23.1
plot(dPK23[,"Ca"], dPK23[,"DV"], xlim=c(10,45), ylim=c(5,30), xlab="Arterial concentration (ug/L)", ylab="Arterio-venous difference (ug/L)", pch=19)
#points(dPK23[,"Ca"],  fPK23a(e$PE), pch=19)
lines(sort(dPK23[,"Ca"]), fPK23a(e$PE)[order(dPK23[,"Ca"])], type="s")
##


fPK23b = function(THETA) # parallel tube model
{
  Ca = e$DATA[,"Ca"]
  Qh = e$DATA[,"Qh"]
  CLint = THETA[1]
  Cav = Ca * (1 - exp(-CLint/Qh))  # eq 23:2
  return(Cav)
}

nlr(fPK23b, dPK23, pNames=c("CLint"), IE=c(5))

fPK23c = function(THETA) # distributed model
{
  Ca = e$DATA[,"Ca"]
  Qh = e$DATA[,"Qh"]
  CLint = THETA[1]
  EPS   = THETA[2]
  A     = CLint/Qh
  Eh = 1 - exp(-(A + 0.5*EPS*EPS*A*A))  # eq 23:3
  Cav = Ca * Eh
  return(Cav)
}

nlr(fPK23c, dPK23, pNames=c("CLint", "EPS"), IE=c(2, 0.5))

fPK23d = function(THETA) # dispersion model
{
  Ca = e$DATA[,"Ca"]
  Qh = e$DATA[,"Qh"]
  CLint = THETA[1]
  Dn    = THETA[2]
  Rn    = CLint/Qh
  a  = (1 + 4*Rn*Dn)^0.5
  BA = 4*a/((1 + a)^2*exp((a - 1)/2/Dn) - (1 - a)^2*exp(-(a + 1)/2/Dn)) # eq 23:6  
  Cav = Ca * (1 - BA)
  return(Cav)
}

nlr(fPK23d, dPK23, pNames=c("CLint", "Dn"), IE=c(5, 10000)) # note CV of Dn

