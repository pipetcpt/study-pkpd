require(wnl)

dPK06 = read.csv("data/PK06.csv", skip=1)
colnames(dPK06) = c("TIME", "DV", "CMT") ; dPK06

Div = 12500
Dpo = 25000

# Numbers in Table 6.2 p502  # Error ? Au,po/Dpo / Au,iv/Div = 1.10 (110 %) Data error ???
# > 586*1204
# [1] 705544
# > 447.5*1128
# [1] 504780
# > 394*863
# [1] 340022
# > 132*1591
# [1] 210012


## Fitting IV and urine data only
dPK06a = dPK06[dPK06[,"CMT"] == 1 | dPK06[,"CMT"] == 3,] ; dPK06a # IV data only
fPK06a = function(THETA)
{
  V  = THETA[1]
  Cl = THETA[2]
  fe = THETA[3]

  T1 = e$DATA[e$DATA[,"CMT"]==1,"TIME"]
  Civ = Div/V*exp(-Cl/V*T1)       # Eq 6:1

  T2 = e$DATA[e$DATA[,"CMT"]==3,"TIME"]
  Xu = fe*Div*(1 - exp(-Cl/V*T2)) # Eq 6:2

  return(c(Civ, Xu))
}

nlr(fPK06a, dPK06a, pNames=c("V", "Cl", "fe"), IE=c(300, 5, 0.2))
wnl5(fPK06a, dPK06a, pNames=c("V", "Cl", "fe"), IE=c(300, 5, 0.2))

## Fitting Oral and urine data only
dPK06b = dPK06[dPK06[,"CMT"] == 2 | dPK06[,"CMT"] == 4,] ; dPK06b # oral data only
fPK06b = function(THETA)
{
  V    = THETA[1]
  Cl   = THETA[2]
  fe   = THETA[3]
  Fo   = THETA[4] # BA for oral
  Ka   = THETA[5]
  Tlag = THETA[6]
  ke   = Cl/V

  T1 = e$DATA[e$DATA[,"CMT"]==2,"TIME"]
  Cpo = Ka*Fo*Dpo/V/(Ka - ke)*(exp(-ke*(T1 - Tlag)) - exp(-Ka*(T1 - Tlag)))   # Eq 6:3
  Cpo[Cpo < 0] = 0

  T2 = e$DATA[e$DATA[,"CMT"]==4,"TIME"]
  Xu = fe*Ka*Fo*Dpo*(1/Ka + exp(-ke*(T2 - Tlag))/(ke - Ka) - ke*exp(-Ka*(T2 - Tlag))/Ka/(ke - Ka)) # Erratum in Eq 6:4, See Gibaldi 2e eq 1.122 p39 or WinNonlin model file

  return(c(Cpo, Xu))
}
#fPK06b(c(300, 5, 0.2, 0.9, 1, 0.4)) # fe & Fo are identifiable?
#fPK06b(c(300, 5, 0.9, 0.2, 1, 0.4))

# Difficult to fit
nlr(fPK06b, dPK06b, pNames=c("V", "Cl", "fe", "F", "Ka", "Tlag"), IE=c(300, 5, 0.2, 0.9, 1, 0.4), LB=c(30, 0.5, 0.02, 0.09, 0.1, 0.04), UB=c(3000, 50, 2, 9, 10, 1))
nlr(fPK06b, dPK06b, pNames=c("V", "Cl", "fe", "F", "Ka", "Tlag"), IE=c(300, 5, 0.2, 0.9, 1, 0.4), LB=c(200, 3, 0.01, 0.09, 0.1, 0), UB=c(400, 10, 0.3, 1, 5, 1))
wnl5(fPK06b, dPK06b, pNames=c("V", "Cl", "fe", "F", "Ka", "Tlag"), IE=c(300, 5, 0.2, 0.9, 1, 0.4), LB=c(200, 3, 0.01, 0.09, 0.1, 0), UB=c(400, 10, 0.3, 1, 5, 1))

## Fitting both IV and oral data
fPK06c = function(THETA)
{
  V    = THETA[1]
  Cl   = THETA[2]
  fe   = THETA[3]
  Fo   = THETA[4] # BA for oral
  Ka   = THETA[5]
  Tlag = THETA[6]
  ke   = Cl/V

  T1 = e$DATA[e$DATA[,"CMT"]==1,"TIME"]
  Civ = Div/V*exp(-Cl/V*T1)   # Eq 6:1

  T2 = e$DATA[e$DATA[,"CMT"]==2,"TIME"]
  Cpo = Ka*Fo*Dpo/V/(Ka - ke)*(exp(-ke*(T2 - Tlag)) - exp(-Ka*(T2 - Tlag)))   # Eq 6:3
  Cpo[Cpo < 0] = 0

  T3 = e$DATA[e$DATA[,"CMT"]==3,"TIME"]
  Xuiv = fe*Div*(1 - exp(-Cl/V*T3)) # Eq 6:2

  T4 = e$DATA[e$DATA[,"CMT"]==4,"TIME"]
  Xupo = fe*Ka*Fo*Dpo*(1/Ka + exp(-ke*(T4 - Tlag))/(ke - Ka) - ke*exp(-Ka*(T4 - Tlag))/Ka/(ke - Ka)) # Erratum in Eq 6:4, See Gibaldi 2e eq 1.122 p39 or WinNonlin model file

  return(c(Civ, Cpo, Xuiv, Xupo))
}

nlr(fPK06c, dPK06, pNames=c("V", "Cl", "fe", "F", "Ka", "Tlag"), IE=c(300, 5, 0.2, 0.9, 1, 0.4))
wnl5(fPK06c, dPK06, pNames=c("V", "Cl", "fe", "F", "Ka", "Tlag"), IE=c(300, 5, 0.2, 0.9, 1, 0.4)) # require wnl_0.3.0

