setwd("D:/Rt/Gab/")
require(wnl)
dPK16Conc = read.csv("PK16_data.csv", skip=1, as.is=TRUE)
colnames(dPK16Conc) = c("TIME", "DV", "CMT") ; dPK16Conc
dPK16Conc = dPK16Conc[dPK16Conc[,"DV"] != "Missing",]
dPK16Conc[,"DV"] = as.numeric(dPK16Conc[,"DV"]) ; dPK16Conc
dPK16Dose = read.csv("PK16_dose.csv") ; dPK16Dose

Dinf1 = 538
Tau1 = 0.983
Dinf2 = 3390
Tau2 = 23.95
Rate1 = Dinf1/Tau1  # 538/0.983, 547.3042
Rate2 = Dinf2/(Tau2 - Tau1)  # 3390/(23.95 - 0.983), 147.6031

require(deSolve)
PKde = function(t, y, p)
{
  if (t < Tau1) {
    RateIn = Rate1
  } else if (t < Tau2) {
    RateIn = Rate2
  } else {
    RateIn = 0
  }

  dy1dt = (RateIn - p["Clm"]*y[1] - p["Clr"]*y[1] - p["Cld"]*y[1] + p["Cld"]*y[2])/p["Vc"] # central, plasma
  dy2dt = (p["Cld"]*y[1] - p["Cld"]*y[2])/p["Vt"] # peripheral, tissue
  dy3dt = p["Clr"]*y[1]
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

Times = sort(c(0, dPK16Conc[, "TIME"]))
iTime = 2:length(Times)

tTime = dPK16Conc[,c("TIME","CMT")] ;
tTime = tTime[order(tTime$TIME),] ; tTime
iTime1 = tTime$CMT == 1
iTime2 = tTime$CMT == 2

y = lsoda(y=c(0, 0, 0), times=Times, func=PKde, parms=c(Vc=1.6, Vt=0.16, Clm=0.05, Clr=0.31, Cld=0.03)) ; y
y2 = y[iTime,]
y3 = y2[iTime1,2] ; y3
y4 = y2[iTime2,4] ; y4

plot(y[,1], y[,4], type="l", lty=2)
lines(y[,1], y[,2])
points(dPK16Conc[dPK16Conc$CMT == 2,"TIME"], dPK16Conc[dPK16Conc$CMT == 2,"DV"], pch=2)
points(dPK16Conc[dPK16Conc$CMT == 1,"TIME"], dPK16Conc[dPK16Conc$CMT == 1,"DV"])

fPK16 = function(THETA)
{
  Vc  = THETA[1]
  Vt  = THETA[2]
  Clm = THETA[3]
  Clr = THETA[4]
  Cld = THETA[5]

  y = lsoda(y=c(0, 0, 0), times=Times, func=PKde, parms=c(Vc=Vc, Vt=Vt, Clm=Clm, Clr=Clr, Cld=Cld))
  y2 = y[iTime,]

  return(c(y2[iTime1,"1"],y2[iTime1,"3"]))
}

nlr(fPK16, dPK16Conc, pNames=c("Vc", "Vt", "Clm", "Clr", "Cld"),
    IE=c(1.5, 0.2, 0.05, 0.32, 0.03),
    LB=c(0.5, 0.1, 0.01, 0.1, 0.01),
    UB=c(3, 0.4, 0.2, 0.5, 0.05)) # fitting failure -> Use NONMEM PK16.CTL

nlr(fPK16, dPK16Conc, pNames=c("Vc", "Vt", "Clm", "Clr", "Cld"),
    IE=c(1.5, 0.2, 0.05, 0.32, 0.03),
    LB=c(0.5, 0.1, 0.01, 0.1, 0.01),
    UB=c(3, 0.4, 0.2, 0.5, 0.05), Error="P") # fitting failure -> Use NONMEM PK16b.CTL

nlr(fPK16, dPK16Conc, pNames=c("Vc", "Vt", "Clm", "Clr", "Cld"),
    IE=c(1.5, 0.2, 0.05, 0.32, 0.03),
    LB=c(0.5, 0.1, 0.01, 0.1, 0.01),
    UB=c(3, 0.4, 0.2, 0.5, 0.05), Error="C") # fitting failure -> Use NONMEM

nlr(fPK16, dPK16Conc, pNames=c("Vc", "Vt", "Clm", "Clr", "Cld"),
    IE=c(1.5, 0.2, 0.05, 0.32, 0.03),
    LB=c(0.5, 0.1, 0.01, 0.1, 0.01),
    UB=c(3, 0.4, 0.2, 0.5, 0.05), Error="POIS") # fitting failure -> Use NONMEM PK16c.CTL

wnl5(fPK16, dPK16Conc, pNames=c("Vc", "Vt", "Clm", "Clr", "Cld"), IE=c(2, 0.2, 0.1, 0.5, 0.1))
wnl5(fPK16, dPK16Conc, pNames=c("Vc", "Vt", "Clm", "Clr", "Cld"), IE=c(2, 0.2, 0.1, 0.5, 0.1), Error="P")
wnl5(fPK16, dPK16Conc, pNames=c("Vc", "Vt", "Clm", "Clr", "Cld"), IE=c(2, 0.2, 0.1, 0.5, 0.1), Error="POIS")
