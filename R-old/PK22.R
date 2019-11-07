setwd("D:/Rt/Gab/")
require(wnl)
dPK22 = read.csv("PK22.csv", skip=1)
colnames(dPK22) = c("TIME", "DV") ; dPK22

dPK22dose = read.csv("PK22_dose.csv") # dPK22dose
dPK22dose[,"DUR"] = dPK22dose[,"AMT"]/dPK22dose[,"RATE"] # dPK22dose
dPK22dose[,"TIME2"] = dPK22dose["TIME"] + dPK22dose[,"DUR"] ; dPK22dose

DoseHist = rbind(dPK22dose[,c("TIME", "RATE")], cbind(TIME=dPK22dose[,"TIME2"],RATE=0)) # DoseHist
DoseHist = DoseHist[order(DoseHist$TIME),] ; DoseHist

#sTIME = seq(0, 73, by=0.1)
#iTIME = findInterval(sTIME, DoseHist[,"TIME"])
#cbind(sTIME, iTIME, DoseHist[iTIME,])

PKde = function(t, y, p)
{
  RateIn = DoseHist[findInterval(t, DoseHist[,"TIME"]),"RATE"]
  CL = p["CLs"]*(1 + y[3])                                           # eq 22:3
  dy1dt = (RateIn - CL*y[1] - p["CLd"]*y[1] + p["CLd"]*y[2])/p["Vc"] # eq 22:1
  dy2dt = (p["CLd"]*y[1] - p["CLd"]*y[2])/p["Vt"]                    # eq 22:2
  dy3dt = p["Kout"]*(p["E0"] + y[1]) - p["Kout"]*y[3]                # eq 22:4
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

# Figure 22.1, p 581
plot(0, 0, type="n", xlim=c(0, 100), ylim=c(0, 700), xlab="Time (h)", ylab="Concentration (ug/L)")
gTIME = seq(0, 100, by=0.1)
y = lsoda(y=c(0, 0, 132.864), times=gTIME, func=PKde, parms=c(Vc=150, CLs=0.04, CLd=97.8, Vt=54, Kout=0.024, E0=132.864))

points(dPK22[,"TIME"], dPK22[,"DV"])
lines(y[,"time"], y[,"1"])
##

Times = c(0, dPK22[, "TIME"])
iTime = 2:length(Times)

fPK22 = function(THETA)
{
  Vc   = THETA[1]
  CLs  = THETA[2]
  CLd  = THETA[3]
  Vt   = THETA[4]
  Kout = THETA[5]
  E0   = THETA[6]

  y = lsoda(y=c(0, 0, 100), times=Times, func=PKde, parms=c(Vc=Vc, CLs=CLs, CLd=CLd, Vt=Vt, Kout=Kout, E0=E0))
  return(y[iTime, "1"])
}
fPK22(c(150, 0.04, 97.8, 54, 0.024, 132.864))
fPK22(c(155, 0.05, 120, 60, 0.03, 100))

nlr(fPK22, dPK22, pNames=c("Vc", "CLs", "CLd", "Vt", "Kout", "E0"),
    IE=c(155, 0.05, 120, 60, 0.03, 100),
    LB=c(100, 0.01, 60, 30, 0.01, 50),
    UB=c(200, 0.09, 200, 120, 0.05, 200)) # fitting failure

nlr(fPK22, dPK22, pNames=c("Vc", "CLs", "CLd", "Vt", "Kout", "E0"),
    IE=c(155, 0.05, 120, 60, 0.03, 100),
    LB=c(100, 0.01, 60, 30, 0.01, 50),
    UB=c(200, 0.09, 200, 120, 0.05, 200), Error="P") # fitting failure -> Use NONMEM



## microconstant model


PKde2 = function(t, y, p)
{
  RateIn = DoseHist[findInterval(t, DoseHist[,"TIME"]),"RATE"]
  CL = p["CLs"]*(1 + y[3])                                           # eq 22:3
  dy1dt = (RateIn - CL*y[1])/p["Vc"] - p["k12"]*y[1] + p["k21"]*y[2]
  dy2dt = p["k12"]*y[1] - p["k21"]*y[2]
  dy3dt = p["Kout"]*(p["E0"] + y[1]) - p["Kout"]*y[3]                # eq 22:4
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

# Figure 22.1, p 581
plot(0, 0, type="n", xlim=c(0, 100), ylim=c(0, 700), xlab="Time (h)", ylab="Concentration (ug/L)")
gTIME = seq(0, 100, by=0.1)
y = lsoda(y=c(0, 0, 137), times=gTIME, func=PKde2, parms=c(Vc=146, CLs=0.04, k12=0.8, k21=1.98, Kout=0.023, E0=137))

points(dPK22[,"TIME"], dPK22[,"DV"])
lines(y[,"time"], y[,"1"])
##

Times = c(0, dPK22[, "TIME"])
iTime = 2:length(Times)

fPK22b = function(THETA)
{
  Vc   = THETA[1]
  CLs  = THETA[2]
  k12  = THETA[3]
  k21   = THETA[4]
  Kout = THETA[5]
  E0   = THETA[6]

  y = lsoda(y=c(0, 0, 100), times=Times, func=PKde2, parms=c(Vc=Vc, CLs=CLs, k12=k12, k21=k21, Kout=Kout, E0=E0))
  return(y[iTime, "1"])
}
fPK22b(c(150, 0.04, 0.8, 1.98, 0.023, 132.864))
fPK22b(c(155, 0.05, 0.5, 1.25, 0.03, 100))

nlr(fPK22, dPK22, pNames=c("Vc", "CLs", "CLd", "Vt", "Kout", "E0"),
    IE=c(155, 0.05, 0.5, 1.25, 0.03, 100),
    LB=c(100, 0.01, 0.1, 0.5, 0.01, 50),
    UB=c(200, 0.09, 1.2, 3, 0.05, 200)) # fitting failure

nlr(fPK22, dPK22, pNames=c("Vc", "CLs", "CLd", "Vt", "Kout", "E0"),
    IE=c(155, 0.05, 0.5, 1.25, 0.03, 100),
    LB=c(100, 0.01, 0.1, 0.5, 0.01, 50),
    UB=c(200, 0.09, 1.2, 3, 0.05, 200), Error="P") # fitting failure -> Use NONMEM
