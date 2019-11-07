setwd("D:/Rt/Gab/")
require(wnl)
dPK51 = read.csv("PK51.csv", skip=1)
colnames(dPK51) = c("TIME", "DV", "Group", "GroupDesc")
dPK51 = dPK51[!is.na(dPK51[,"DV"]),] ; dPK51

ivAMT = 5000
TAU = 5
infRate = ivAMT/TAU # 1000
poAMT = 8000

require(deSolve)
PKde1 = function(t, y, p)
{
  if (t < TAU) {
    RateIn = infRate
  } else {
    RateIn = 0
  }
  dy1dt = (RateIn - p["CL"]*y[1] - p["CLd"]*y[1] + p["CLd"]*y[2])/p["Vc"]
  dy2dt = (p["CLd"]*y[1] - p["CLd"]*y[2])/p["Vt"]
  dy3dt = (p["CL"]*y[1] - p["CL2"]*y[3] - p["CLd2"]*y[3] + p["CLd2"]*y[4])/p["Vc2"]
  dy4dt = (p["CLd2"]*y[3] - p["CLd2"]*y[4])/p["Vt2"]
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt)))
}

PKde2 = function(t, y, p)
{
  if (t < p["Tlag"]) {
    Ka = 0
  } else {
    Ka = p["Ka"]
  }
  dy1dt = -Ka*y[1]
  dy2dt = (p["Fa"]*Ka*y[1] - p["CL"]*y[2] - p["CLd"]*y[2] + p["CLd"]*y[3])/p["Vc"]
  dy3dt = (p["CLd"]*y[2] - p["CLd"]*y[3])/p["Vt"]
  dy4dt = ((1 - p["Fa"])*Ka*y[1] +  p["CL"]*y[2] - p["CL2"]*y[4] - p["CLd2"]*y[4] + p["CLd2"]*y[5])/p["Vc2"]
  dy5dt = (p["CLd2"]*y[4] - p["CLd2"]*y[5])/p["Vt2"]
  return(list(c(dy1dt, dy2dt, dy3dt, dy4dt, dy5dt)))
}

T1 = sort(unique(c(0, dPK51[dPK51$Group <= 2, "TIME"], 1440)))
iTime1 = T1 %in% dPK51[dPK51$Group == 1,"TIME"]
iTime2 = T1 %in% dPK51[dPK51$Group == 2,"TIME"]

y1 = lsoda(y=c(0,0,0,0), times=T1, func=PKde1, parms=c(Vc=19.4177, Vt=10.287, Vc2=4.935, Vt2=57.3215, CL=0.564672, CLd=0.0698416, CL2=0.0806763, CLd2=0.550319))
y1
y1[iTime1, "1"]
y1[iTime2, "3"]
dPK51[dPK51$Group <= 2,]

T2 = sort(unique(c(0, dPK51[dPK51$Group >= 3, "TIME"], 1440)))
iTime3 = T2 %in% dPK51[dPK51$Group == 3,"TIME"]
iTime4 = T2 %in% dPK51[dPK51$Group == 4,"TIME"]

y2 = lsoda(y=c(poAMT,0,0,0,0), times=T2, func=PKde2, parms=c(Ka=0.032, Vc=19.4177, Vt=10.287, Vc2=4.935, Vt2=57.3215, CL=0.564672, CLd=0.0698416, CL2=0.0806763, CLd2=0.550319, Fa=0.244834, Tlag=20.718))
y2
y2[iTime3, "2"]
y2[iTime4, "4"]
dPK51[dPK51$Group >= 3, ]

# Figure 51.1, p 712
plot(0, 0.01, type="n", xlim=c(0, 1500), ylim=c(0.01, 1000), log="y", xlab="Time (min)", ylab="Concentration")
points(dPK51[dPK51$Group == 1,"TIME"], dPK51[dPK51$Group == 1,"DV"], pch=19) # iv parent
points(dPK51[dPK51$Group == 2,"TIME"], dPK51[dPK51$Group == 2,"DV"], pch=15) # iv metabolite
points(dPK51[dPK51$Group == 3,"TIME"], dPK51[dPK51$Group == 3,"DV"], pch=1)  # po parent
points(dPK51[dPK51$Group == 4,"TIME"], dPK51[dPK51$Group == 4,"DV"], pch=0)  # po metabolite
lines(y1[,"time"], y1[,"1"]) # iv parent
lines(y1[,"time"], y1[,"3"]) # iv metabolite
lines(y2[,"time"], y2[,"2"]) # po parent
lines(y2[,"time"], y2[,"4"]) # po metabolite
###

T1 = sort(unique(c(0, dPK51[dPK51$Group <= 2, "TIME"])))
iTime1 = T1 %in% dPK51[dPK51$Group == 1,"TIME"]
iTime2 = T1 %in% dPK51[dPK51$Group == 2,"TIME"]
T2 = sort(unique(c(0, dPK51[dPK51$Group >= 3, "TIME"])))
iTime3 = T2 %in% dPK51[dPK51$Group == 3,"TIME"]
iTime4 = T2 %in% dPK51[dPK51$Group == 4,"TIME"]

fPK51 = function(THETA)
{
  Ka   = THETA[1]
  Vc   = THETA[2]
  Vt   = THETA[3]
  Vc2  = THETA[4]
  Vt2  = THETA[5]
  CL   = THETA[6]
  CLd  = THETA[7]
  CL2  = THETA[8]
  CLd2 = THETA[9]
  Fa   = THETA[10]
  Tlag = THETA[11]
  y1 = lsoda(y=c(0,0,0,0), times=T1, func=PKde1, parms=c(Vc=Vc, Vt=Vt, Vc2=Vc2, Vt2=Vt2, CL=CL, CLd=CLd, CL2=CL2, CLd2=CLd2))
  y2 = lsoda(y=c(poAMT,0,0,0,0), times=T2, func=PKde2, parms=c(Ka=Ka, Vc=Vc, Vt=Vt, Vc2=Vc2, Vt2=Vt2, CL=CL, CLd=CLd, CL2=CL2, CLd2=CLd2, Fa=Fa, Tlag=Tlag))
  return(c(y1[iTime1, "1"], y1[iTime2, "3"], y2[iTime3, "2"], y2[iTime4, "4"]))
}
fPK51(c(0.032, 19.4177, 10.287, 4.935, 57.3215, 0.564672, 0.0698416, 0.0806763, 0.550319, 0.244834, 20.718))

nlr(fPK51, dPK51, pNames=c("Ka", "Vc", "Vt", "Vc2", "Vt2", "CL", "CLd", "CL2", "CLd2", "Fa", "Tlag"),
    IE=c(0.04, 20, 11, 5, 60, 1, 0.1, 0.1, 1, 0.5, 20.72), 
    LB=c(0.01, 10, 5, 3, 30, 0.5, 0.05, 0.05, 0.5, 0.1, 20),
    UB=c(0.1, 30, 20, 7, 90, 2, 0.2, 0.2, 2, 1, 21), Error="P") # COV failure, see NONMEM
e$PE # similar

# remove one awkward record (group==3 & time==21, dv==1)
dPK51b = dPK51[!(dPK51$Group==4 & dPK51$TIME==21),]

iTime4 = T2 %in% dPK51b[dPK51b$Group == 4,"TIME"]

nlr(fPK51, dPK51b, pNames=c("Ka", "Vc", "Vt", "Vc2", "Vt2", "CL", "CLd", "CL2", "CLd2", "Fa", "Tlag"),
    IE=c(0.04, 20, 11, 5, 60, 1, 0.1, 0.1, 1, 0.5, 20.8), 
    LB=c(0.01, 10, 5, 3, 30, 0.5, 0.05, 0.05, 0.5, 0.1, 10),
    UB=c(0.1, 30, 20, 7, 90, 2, 0.2, 0.2, 2, 1, 30), Error="P") # COV Failure
e$PE # similar

# remove one more awkward recrod
dPK51c = dPK51[dPK51$TIME != 21,] # Only Group 3 and 4 are affected.

T2 = sort(unique(c(0, dPK51b[dPK51c$Group >= 3, "TIME"])))
iTime3 = T2 %in% dPK51c[dPK51c$Group == 3,"TIME"]
iTime4 = T2 %in% dPK51c[dPK51c$Group == 4,"TIME"]

nlr(fPK51, dPK51c, pNames=c("Ka", "Vc", "Vt", "Vc2", "Vt2", "CL", "CLd", "CL2", "CLd2", "Fa", "Tlag"),
    IE=c(0.04, 20, 11, 5, 60, 1, 0.1, 0.1, 1, 0.5, 20.8),
    LB=c(0.01, 10, 5, 3, 30, 0.5, 0.05, 0.05, 0.5, 0.1, 10),
    UB=c(0.1, 30, 20, 7, 90, 2, 0.2, 0.2, 2, 1, 30), Error="P") # Success

