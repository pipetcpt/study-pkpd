# PD 3
require(wnl)
setwd("D:/Rt/PD")

dPD3a = read.csv("PD3n0800.csv") # PD3no800.csv
colnames(dPD3a) = c("Conc", "DV")

dPD3b = read.csv("PD3with800.csv")
colnames(dPD3b) = c("Conc", "DV")

fPD3a = function(THETA)
{
  E0   = THETA[1]
  Imax = THETA[2]
  IC50 = THETA[3]
  Conc = e$DATA[,"Conc"]
  Resp = E0 - Imax*Conc/(IC50 + Conc)
  return(Resp)
}

fPD3b = function(THETA)
{
  E0   = THETA[1]
  Imax = THETA[2]
  IC50 = THETA[3]
  n    = THETA[4]
  Conc = e$DATA[,"Conc"]
  Resp = E0 - Imax*Conc^n/(IC50^n + Conc^n)
  return(Resp)
}

r1da = nlr(fPD3a, dPD3a, c("E0", "Imax", "IC50"), c(175, 35, 120)) ; r1da
r1daw = wnl5(fPD3a, dPD3a, c("E0", "Imax", "IC50"), c(175, 35, 120)) ; r1daw

r1db = nlr(fPD3a, dPD3b, c("E0", "Imax", "IC50"), c(175, 35, 120)) ; r1db
r1dbw = wnl5(fPD3a, dPD3b, c("E0", "Imax", "IC50"), c(175, 35, 120)) ; r1dbw

r2da = nlr(fPD3b, dPD3a, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; r2da
r2daw = wnl5(fPD3b, dPD3a, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; r2daw

r2db = nlr(fPD3b, dPD3b, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; r2db
r2dbw = wnl5(fPD3b, dPD3b, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; r2dbw

# p735 eq 3.3
r1dbw
r2dbw
nRec = dim(dPD3b)[1]
WRSS1 = r1dbw$WRSS 
WRSS2 = r2dbw$WRSS
df1 = nRec - length(r1dbw$PE)
df2 = nRec - length(r2dbw$PE)

F.val = abs(WRSS1 - WRSS2)/abs(df1 - df2)/(WRSS2/df2) ; F.val
qf(1 - 0.05, df1 - df2, df2)
1 - pf(F.val, df1 - df2, df2)

# Figure 3.1
plot(dPD3b$Conc, dPD3b$DV, ylim=c(135, 180), xlab="Concentration (ug/L)", ylab="Response (mmHg)", pch=16)
abline(h=seq(140, 180, 5), lty=3)
x = 0:800
e$DATA = cbind(DV=0, Conc=x)
y1 = fPD3a(r1db$Est["PE", c("E0", "Imax", "IC50")])
y2 = fPD3b(r2db$Est["PE", c("E0", "Imax", "IC50", "n")])
lines(x, y1, col="red")
lines(x, y2)

# Figure 3.3
e$DATA = dPD3b
resid1 = dPD3b[,"DV"] - fPD3a(r1db$Est["PE", c("E0", "Imax", "IC50")])
resid2 = dPD3b[,"DV"] - fPD3b(r2db$Est["PE", c("E0", "Imax", "IC50", "n")])

dev.new()
plot(dPD3b$Conc, resid1, col="red", type="l", xlab="Concentration (ug/L)", ylab="Residual")
abline(h=seq(-3,4), lty=3)
abline(h=0)
lines(dPD3b$Conc, resid2, col="blue")

# Figure 3.4
dev.new()
plot(x, (y1 - y2)^2, ylim=c(0, 10), col="red", type="l", xlab="Concentration (ug/L)", ylab="Delta function")
abline(h=seq(2, 10, 2), lty=3)

# Figure 3.5
Fx   = expression(E0 - Imax*Conc^n/(IC50^n + Conc^n))
Conc = dPD3b[,"Conc"]
E0   = r2db$Est[1, "E0"]
Imax = r2db$Est[1, "Imax"]
IC50 = r2db$Est[1, "IC50"]
n    = r2db$Est[1, "n"]

Fd1 = deriv(Fx, "n", function.arg=c("E0", "Imax", "IC50", "n", "Conc"), func=TRUE)
g1 = attr(Fd1(E0, Imax, IC50, n, Conc), "gradient")

Fd2 = deriv(Fx, "IC50", function.arg=c("E0", "Imax", "IC50", "n", "Conc"), func=TRUE)
g2 = attr(Fd2(E0, Imax, IC50, n, Conc), "gradient")

# Left
dev.new()
plot(dPD3b[,"Conc"], dPD3b[,"DV"], xlab="Concentration (ug/L)", ylab="Response & derivative", ylim=c(0,250), pch=16)
lines(x, y1)
lines(Conc, g1*22 + E0)

# Right
dev.new()
plot(dPD3b[,"Conc"], dPD3b[,"DV"], xlab="Concentration (ug/L)", ylab="Response & derivative", ylim=c(0,250), pch=16)
lines(x, y1)
lines(Conc, g2*1000)

# Figure 3.6
Conc = dPD3b[,"Conc"]
nRec = length(Conc)
Pred = 171.4 - 34.66*Conc^2.031/(139.6^2.031 + Conc^2.031) ; Pred

s1 = cbind(Conc, DV=Pred*(1 + 0.005*rnorm(nRec))) 
s2 = cbind(Conc, DV=Pred*(1 + 0.02*rnorm(nRec)))
s3 = cbind(Conc, DV=Pred*(1 + 0.03*rnorm(nRec)))
s4 = cbind(Conc, DV=Pred*(1 + 0.05*rnorm(nRec)))

rs1 = nlr(fPD3b, s1, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; rs1
rs2 = nlr(fPD3b, s2, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; rs2
rs3 = nlr(fPD3b, s3, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; rs3
rs4 = nlr(fPD3b, s4, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5)) ; rs4

rs1a = nlr(fPD3b, s1, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs1a
rs2a = nlr(fPD3b, s2, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs2a
rs3a = nlr(fPD3b, s3, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs3a
rs4a = nlr(fPD3b, s4, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs4a

s1b = s1
s1b[,"DV"] = 1e9*s1b[,"DV"]
rs1b = nlr(fPD3b, s1b, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2), Error="P") ; rs1b

set.seed(3)
noise = rnorm(nRec)
s1 = cbind(Conc, DV=Pred*(1 + 0.006*noise)) 
s2 = cbind(Conc, DV=Pred*(1 + 0.02*noise))
s3 = cbind(Conc, DV=Pred*(1 + 0.03*noise))
s4 = cbind(Conc, DV=Pred*(1 + 0.05*noise))

rs1 = nlr(fPD3b, s1, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs1
rs2 = nlr(fPD3b, s2, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs2
rs3 = nlr(fPD3b, s3, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs3
rs4 = nlr(fPD3b, s4, c("E0", "Imax", "IC50", "n"), c(175, 35, 120, 2.5), Error="P") ; rs4

fPred = function(Res, Conc)
{
  E0   = Res$Est["PE","E0"]
  Imax = Res$Est["PE","Imax"]
  IC50 = Res$Est["PE","IC50"]
  n    = Res$Est["PE","n"]
  return(E0 - Imax*Conc^n/(IC50^n + Conc^n))
}

DefPar = par(no.readonly=TRUE)
dev.new()
par(mfrow=c(2,2))
plot(Conc, s1[,"DV"], pch=16)
lines(x, fPred(rs1, x))
plot(Conc, s2[,"DV"], pch=16)
lines(x, fPred(rs2, x))
plot(Conc, s3[,"DV"], pch=16)
lines(x, fPred(rs3, x))
plot(Conc, s4[,"DV"], pch=16)
lines(x, fPred(rs4, x))
par(DefPar)

# Figure 3.7
CV = rs1$Est["RSE", c("E0", "Imax", "IC50", "n")]
CV = rbind(CV, rs2$Est["RSE", c("E0", "Imax", "IC50", "n")])
CV = rbind(CV, rs3$Est["RSE", c("E0", "Imax", "IC50", "n")])

dev.new()
ErrLevel = c(0.6, 2, 3)
plot(ErrLevel, CV[,"n"], xlim=c(0, 4), ylim=c(0, 60), type="o", pch=17, xlab="Error Level (%)", ylab="Precision (CV %)")
text(3.2, 41, "n")
lines(ErrLevel, CV[,"IC50"], type="o", pch=15)
text(3.2, 23, "IC50")
lines(ErrLevel, CV[,"Imax"], type="o", pch=18)
text(3.2, 19, "Imax")
lines(ErrLevel, CV[,"E0"], type="o", pch=16)
text(3.2, 2, "E0")
 

# Table 3.4
r2da
r2db
Cols = c("Imax", "IC50", "n", "E0")
Data500 = diag(r2da$Cov)[Cols]/r2da$Est["PE","AddErrVar"]
Data800 = diag(r2db$Cov)[Cols]/r2db$Est["PE","AddErrVar"]
VIF.ratio = Data500/Data800

Tab3.4 = cbind(Estimate=r2db$Est["PE",Cols], Data500, Data800, VIF.ratio)
Tab3.4





