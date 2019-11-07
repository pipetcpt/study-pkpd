setwd("D:/Rt/Gab/")
require(wnl)
require(NonCompart)

dPK07 = read.csv("PK07.csv", skip=1) # Erratum: time unit hr -> min
colnames(dPK07) = c("TIME", "DV"); dPK07

x = dPK07[,"TIME"]
y = dPK07[,"DV"]

## Plotting
plot(x, y, type="o")
plot(x, y, type="o", log="y")

## NCA
ResNCA = sNCA(x, y, adm="Bolus", dose=100, doseUnit="ug", timeUnit="min") ; ResNCA
ResNCA["CLO"]      # CL observed      0.3711596 L/min
ResNCA["MRTIVIFO"] # MRT iv observed  298.4924 min
ResNCA["VSSO"]     # Vss observed     110.7883 L

## Curve stripping
B = exp(ResNCA["b0"]) ; B # 0.7435027
beta = ResNCA["LAMZ"] ; beta # 0.003035501 

logy = log(y)
logyhat = B - beta*x ; logyhat
x2 = x[x < ResNCA["LAMZLL"]] ; x2
logy2 = (logy - logyhat)[x < ResNCA["LAMZLL"]] ; logy2
reslm = lm(logy2 ~ x2) ; reslm # str(reslm)

A = exp(reslm$coefficients[1]) ; A # 0.7514636 
alpha = -reslm$coefficients[2] ; alpha # 0.01097179

fx = function(Time)
{
  return(A*exp(-alpha*Time) + B*exp(-beta*Time))
}
cbind(y, fx(x))

## nonlinear regression
# monoexponential
fPK07a = function(THETA)
{
  A = THETA[1]
  K = THETA[2]
  TIME = e$DATA[,"TIME"]
  Cp = A*exp(-K*TIME)
  return(Cp)
}

r1 = nlr(fPK07a, dPK07, pNames=c("A", "K"), IE=c(2, 0.1), Error="POIS") ; r1
RES1 = e$Residual
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK07a, dPK07, pNames=c("A", "K"), IE=c(1.1, 1.0), Error="POIS")

# biexponential
fPK07b = function(THETA)
{
  A     = THETA[1]
  alpha = THETA[2]
  B     = THETA[3]
  beta  = THETA[4]
  TIME = e$DATA[,"TIME"]
  Cp = A*exp(-alpha*TIME) + B*exp(-beta*TIME)
  return(Cp)
}

r2 = nlr(fPK07b, dPK07, pNames=c("A", "alpha", "B", "beta"), IE=c(1.1, 1.0, 0.2, 0.01), Error="POIS") ; r2
RES2 = e$Residual
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK07b, dPK07, pNames=c("A", "alpha", "B", "beta"), IE=c(1.1, 1.0, 0.2, 0.01), Error="POIS")

# triexponential
fPK07c = function(THETA)
{
  A     = THETA[1]
  alpha = THETA[2]
  B     = THETA[3]
  beta  = THETA[4]
  Cc    = THETA[5]
  gamma = THETA[6]
  TIME = e$DATA[,"TIME"]
  Cp = A*exp(-alpha*TIME) + B*exp(-beta*TIME) + Cc*exp(-gamma*TIME)
  return(Cp)
}

r3 = nlr(fPK07c, dPK07, pNames=c("A", "alpha", "B", "beta", "C", "gamma"), 
         IE=c(1.0, 0.1,  0.7, 0.03, 0.6, 0.005), 
         LB=c(0.5, 0.05, 0.2, 0.01, 0.1, 0.001), 
         UB=c(2.0, 0.2,  1.0, 0.05, 1.0, 0.01), Error="POIS") ; r3
e$r
RES3 = e$Residual
dev.new() ; plot(e$Residual, type="o") ; abline(h=0)
wnl5(fPK07c, dPK07, pNames=c("A", "alpha", "B", "beta", "C", "gamma"), 
         IE=c(1.0, 0.1,  0.7, 0.03, 0.6, 0.005), 
         LB=c(0.5, 0.05, 0.2, 0.01, 0.1, 0.001), 
         UB=c(2.0, 0.2,  1.0, 0.05, 1.0, 0.01), Error="POIS")


dev.new()
matplot(x, cbind(RES1, RES2, RES3), type="o") # Figure 7.2 p510
abline(h=0)
