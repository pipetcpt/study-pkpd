setwd("D:/Rt/Gab")

PRED1i = function(V, K, Dose=10000, TIMEi)
{
  return(Dose/V * exp(-K*TIMEi))
}

PRED2i = function(Cl, V, Dose=10000, TIMEi)
{
  return(Dose/V * exp(-Cl/V*TIMEi))
}

PRED3i = function(AUC, AUMC, Dose=10000, TIMEi)
{
  K = AUC / AUMC
  V = Dose / AUC / K
  return(Dose/V * exp(-K*TIMEi)) 
}

OBJi = function(UCP)
{
  if (e$STEP == "EST") {
    P = DECN(UCP[1:e$nTheta])
  }
  Fi = PRED1i(P[1], P[2], Dose=10000, TIMEi=TIMEi)
  return(nRec*log(UCP[3]) + sum((Yi - Fi)^2)/UCP[3])
}



Data = read.csv("PK01.csv")
Dose = 10000

Subj = unique(Data[,"Subject"])
nSubj = length(Subj)

Result = matrix(nrow=nSubj, ncol=12)
colnames(Result) = c("V_PE", "V_SE", "V_CV", "k_PE", "k_SE", "k_CV", "sgsq_PE", "sgsq_SE", "sgsq_CV", "CL_PE", "CL_SE", "CL_CV")   

#require(NonCompart)
#resNCA = NCA(Data, "Subject", "Time", "Conc", AdmMode="Bolus", Dose=10000)
#resNCA

require(numDeriv)

Lower = double(length=3)
Upper = double(length=3)

grCL = deriv(~V*k,
            c("V","k","sigsq"),
            function.arg=c("V", "k", "sigsq"),
            func=TRUE, hessian=TRUE)

for (i in 1:nSubj) {
  cSubj = Subj[i]
  DATi = Data[Data[,"Subject"]==cSubj, ]
  Yi = DATi[,"Conc"]
  TIMEi = DATi[,"Time"]
  nRec = length(TIMEi)
  Res = optim(par=c(0.1, 0.1, 0.1), fn=OBJi, method="L-BFGS-B")
  Theta = DECN(Res$par[1:e$nTheta])
  Hess = hessian(OBJi, c(Theta, Res$par[3]))
  Ci = solve(Hess/2)
  SE = sqrt(diag(Ci))
  Result[i,1] = Theta[1]
  Result[i,2] = SE[1]
  Result[i,3] = Result[i,2] / Result[i,1] * 100
  Result[i,4] = Theta[2]
  Result[i,5] = SE[2]
  Result[i,6] = Result[i,5] / Result[i,4] * 100
  Result[i,7] = Res$par[3]
  Result[i,8] = SE[3]
  Result[i,9] = Result[i,8] / Result[i,7] * 100

  gr1 = attr(grCL(Res$par[1], Res$par[2]), "gradient")
  Result[i,10] = Res$par[1]*Res$par[2]
  Result[i,11] = sqrt(gr1 %*% Ci %*% t(gr1)) ; SE
  Result[i,12] = Result[i,11] / Result[i,10] * 100
} ; Result


OBJi(Init)

CL = function(V, k)
{
  return(V, k)
}








