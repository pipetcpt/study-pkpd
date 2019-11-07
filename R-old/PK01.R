require(wnl)
dPK01 = read.csv("data/PK01.csv", skip=1)
colnames(dPK01) = c("TIME", "DV", "ID") ; dPK01

## NCA
library(NonCompart)
tblNCA(dPK01, key="ID", colTime="TIME", colConc="DV", dose=10, adm="Bolus")
###


IDs = unique(dPK01[,"ID"])
nID = length(IDs)

DOSE = 10000 # ug

fPK01 = function(THETA) # Prediction function
{
  V  = THETA[1]
  K  = THETA[2]
  Cp = DOSE/V*exp(-K*TIME)  # External DOSE, TIME, eq 1:2
  return(Cp)
}

Result = vector()
for (i in 1:nID) {
  cID = IDs[i]
  Data = dPK01[dPK01$ID == cID,]
  TIME = dPK01[dPK01$ID == cID,"TIME"]
  Res = nlr(fPK01, Data, pNames=c("V", "k"), IE=c(20, 0.2),
            SecNames=c("CL", "AUC", "AUMC" , "Thalf", "MRT"), 
            SecForms=c(~V*k, ~DOSE/V/k, ~DOSE/V/k/k, ~log(2)/k, ~1/k))
  print(paste0("### ID = ", cID, " ###"))
  print(Res)
  Result = rbind(Result, cbind(ID=cID, Res$Est))
} ; Result

# Figure 1.1, p 470
plot(0, 1, type="n", xlim=c(0, 160), ylim=c(10, 1000), log="y", xlab="Time (min)", ylab="Concentration (ug/L)")
for (i in 1:nID) {
  cID = IDs[i]
  TIME = dPK01[dPK01$ID == cID,"TIME"]
  points(TIME, dPK01[dPK01$ID == cID,"DV"], pch=14+i)
  cTHETA = Result[Result[,"ID"]==cID & rownames(Result)=="PE", c("V", "k")]
  lines(TIME, fPK01(cTHETA))   
}
###
