setwd("D:/Rt/Gab/")
require(NonCompart)

dPK15 = read.csv("PK15.csv", skip=1, as.is=TRUE)
colnames(dPK15) = c("DoseGrp", "Conc", "Period", "ID", "Time", "Dose", "CoD", "Sex") ; dPK15
dPK15 = dPK15[dPK15$Conc != "missing" & dPK15$DoseGrp != 0,]
dPK15[,"Conc"] = as.numeric(dPK15[,"Conc"])
dPK15 = dPK15[order(dPK15$DoseGrp, dPK15$ID, dPK15$Period, dPK15$Time),] ; dPK15

DoseGrps = unique(dPK15$DoseGrp) ; DoseGrps
nDoseGrp = length(DoseGrps) ; nDoseGrp

IDs = unique(dPK15$ID) ; IDs
nID = length(IDs) ; nID

Periods = unique(dPK15$Period) ; Periods
nPeriod = length(Periods) ; nPeriod

oldpar = par(no.readonly=TRUE)

dev.new()
par(mfrow = c(3,2))

cSex = c("male", "female")
nSex = 2
cPeriod = c(11, 26, 52)
nPer = 3

Col = c("red", "blue", "green")
names(Col) = c(10, 56, 320)

NCAs = unique(dPK15[dPK15$Period %in% c(11, 26, 52),c("DoseGrp", "ID", "Sex", "Period")]) ; NCAs
resNCA = cbind(NCAs, AUClast=0, Cmax=0) ; resNCA

for (j in 1:nPer) {
  for (i in 1:nSex) {
    Data = dPK15[dPK15$Sex==cSex[i] & dPK15$Period == cPeriod[j], ]
    IDs = unique(Data$ID)
    nID = length(IDs)
    xlm = c(0, 24)
    ylm = c(0.01, max(Data[,"Conc"], na.rm=TRUE))
    plot(xlm, ylm, type="n", ylab="Concentration (uM)", log="y", main=paste0("Period=",cPeriod[j],"w, ", cSex[i]))
    for (k in 1:nID) {
      iData = Data[Data$ID == IDs[k],]
      x = iData[,"Time"]
      y = iData[,"Conc"]
      lines(x, y, col=Col[as.character(unique(iData[,"DoseGrp"]))])
      rNCA = sNCA(x, y)
      resNCA[resNCA$ID == IDs[k] & resNCA$Period==cPeriod[j],"AUClast"] = rNCA["AUCLST"]
      resNCA[resNCA$ID == IDs[k] & resNCA$Period==cPeriod[j],"Cmax"] = rNCA["CMAX"]
    }
  }
} ; resNCA

tSum = unique(resNCA[,c("DoseGrp", "Sex", "Period")])
tSum = cbind(tSum, AUC.mean=0, AUC.sd=0, Cmax.mean=0, Cmax.sd=0) ; tSum
for (i in 1:nrow(tSum)) {
  tDoseGrp = tSum[i, "DoseGrp"]
  tSex     = tSum[i,"Sex"]
  tPeriod  = tSum[i, "Period"]
  tSum[i,"AUC.mean"]  = mean(resNCA[resNCA$DoseGrp == tDoseGrp & resNCA$Sex == tSex & resNCA$Period == tPeriod, "AUClast"])
  tSum[i,"AUC.sd"]    = sd(resNCA[resNCA$DoseGrp == tDoseGrp & resNCA$Sex == tSex & resNCA$Period == tPeriod, "AUClast"])
  tSum[i,"Cmax.mean"] = mean(resNCA[resNCA$DoseGrp == tDoseGrp & resNCA$Sex == tSex & resNCA$Period == tPeriod, "Cmax"])
  tSum[i,"Cmax.sd"]   = sd(resNCA[resNCA$DoseGrp == tDoseGrp & resNCA$Sex == tSex & resNCA$Period == tPeriod, "Cmax"])
} ; tSum


dev.new()
par(mfrow = c(2,1))

cPlot = c("AUC.mean", "Cmax.mean")
nPlot = length(cPlot)
for (i in 1:nPlot) {
  xlm = c(0, max(DoseGrps))
  ylm = c(0, max(tSum[,cPlot[i]]))
  plot(xlm, ylm, type="n", xlab="Dose", ylab=strsplit(cPlot[i],"[.]")[[1]][1])
  for (j in 1:nPer) {
    for (k in 1:nSex) {
      x = tSum[tSum$Period == cPeriod[j] & tSum$Sex == cSex[k], "DoseGrp"]
      y = tSum[tSum$Period == cPeriod[j] & tSum$Sex == cSex[k], cPlot[i]]
      lines(x, y, lty=k, col=j)
    }
  }
  legend(0, ylm[2], t(outer(cPeriod, cSex, paste, sep="w-")), col=sort(rep(1:nPer,nSex)), lty=1:nSex, cex=0.8)
}

## Second version for plot
dev.new()
par(mfrow = c(2,1))

cPlot = c("AUC.mean", "Cmax.mean")
nPlot = length(cPlot)
for (i in 1:nPlot) {
  xlm = c(0, max(DoseGrps))
  ylm = c(0, max(tSum[,cPlot[i]]))
  plot(xlm, ylm, type="n", xlab="Dose", ylab=strsplit(cPlot[i],"[.]")[[1]][1])
  vLegChr = vector()
  vLegCol = vector()
  vLegLty = vector()
  for (j in 1:nPer) {
    for (k in 1:nSex) {
      x = tSum[tSum$Period == cPeriod[j] & tSum$Sex == cSex[k], "DoseGrp"]
      y = tSum[tSum$Period == cPeriod[j] & tSum$Sex == cSex[k], cPlot[i]]
      lines(x, y, lty=k, col=j)
      vLegChr = c(vLegChr, paste0(cPeriod[j],"w-",cSex[k]))
      vLegCol = c(vLegCol, j)
      vLegLty = c(vLegLty, k)
    }
  }
  legend(0, ylm[2], vLegChr, col=vLegCol, lty=vLegLty, cex=0.8)
}


par(oldpar)
