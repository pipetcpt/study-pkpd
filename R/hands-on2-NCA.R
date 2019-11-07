require(wnl)
dPK01 = read.csv("data/hands-on2.csv", skip=1)
colnames(dPK01) = c("TIME", "DV", "ID")
NonCompart::tblNCA(dPK01, key="ID", colTime="TIME", colConc="DV", dose=10, adm="Bolus")

library(ggplot2)
ggplot(dPK01, aes(TIME, DV, group = ID, color = as.factor(ID))) +
  geom_line() + geom_point() + #scale_y_log10() +
  labs(color = 'Subject ID', x = 'Time (hour)', y = 'DV (mg/L)') +
  scale_y_continuous(limits=c(0,1000)) +
  theme_bw()
