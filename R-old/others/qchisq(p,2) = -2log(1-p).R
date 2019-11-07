p = seq(0, 1, length=101)
res = matrix(nrow=101, ncol=2)
for (i in 1:101) {
  res[i,1] = -2*log(1 - p[i])
  res[i,2] = qchisq(p[i], 2)
} ; res

