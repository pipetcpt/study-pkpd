# qchisq(p, df) becomes asymptocially df*qf(p, df, df2) as df2 goes to infinity

compChisqF = function(dfs, dimR=2, ci=0.95)
{
  qchi = qchisq(ci, dimR)
  ndf = length(dfs)
  res = matrix(nrow=ndf, ncol=3)
  for (i in 1:ndf) {
    res[i, 1] = dimR*qf(ci, dimR, dfs[i])
    res[i, 2] = res[i, 1] - qchi
    res[i, 3] = res[i, 2] / res[i, 1] * 100 
  }
  res = cbind(dfs,res)
  colnames(res) = c("df", "dimR*qf", "difference", "pct difference")
  return(res)
}

dimR = 2
ci = 0.95
dfs = c(1:30, 100, 1000, 10000, 100000)

qchisq(ci, dimR)
compChisqF(dfs)

qchisq(ci, 3)
compChisqF(dfs, dimR=3)

qchisq(ci, 10)
compChisqF(dfs, dimR=10)

qchisq(ci, 100)
compChisqF(dfs, dimR=100)
