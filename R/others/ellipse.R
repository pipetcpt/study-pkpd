ellipse = function(center=c(0,0), radius=c(2,1), alpha=0, npoints=100)
{
  theta = seq(0, 2*pi, length=npoints + 1)
  x0 = radius[1]*cos(theta)
  y0 = radius[2]*sin(theta)
  x = x0 * cos(alpha) - y0 * sin(alpha) + center[1]
  y = x0 * sin(alpha) + y0 * cos(alpha) + center[2]
  LongAxis = max(radius)
  xlm = center[1] + c(-1,1)*LongAxis 
  ylm = center[2] + c(-1,1)*LongAxis
  plot(x, y, xlim=xlm, ylim=ylm, type="l")
  return(cbind(x, y))
}

dev.new()
ellipse()

dev.new()
ellipse(c(0, 0), c(3, 2), pi/4)

dev.new()
ellipse(c(1, 1), c(3, 2), pi/4)


ellipRange = function(center=c(0,0), radius=c(2,1), alpha=0)
{
  x0max = sqrt(radius[1]^2*cos(alpha)^2 + radius[2]^2*sin(alpha)^2) 
  y0max = sqrt(radius[1]^2*sin(alpha)^2 + radius[2]^2*cos(alpha)^2)
  xrange = center[1] + c(-1,1)*x0max
  yrange = center[2] + c(-1,1)*y0max
  return(cbind(xrange, yrange))
}

ellipRange()
ellipRange(c(0, 0), c(3, 2), pi/4)
ellipRange(c(1, 1), c(3, 2), pi/4)
ellipRange(c(1, 1), c(3, 2), pi/6)


############# Drawing using sampled data and distribution
ci = 0.95
npoints = 100
dimR = 2 # plotting dimension
npara = 5 # mu: 2, mCov: 3

mu0 = c(0, 0)
mCov0 = matrix(c(1,0.5,0.5,1), nrow=2)
eg0 = eigen(mCov0)
alpha0 = atan(eg0$vectors[2,1]/eg0$vectors[1,1])
radius0 = sqrt(qchisq(ci, dimR)*eg0$values) # normal distribution, when plot with mu0, mCov0

theta = seq(0, 2*pi, length=100 + 1)

x0 = radius0[1]*cos(theta)
y0 = radius0[2]*sin(theta)
x1 = x0 * cos(alpha0) - y0 * sin(alpha0) + mu0[1]
y1 = x0 * sin(alpha0) + y0 * cos(alpha0) + mu0[2]
LongAxis = max(radius0)
xlm = mu0[1] + c(-1,1)*LongAxis
ylm = mu0[2] + c(-1,1)*LongAxis
plot(x1, y1, xlim=xlm, ylim=ylm, type="l", col="red")
points(mu0[1], mu0[2], pch="+", col="red")

library(MASS)
Data = mvrnorm(npoints, mu0, mCov0) ; Data
points(Data)

mu = colMeans(Data) ; mu
mCov = cov(Data) ; mCov
eg = eigen(mCov) # sorted as descening order
alpha = atan(eg$vectors[2,1]/eg$vectors[1,1]) ; alpha*180/pi
radius = sqrt(dimR*qf(ci, dimR, npoints - npara)*eg$values) # for sample data, when plot with estimated mu, mCov

x2 = radius[1]*cos(theta)
y2 = radius[2]*sin(theta)
x3 = x2 * cos(alpha) - y2 * sin(alpha) + mu[1]
y3 = x2 * sin(alpha) + y2 * cos(alpha) + mu[2]
lines(x3, y3, type="l")
points(mu[1], mu[2], pch="+")

x0max = sqrt(radius[1]^2*cos(alpha)^2 + radius[2]^2*sin(alpha)^2) 
y0max = sqrt(radius[1]^2*sin(alpha)^2 + radius[2]^2*cos(alpha)^2)
xrange = mu[1] + c(-1,1)*x0max
yrange = mu[2] + c(-1,1)*y0max
cbind(xrange, yrange)


## using rgl
library(rgl)
npoints = 101
x = seq(-3, 3, length=npoints)
y = x

mvdnorm = function(v, mu=c(0,0), mCov=matrix(c(1,0,0,1),nrow=2))
{
  n = length(mu)
  return(as.numeric((2*pi*det(mCov))^(-n/2)*exp(-0.5 * t(v - mu) %*% solve(mCov) %*% (v - mu))))
}

mvdnorm(c(1,1), mu0, mCov0)

z = matrix(nrow=npoints, ncol=npoints)

for (i in 1:101) {
  for (j in 1:101) {
    z[i,j] = mvdnorm(c(x[i],y[j]), mu=mu0, mCov=mCov0)
  }
}
open3d()
clear3d("all")
bg3d(color="#887777")
light3d()
zscale = 10
surface3d(x, y, z*zscale, color="#FF2222",alpha=0.5)

z3 = matrix(rep(0, npoints*npoints), nrow=npoints)
spheres3d(Data[,1], Data[,2], z3, radius=0.1,color="#CCCCFF")

for (i in 1:npoints) {
 for (j in 1:npoints) {
   z3[i,j] = mvdnorm(c(x1[i],y1[j]), mu=mu0, mCov=mCov0)
  }
}
lines3d(x1, y1, z3)


z4 = matrix(rep(0.1, npoints*npoints), nrow=npoints)
surface3d(x, y, z4, color="#000000")


