
ellipsem = function (mu, amat, c2, npoints = 100, showcentre = TRUE, ...)
{
  if (all(dim(amat) == c(2, 2))) {
    eamat = eigen(amat)
    hlen = sqrt(c2/eamat$val)
    theta = angle(eamat$vec[1, 1], eamat$vec[2, 1])
    ellipse(hlen[1], hlen[2], theta, mu[1], mu[2], npoints = npoints, ...)
    if (showcentre) points(mu[1], mu[2], pch = 3)
  }
#  invisible()
}

ellipse = function (hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, yc = 0, newplot = TRUE, npoints = 100, ...)
{
  a = seq(0, 2 * pi, length = npoints + 1)
  x = hlaxa * cos(a)
  y = hlaxb * sin(a)
  alpha = angle(x, y)
  rad = sqrt(x^2 + y^2)
  xp = rad * cos(alpha + theta) + xc
  yp = rad * sin(alpha + theta) + yc
  if (newplot) plot(xp, yp, type = "l", ...)
  else lines(xp, yp, ...)
#  invisible()
}

angle = function (x, y)
{
  angle2 = function(xy) {
    x <- xy[1]
    y <- xy[2]
    if (x > 0) {
      atan(y/x)
    } else {
      if (x < 0 & y != 0) {
        atan(y/x) + sign(y) * pi
      } else {
        if (x < 0 & y == 0) {
          pi
        } else {
          if (y != 0) {
            (sign(y) * pi)/2
          } else {
            NA
          }
        }
      }
    }
  }
  apply(cbind(x, y), 1, angle2)
}

x = rnorm(100)
y = rnorm(100)

Data = cbind(x, y)
mu = colMeans(Data)
amat = cov(Data)
ci = 0.95
c2 = -2*log(1 - ci) ; c2
dev.new()
ellipsem(mu, amat, c2)
points(x,y)

 