require(lokern)
data(xSim)
n <- length(xSim)
tt <- ((1:n) - 1/2)/n # equidistant x

str(gk <- glkerns(tt, xSim))
summary(gk$est)
gk$bandwidth
glkerns(tt,xSim, deriv = 1)$bandwidth
glkerns(tt,xSim, deriv = 2)$bandwidth

