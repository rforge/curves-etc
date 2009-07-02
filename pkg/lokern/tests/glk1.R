require(lokern)
data(xSim)
xSim # `test' for the dataset

n <- length(xSim)
tt <- ((1:n) - 1/2)/n # equidistant x

str(gk <- glkerns(tt, xSim))
summary(gk$est)
gk$bandwidth
glkerns(tt,xSim, deriv = 1)$bandwidth
glkerns(tt,xSim, deriv = 2)$bandwidth

demo("glk-derivs")

p.3glks(tt, xSim, kord = 3)

p.3glks(tt, xSim, kord = 4, useB = 0.15)

str(p.3glks(tt, xSim, kord = 5, useB = 0.12) ) # k.ord = (4,5,4) => less sensiacl?

p.3glks(tt, xSim, kord = 6, useB = 0.2, derivs = 0:3) # k.ord = (6,5,6, 5)

## "FIXME" visually compare with numerical derivatives (e.g. from splines).
