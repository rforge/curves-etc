library(cobs99)
options(digits = 6)
postscript("ex3.ps")

data(women) # 15 obs.
attach(women)

## Interpolation! very easy problem (very smooth data!),   BUT :
try( ## gives  "ifl = 5" !!!!!
cobw <- cobs(weight, height, knots = weight, nknots = length(weight))
)

## scobs() "works"  (gives solution with extraneous "waves")

## this works fine
cobw. <- cobs(weight, height, nknots = length(weight), constraint="concave")
cobw.
plot(weight, height); lines(predict(cobw.), col=2); rug(cobw.$knot, col="blue")

## Simple smoothing spline with "full" knots -- gave bug in earliier versions of cobs99:
(con <- rbind(c( 1, weight[1], height[1]),
              c(-1, weight[3], height[3]),
              c( 0, weight[9], height[9])))

n <- length(weight)
coS   <- cobs(weight, height, nknots=n, n.sub=n, lambda = 1)
coS.  <- cobs(weight, height, nknots=n, n.sub=n, lambda = 1,
              pointwise = con)
coSc  <- cobs(weight, height, nknots=n, n.sub=n, lambda = 1, constraint="concave")
coSc. <- cobs(weight, height, nknots=n, n.sub=n, lambda = 1, constraint="concave",
              pointwise = con)

co1S  <- cobs(weight, height, degree=1, nknots=n, n.sub=n, lambda = 1)
## failed in earlier cobs99 versions:
co1S. <- cobs(weight, height, degree=1, nknots=n, n.sub=n, lambda = 1,
             pointwise = con)
## both these two give non-convergence in 300 iterations
co1Sc <- cobs(weight, height, degree=1, nknots=n, n.sub=n, lambda = 1, constraint="concave")
co1Sc.<- cobs(weight, height, degree=1, nknots=n, n.sub=n, lambda = 1, constraint="concave",
              pointwise = con)
## but with large enough 'maxiter', they do converge:
co1Sc <- cobs(weight, height, degree=1, nknots=n, n.sub=n, lambda = 1,
              constraint="concave", maxiter = 50000)
co1Sc.<- cobs(weight, height, degree=1, nknots=n, n.sub=n, lambda = 1,
              constraint="concave", maxiter = 10000, pointwise = con)

round( cbind(resid(coS ), resid(coS. ),
             resid(coSc), resid(coSc.),
             resid(co1S ), resid(co1S. ),
             resid(co1Sc), resid(co1Sc.)
             ), 3)

