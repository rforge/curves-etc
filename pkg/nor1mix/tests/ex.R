## Marron Wand examples are defined as norMix() calls in  ../R/zMarrWand-dens.R
library("nor1mix")

ii <- c(32, 39, 40, 40, 48, 48, 49, 56, 57, 58, 65)
pp <- 0.486 + ii / 100000 # very constrained set

e <- norMix(mu = c(-0.825,0.275), sig2 = c(0.773,0.773), w = c(1, 3)/4)

qnorMix(pp, e, trace = TRUE)
## failed for version <= 1.0-5

q.pp <- -c(7.78151762922529, 7.60511100150266, 7.57991031275271, 7.57991031275271,
           7.37830712226037, 7.37830712226037, 7.35310701314534, 7.17670804948685,
           7.15150845450588, 7.12630892371882, 6.94991400429199) / 1000

for (m in eval(formals(qnorMix)$method)) {
    cat("method ", m,":")
    stopifnot(all.equal(q.pp, qnorMix(pp, e, method = m, tol = 1e-14),
			tol = 1e-13),# 1.022e-14 (32-bit)
              abs(qnorMix(rep(1/2, 8), MW.nm10, method = m)) < 1e-14
              )
    cat("\n")
}

### a "nasty" example (for the newton steps):
ip <- c(0, 1e-11, 3.33e-09, 7.705488e-05, 0.0001670041, 0.00125378934,
        0.00141169953, 0.00357409109, 0.00644073795, 0.00853238955, 0.01361442421,
        0.01672761627, 0.02067755849, 0.02124092026, 0.03327537558, 0.03527226553,
        0.05365983941, 0.05482289811, 0.05669602608, 0.05982167629)
qv <- qnorMix(1-ip, MW.nm12, trace=1)
## now ok
##--> divergence in Newton
## Error in if (relErr < tol) break : missing value where TRUE/FALSE needed
## qv <- qnorMix(1-ip, MW.nm12, trace=2)

### FIXME --- lower.tail=FALSE  does not work correctly (now *somewhat* works)
qv. <- qnorMix(ip, MW.nm12, lower.tail=FALSE, trace=2, maxiter=50)
## -> infinite loop! --- now "just wrong"
stopifnot(all.equal(qv, qv.))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
