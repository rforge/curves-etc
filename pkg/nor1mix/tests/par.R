library(nor1mix)

## Check  nM2par(), par2norMix() and llnorMix() :

nms <- paste("MW.nm", 1:16, sep="")
for(n in nms) {
    cat(n,":")
    obj <- get(n, envir = as.environment("package:nor1mix"))
    xx <- rnorMix(1000, obj)
    logLik.x <- sum(dnorMix(xx, obj, log = TRUE))
    pp   <- nM2par(obj) # use "current default"
    pp.l <- nM2par(obj, trafo = "logit")
    pp.c <- nM2par(obj, trafo = "clr1")
    nm   <- par2norMix(pp) # (current) default
    nm.l <- par2norMix(pp.l, trafo= "logit")
    nm.c <- par2norMix(pp.c, trafo=  "clr1")
    stopifnot(exprs = {
        all.equal(pp,   nM2par(nm),                   tol= 1e-15)
        all.equal(pp.l, nM2par(nm.l, trafo= "logit"), tol= 1e-15)
        all.equal(pp.c, nM2par(nm.c, trafo=  "clr1"), tol= 1e-15)
        all.equal(obj, nm,   check.attributes=FALSE, tol=4e-15)
        all.equal(obj, nm.l, check.attributes=FALSE, tol=4e-15)
        all.equal(obj, nm.c, check.attributes=FALSE, tol=4e-15)
        ## xx
        all.equal(llnorMix(pp  , xx),                 logLik.x, tol = 1e-15)
        all.equal(llnorMix(pp.l, xx, trafo= "logit"), logLik.x, tol = 1e-15)
        all.equal(llnorMix(pp.c, xx, trafo=  "clr1"), logLik.x, tol = 1e-15)
    })
    cat(" [ok]\n")
}
