### Since we had an undetected bug in  rnorMix()...

## These are defined as norMix() calls in  ../R/zMarrWand-dens.R
library("nor1mix")
ppos <- which("package:nor1mix" == search())
nms <- ls(pat="^MW.nm", pos = ppos)
nms <- nms[order(as.numeric(substring(nms,6)))] # warning <== "MW.nm2.old"

set.seed(123)
for(n in nms) {
    obj <- get(n, pos = ppos)
    cat("\n",n,":\n"); print(obj)
    cat("4 random X from", n,":")
    print(rnorMix(4, obj))
}

## Testing of sort.norMix():
stopifnot(sapply(nms, function(n) {
    o <- get(n)
    if(is.unsorted(o[,"mu"]))
        o <- sort(o)
    identical(o, sort(o))
}))
