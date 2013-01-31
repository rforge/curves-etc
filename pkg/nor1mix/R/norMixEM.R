#
#  Copyright (C) 2010 Friedrich Leisch
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


norMixEM <- function(x, m, name=NULL,
                     iter.max=100, tol=sqrt(.Machine$double.eps))
{
    MYCALL <- match.call()
    if(is.null(name)) name <- deparse(MYCALL)
    nx <- length(x)

    m <- as.integer(m)
    if(max(m)<2) stop("Number of clusters must be larger or equal to 2")
    
    if(length(m)>1){
        init <- rep(m, length=nx)
        m <- max(init)
    }
    else{
        q <- quantile(x, seq(0, 1, by=1/m))
        q[1] <- q[1]-0.1
        init <- as.integer(cut(x, q))
    }
        
    iter.max <- as.integer(iter.max)[1]
    
    post <- matrix(0, nrow=nx, ncol=m)
    post[cbind(1:nx, init)] <- 1

    mu <- sd <- double(m)

    llh <- -Inf
    for(n in 1:iter.max)
    {
        prior <- colMeans(post)
        llh.old <- llh
        for(l in 1:m)
        {
            sp <- sum(post[,l])
            mu[l] <- sum(post[,l]*x)/sp
            sd[l] <- pmax(sqrt( sum(post[,l] * (x-mu[l])^2 ) / sp),
                          tol)
            post[,l] <- prior[l] * dnorm(x, mu[l], sd[l])
        }
        llh <- sum(log(rowSums(post)))
        if(abs(llh-llh.old)/(abs(llh)+0.1) < tol) break
        
        post <- post/rowSums(post)
    }
    norMix(mu=mu, sig2=sd^2, w=prior, name=name)
}

