## From: Ravi Varadhan <rvaradhan@jhmi.edu>
## To: maechler@stat.math.ethz.ch
## Subject: A potential "bug" report - An issue with "lokern" package
## Date: Tue, 03 Mar 2009 12:43:36 -0500

## Dear Martin,

## I was trying to run the glkerns() and lokerns() functions in the "lokern" package for the Old Faithful Geyser data from "MASS".  Both functions just freeze-up.  I think that the problem lies in "x" values being not unique for this data (some values are replicated).

## However, if run with is.rand = FALSE, then everything is ok.  It would useful, I think, if you have a check in the code that checks for dupicity in X-values and then either provides a warning to change the "is.rand" option or just go ahead and change that inside the code.

data(geyser, package = "MASS")

x <- geyser$duration
y <- geyser$waiting
plot(x, y)

require(lokern)
ord <- order(x)
x <- x[ord]
y <- y[ord]
glkerns(x, y, is.rand=FALSE) # is fine
glkerns(x, y) # freezes-up


## Thank you,
## Ravi.

## ____________________________________________________________________

## Ravi Varadhan, Ph.D.
## Assistant Professor,
## Division of Geriatric Medicine and Gerontology
## School of Medicine
## Johns Hopkins University

## Ph. (410) 502-2619
## email: rvaradhan@jhmi.edu

