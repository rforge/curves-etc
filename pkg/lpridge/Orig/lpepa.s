## lpepa	local polynomials with Epanechnikov weights 
##		for regression functions and derivatives 

lpepa_function(x,y,bandwidth,deriv=0,n.out=100,x.out=NULL,order=1,mnew=100,var=F)

{
## x		inputgrid
## y		data

## control and sort inputgrid and data
  n_length(x)
  if (length(y)!=n) stop("Input grid and data must have the same length.")
  sorvec_sort.list(x)
  x_x[sorvec]
  y_y[sorvec]

## bandwidth	bandwidth for estimation

## deriv	derivative of regression function to be estimated

## n.out	length of outputgrid
## x.out	outputgrid

## compute and control outputgrid
  if (is.null(x.out)) {
    n.out_as.integer(n.out)
    x.out_seq(min(x),max(x),length=n.out)
  }
  else {
    n.out_length(x.out)
  }

## compute vector of bandwiths
  if (length(bandwidth)==1)
    bandwidth_as.double(rep(bandwidth,n.out))
  if (length(bandwidth)!=n.out)
    stop("Length of bandwith is not equal to length of output grid.")

## sort outputgrid and bandwidth
  sorvec_sort.list(x.out)
  x.out_x.out[sorvec]
  bandwidth_bandwidth[sorvec]

## order	order of local polynomial approximation
## check order, deriv
  if (order<0) stop("Polynomial order is negative.")
  if (deriv>order) stop("Order of derivative is larger than polynomial order.")

## mnew		force of restart

## var		switch for variance estimation

## internal parameters and arrays (see fortran routine)
  leng_10
  nmoms_as.integer(length(x)/leng+1)
  imoms_integer(nmoms)
  moms_double(nmoms*4*(2+order+as.integer(var)))

## check internal limitations from fortran routine
  if ((2+order)>12) stop("Polynomial order exceeds 10.")
  if ((as.integer(var)<0)||(as.integer(var)>1)) stop("var has wrong values")


  res_.Fortran("lpepa",
	   x=as.double(x),
	   y=as.double(y),
	   as.integer(n),
	   bandwidth=as.double(bandwidth),
	   deriv=as.integer(deriv),
	   order=as.integer(order),
	   x.out=as.double(x.out),
	   as.integer(n.out),
	   mnew=as.integer(mnew),
	   as.integer(imoms),
	   as.double(moms),
	   est=double(n.out),
	   as.integer(leng),
	   as.integer(nmoms),
	   var=as.integer(var),
	   est.var=double(n.out)
	   )

  return(x=x,y=y,bandwidth=res$bandwidth,deriv=deriv,x.out=x.out,
     order=order,mnew=mnew,var=var,est=res$est,est.var=res$est.var)
}
	    
