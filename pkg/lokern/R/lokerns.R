## lokerns   kernel regression smoothing with local bandwidth selection

lokerns_function(x,y,deriv=0,n.out=300,x.out=NULL,korder=NULL,ihetero=F,
	irnd=T,inputb=F,m1=400,xl=NULL,xu=NULL,s=NULL,sig=NULL,bandwidth=NULL)

{
## x		inputgrid

## y		data

## control and sort inputgrid and data
  n_length(x)
  if (length(y)!=n) stop("Input grid and data must have the same length.")
  sorvec_sort.list(x)
  x_x[sorvec]
  y_y[sorvec]

## deriv          derivative of regression function to be estimated

## n.out length of outputgrid

## compute and sort outputgrid
  if (is.null(x.out)) 
          { 
           n.out_as.integer(n.out)
           x.out_seq(min(x),max(x),length=n.out)
          }
     else {
           sorvec_sort.list(x.out)
           x.out_x.out[sorvec]
           n.out_length(x.out)
          }

## korder         kernel order

## check deriv and korder
   if (is.null(korder)) korder_deriv+2
   if (deriv<0) stop("Order of derivative is negative.")
   if (deriv>4) stop("Order of derivative is too large.")
   if (korder>6) korder_deriv+2

## ihetero 	homo- or heteroszedasticity of error variables

## irnd     	random or non-random t-grid

## inputb 	input bandwidth or estimation of plug-in bandwidth

## m1	  	discretization for integral functional estimation
   if (m1<10) stop("number of discretizations m1 is too small")

## xl           lower bound for integral approximation and variance estimation
## xu           upper bound for integral approximation and variance estimation

  if (is.null(xl)||is.null(xu)) 
       {
        xl_as.double(1)
        xu_as.double(0) 
       }
    
## s		mid-point grid
  if (length(s)!=length(x)+1)
       s_as.double(rep(0,n+1))
  if (is.null(s)) 
       s_as.double(rep(0,n+1))

## sig          input variance
  if (is.null(sig)) sig_as.double(0)

## bandwidth    input bandwidth function
  if (is.null(bandwidth)) 
       {
        inputb_as.integer(0)
        bandwidth_as.double(rep(0,n.out))
       }
  if (length(bandwidth)!= length(x.out)) 
       {
        inputb_as.integer(0)
        bandwidth_as.double(rep(0,n.out))
       }

   if (deriv>2 & as.integer(inputb)==0) stop("Order of derivative is too large.")
   if (korder>4 & as.integer(inputb)==0) korder_deriv+2

## internal parameters and arrays for fortran routine
  len1_as.integer(length(x)+1)
  work1_double(len1*5)
  work2_double(as.integer(m1)*3)
  work3_double(n.out)
  est_double(n.out)
  irnd1_as.integer(1-irnd)


## calling fortran routine
res_.Fortran("lokerns",
           x=as.double(x),
           y=as.double(y),
           as.integer(n),
           x.out=as.double(x.out),
           as.integer(n.out),
           deriv=as.integer(deriv),
           korder=as.integer(korder),
           ihetero=as.integer(ihetero),
           irnd=as.integer(irnd1),
           as.integer(inputb),
           as.integer(m1),
           xl=as.double(xl),
           xu=as.double(xu),
           s=as.double(s),
           sig=as.double(sig),
           as.double(work1),
           as.double(work2),
           as.double(work3),
           bandwidth=as.double(bandwidth),         
           est=as.double(est)
           )

  return(x=x,y=y,bandwidth=res$bandwidth,x.out=x.out,est=res$est,sig=res$sig,
	deriv=deriv,korder=korder,xl=res$xl,xu=res$xu,s=res$s)
}
