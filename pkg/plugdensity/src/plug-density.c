#include <Rmath.h>

void plugin(double *x, int *n, double *z, int *m, double *f, double *h)
{

/************************************************************************
 *	 Version: 1995
 *
 *	 Purpose:
 *
 *	 Simple	 Subroutine for kernel density estimation
 *	 with iterative plug-in bandwidth selection
 *
 *	 This version only uses the gauss kernel and estimates only
 *	 the density itself and not its derivatives.
 *
 *  INPUT:
 *	x[]  double    sorted data array
 *	n    int       length of  x
 *	z[]  double    output grid (sorted array)
 *	m    int       length of z
 *  OUTPUT:
 *	f[]  double    estimated density (array of length m)
 *	h    double    estimated iterative plugin bandwidth
 *
 *
 **************************************************************************/

    const double I_7 = 1./7.;
    const double rt2pi = sqrt(2 * M_PI);
    const double rtpi2 = 2. * M_SQRT_PI;

    const int iter = 5;

    int nx,n2,i,j,it,jbegin,jend;
    double xiqr,h2,h3,s2,s3,d2,d3, rhat2,rhat3,co1,co2,a,s,t;

    /* initializations */

    nx=(*n);
    n2 = nx * nx;
    xiqr=x[(3 * nx)/4 - 1] - x[nx / 4];/* = IQR(x[]) */

    /* estimate inflation constant c */

    h2=(0.920 * xiqr)/ pow(nx, I_7);
    h3=(0.912 * xiqr)/ pow(nx, 1./9.);

    s2 = s3 = 0.;
    for (i = 0; i <= nx-2 ; i++) {
	for (j=i+1 ; j <= nx-1 ; j++) {
	    t = x[i] - x[j];
	    a = t / h2;	    d2= a*a;
	    a = t / h3;	    d3= a*a;
	    if(d2 > 50 && d3 > 60) break;

	    s2 += exp(-d2/2.)*(3. + d2*(-6.+d2));
	    s3 += exp(-d3/2.)*(-15. + d3*(45.+d3*(-15.+d3)));
	}
    }
    rhat2 = 2.*s2/(rt2pi*n2*pow(h2,5)) + 3./(rt2pi*nx*pow(h2,5));
    rhat3 =-2.*s3/(rt2pi*n2*pow(h3,7)) + 15./(rt2pi*nx*pow(h3,7));
    co1= 1.357 * pow(rhat2/rhat3, I_7);
    co2= 1./rtpi2;
    a = 1.132795764/(pow(rhat3, I_7)* sqrt(nx));

/* MM: FIXME?  below we drop all terms  exp(-t^2 / 2) as soon as |t| > 5;
       -----   where exp(- 25/2) is "only" 3.727e-6
 */

    /* loop over iterations */

    for (it = 1; it <= iter; it++) {
	s2=0.;
	for (i = 0; i <= nx-2; i++) {
	    for (j=i+1; j <= nx-1; j++) {
		t = (x[i] - x[j])/a;
		d2= t*t;
		if (d2 > 50) break;
		s2 += exp(-d2/2.)*(3.+d2*(-6.+d2));
	    }
	}
	a = rt2pi * nx * pow(a,5);
	rhat2= 2.*s2/(nx * a) + 3./a;

	/* estimate bandwidth by asymptotic formula */

	a=co1*pow(co2/(rhat2*nx), I_7);
    }
    *h = pow(co2/(rhat2*nx),0.2);

    /* Estimate density f[i] = f_h(z[i]) with plugin bandwidth h : */

    jbegin = jend = 0;
    for (i = 0; i < *m; i++) {
	s = 0.;
	for(j = jbegin; j <= jend ; j++) {
	    t=(z[i] - x[j])/(*h);
	    if(t > 5. && j < nx-1) {
		jbegin++;
		continue;
	    }
	    s += exp(-t*t/2.);
	}
	for (jend = j; jend <= nx-1 ; jend++) {
	    t=(z[i] - x[jend])/(*h);
	    if(t < -5.) break;
	    s += exp(-t*t/2.);
	}
	f[i] = s/(nx*(*h)*rt2pi);
	jend--;
    }
}
