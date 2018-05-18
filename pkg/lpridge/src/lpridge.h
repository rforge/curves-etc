/* .Fortran calls */

#define FDEF(name)  {#name, (DL_FUNC) &F77_NAME(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

/*     {"lpepa_s",   (DL_FUNC) &F77_NAME(lpepa_s),   16},

      subroutine lpepa_s(t,x,n, b,nue,p, tt,m,mnew,imoms,moms,y,
     .              leng,nmoms,nvar,var)
      integer n,nue,p,m,mnew,imoms(*),leng,nmoms,nvar
      double precision t(*),x(*),b(*),tt(*),moms(nmoms,0:*),y(*),var(*)
*/

// lpepa.f :
void F77_NAME(lpepa_s)(
    const double t[], const double x[], const int *n,
    double b[], int *nue, int *p,
    const double tt[], const int *m, const int *mnew,
    int imoms[], double moms[], double y[],
    int *leng, int *nmoms,
    int *nvar, double var[]);

// lpridge.f :
void F77_NAME(lpridge_s)(
    const double t[], const double x[], const int *n,
    double b[], int *nue, int *p, int *kord,
    const double wk[], const double tt[], const int *m, const int *mnew,
    int imoms[], double moms[], double y[], int *leng, int *nmoms,
    int *nvar, double var[], double *ridge, int *nsins);
