#include <math.h>

double plugin(double x[],int n,double z[],int m,double f[])
{

/************************************************************************
*       Version: 1995
*                                                                       
*       Purpose:                                                        
*                                                                       
*       Simple  Subroutine for kernel density estimation 
*       with iterative plug-in bandwidth selection
*       
*       This version only uses the gauss kernel and estimates only
*       the density itself and not its derivatives.
*                                                                       
*  INPUT    x       DOUBLE    sorted data array
*  INPUT    n       INT       length of  x
*  INPUT    z       DOUBLE    output grid (sorted array)
*  INPUT    m       INT       length of z
*  OUTPUT   f       DOUBLE    estimated density (array of length m)
*  OUTPUT   plugin  DOUBLE    estimated iterative plugin bandwidth
*                                                                       
*                                                                       
**************************************************************************/

int i,j,iter=5,it,jbegin=0,jend=0;
double xn,xiqr,rtpi2,rt2pi,h2,h3,s2=0.0,s3=0.0,d2,d3,e2,e3,
       rhat2,rhat3,co1,co2,a,s,t,h;

/* initializations */

xn=n;
i=floor(0.75*n)-1;
xiqr=x[i];
i=floor(0.25*xn);
xiqr-=x[i];
rt2pi=sqrt(6.28318529);
rtpi2=sqrt(3.141592645)*2.0;

/* estimate inflation constant c */

h2=(0.920*xiqr)/pow(xn,1.0/7.0);
h3=(0.912*xiqr)/pow(xn,1.0/9.0);

for (i=0;i<=n-2;i++)
{  for (j=i+1;j<=n-1;j++)
   {  d2=pow((x[i]-x[j])/h2,2);
      d3=pow((x[i]-x[j])/h3,2);
      if((d2>50)&&(d3>60)) break;
      e2=exp(-d2/2.0);
      e3=exp(-d3/2.0);
      s2+=e2*(3.0+d2*(-6.0+d2));
      s3+=e3*(-15.0+d3*(45.0+d3*(-15.0+d3)));
   }
}
rhat2=2.0*s2/(rt2pi*pow(xn,2)*pow(h2,5));
rhat2+=3.0/(rt2pi*xn*pow(h2,5));
rhat3=-2.0*s3/(rt2pi*pow(xn,2)*pow(h3,7));
rhat3+=15.0/(rt2pi*xn*pow(h3,7));
co1=1.357*pow(rhat2/rhat3,1.0/7.0);
co2=1.0/rtpi2;
a=1.132795764/(pow(rhat3,1.0/7.0)*pow(xn,0.5));

/* loop over iterations */

for (it=1;it <=iter;it++)
{  s2=0.0;
   for (i=0;i<=n-2;i++)
   {  for (j=i+1;j<=n-1;j++)
      {  d2=pow((x[i]-x[j])/a,2);
         if (d2>50) break;
         e2=exp(-d2/2.0);
         s2+=e2*(3.0+d2*(-6.0+d2));
      }
    }
   rhat2=2.0*s2/(rt2pi*pow(xn,2)*pow(a,5));
   rhat2+=3.0/(rt2pi*xn*pow(a,5));

   /* estimate bandwidth by asymptotic formula */

   a=co1*pow(co2/(rhat2*xn),1.0/7.0);
}
h=pow(co2/(rhat2*xn),0.2);

/* estimate density with plugin bandwidth */

for (i=0;i<=m-1;i++)
{  s=0.0;
   for(j=jbegin;j<=jend;j++)
   {   t=(z[i]-x[j])/h;
       if((t>5.0)&&(j<n-1))
       {  jbegin++;
          continue;
       }
       s+=exp(-t*t/2.0);
   }
   for (jend=j;jend <=n-1;jend++)
   {  t=(z[i]-x[jend])/h;
      if(t<-5.0) break;
      s+=exp(-t*t/2.0);
   }
   f[i]=s/(xn*h*rt2pi);
   jend--;
}
return h;
}
