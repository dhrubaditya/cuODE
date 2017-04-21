#include <stdio.h>
#include <math.h>
#include "SpecialFunctions.h"
/*---------------------------------------*/
__device__ double LnGamma(double zp1){
  /* Calculates the logarithm of Gamma function using Lancsoz algorithm, 
     following the code given in Numerical Recepies (fortran) page 206-207 
     \Gamma(zz+1) = (zz+gamma+1/2)^{z+1/2}e^{-(z+gamma+1/2)}\sqrt{2\pi}[c_0+\sum_{j=1,6} \frac{c_j}{z+j}]
     ln\Gamma(zz+1) = (zz+1/2)ln(zz+gamma+1/2)+[-(z+gamma+1/2)] + (0.5)\ln(2\pi)+\ln[c_0+\sum_{j=1,6} \frac{c_j}{z+j}]]
     ln\Gamma(zp1) = (zp1-1/2)ln(zp1-1/2+gamma)+[-(zp1-1/2+gamma)] + (0.5)\ln(2\pi)+\ln[c_0+\sum_{j=1,6} \frac{c_j}{zp1+j-1}]]
*/
  double LG;
  int const N=6;
  double zz=zp1-1;
  if (zz >= 0 ){
    double pi;
    pi = 4.*atan(1.);
    double sqrt2pi=sqrt(2.*pi);
    double cof[N+1];
    cof[0]=1.000000000190015;
    cof[1]=76.18009171947146;
    cof[2]=-86.50532032941677;
    cof[3]=24.01409824083091;
    cof[4]=-1.231739572450155;
    cof[5]=0.1208650973866179e-2;
    cof[6]=-0.5395239384953e-5;
    double gamma=5.;
    double sum=0.;
    for (int j=1; j<=N; j++){
      sum = sum+cof[j]/(zz+(double)(j));
    }
    sum = sum + cof[0];
    double lnsum=log(sqrt2pi*sum);
    LG= (zz+0.5)*log(zz+gamma+0.5)-(zz+gamma+0.5) +lnsum ;
  }else{
    printf("Routine calculates Log Gamma function only for positive argument\n");
    LG=0.;
  }
    return LG;
}
/*---------------------------------------*/

