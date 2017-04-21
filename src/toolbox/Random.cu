#include <stdio.h>
#include <math.h>
#include "SpecialFunctions.h"
#include "Random.h"
/*------------ Uses CUDA random number generator -------- */
__global__ void init_random(unsigned long long *seed, curandState  *global_state){
  int tid = blockIdx.x;
  unsigned long long local_seed = seed[tid];
  curandState local_state;
  local_state = global_state[tid];
  curand_init(local_seed,tid,0, &local_state);
  global_state[tid] = local_state;
}
/*--------------------------------------*/
__global__ void random(double *x, curandState *global_state){
  int tid =  blockIdx.x;
  curandState local_state;
  local_state = global_state[tid];
  x[tid] = (double) curand(&local_state);
  global_state[tid] = local_state;
}
/*--------------------------------------*/
__global__ void UniformRandom(double *x, curandState *global_state){
  int tid =  blockIdx.x;
  curandState local_state;
  local_state = global_state[tid];
  x[tid] = (double) curand_uniform(&local_state);
  global_state[tid] = local_state;
}
/*--------------------------------------*/
__device__ double Gaussian(double mean, double sigma, curandState *mystate){
  double xx= (double) curand_normal(mystate);
  double yy=mean+sigma*xx;
  return yy;
}
/*--------------------------------------*/
__device__ double Poisson(double xmean, curandState *mystate){
  double reject_factor=0.9,reject;
  double pi;
  pi = 4.*atan(1.);
  double x,xcomp;
  if (xmean < 12.){
    x=-1.;
    double exp_nxm=exp(-xmean);
    double uni_var_product=1.;
    while(uni_var_product > exp_nxm){
      x=x+1.;
      double rand = (double) curand_uniform(mystate);
      uni_var_product=uni_var_product*rand;
    }
  }else{
   double sq = sqrt(2.0*xmean);
   double log_xmean = log(xmean);
   double GG = xmean*log_xmean - LnGamma(xmean+1.0);
   do {
     do {
       double rand = (double) curand_uniform(mystate);
       xcomp = tan(pi*rand);
       x = sqrt(2.*xmean)*xcomp + xmean;
     } while (x < 0.0);
     x = floor(x);
     reject = reject_factor*(1.0 + xcomp*xcomp)*exp(x*log_xmean - LnGamma(x+1.0) - GG);
   } while (curand_uniform(mystate) > reject);
  }
  return x;
}
