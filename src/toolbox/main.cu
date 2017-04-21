#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "CUDA.h"
#include "SpecialFunctions.h"
#include "Random.h"
#define pen 2
/* ------------------------ */
__global__ void evolve(double *a,double *pr,
                       double *tt,double *deltat, 
                       double *tnxt,double *nuu,
                       double *a0, 
                       curandState *global_state){
    int tid = blockIdx.x;
    curandState mystate=global_state[tid];
    if (tid < pen) {
      while(*tt < *tnxt){
        double poisson_mean=(*nuu)*(*tt);
        pr[tid]=Poisson(poisson_mean,&mystate);
        int ppower=(int) fmod(pr[tid],2.);
        a[tid] = (*a0)*powf(-1,ppower);
        *tt=*tt+*deltat;
      }
    }
    global_state[tid] = mystate;
}
/* ------------------------ */
int main(void){
  double alpha[pen],prand[pen];
  double *dev_alpha,*dev_prand;
  unsigned long long seed[pen];
  unsigned long long *dev_seed;
  curandState *dev_global_state;
  double dt=1e-4;
  double *dev_dt;
  double TMAX=4.;
  double tnext;
  double *dev_tnext;
  double t=0.;
  double *dev_t;
  double alpha0=1.;
  double *dev_alpha0;
  double tauc=0.01;
  double nu=1./tauc;
  double *dev_nu;
  int Ndiag=10;
/* ------- set host values ----------- */
  seed[0]=37;
  seed[1]=53;
  dev_seed =  host2dev(pen,seed);
  cudaMalloc( (void**)&dev_global_state, pen*sizeof(curandState) );
  init_random<<<pen,1>>>(dev_seed,dev_global_state);
  alpha[0]=1.;
  alpha[1]=-1.;
  dev_alpha = host2dev(pen,alpha);
  prand[0]=0.;
  prand[1]=1.;
  dev_prand = host2dev(pen,prand);
/*---------------------------------*/
  dev_t = host2dev(1,&t);
  dev_dt = host2dev(1,&dt);
  double dtnext=TMAX/Ndiag;
  tnext=dtnext;
  dev_tnext = host2dev(1,&tnext);
  dev_nu = host2dev(1,&nu);
  dev_alpha0 = host2dev(1,&alpha0);
/*----------------------------------------------------------------*/
  while (t < TMAX){
    evolve<<<pen,1>>>(dev_alpha,dev_prand,dev_t,dev_dt,dev_tnext,dev_nu,dev_alpha0,dev_global_state);
   dev2host(alpha,pen,dev_alpha); 
   dev2host(prand,pen,dev_prand); 
   t=t+dtnext;
   h2d(dev_t,1,&t);
   tnext=t+dtnext;
   h2d(dev_tnext,1,&tnext);
   printf("%lf\t%lf\t%lf\n",t,prand[0],alpha[0]);
   printf("Next diagnostic at %f\n",tnext);
 }
  
}
