#include <iostream>
#include<fstream>
#include "CUDA.h"
#include "Random.h"
#include "model.h"
using namespace std;
/**************************/
__device__ double const omega2=1.;
/* ----------------------------------------*/
__device__ void stochastic(double yy[],curandState dev_rand_state[], double tlocal,
         double deltat,int lindex){}
/* ----------------------------------------*/
__device__ void eval_rhs(double rhs[],double tt,double yy[],int lindex){
  /* we solve:
        (d/dt)x = v;  (d/dt)v = -\omega^2 x
 */
  double xx=yy[lindex];
  double vv=yy[lindex+1];
  rhs[0]=vv;
  rhs[1]=-omega2*xx;
}
/* ----------------------------------------*/
__host__ void iniconf(double y[],int Nensemble, curandState rand_state[]){
  double rand[Nensemble];
  double *dev_rand;
  curandState *dev_iniran_state;
  unsigned long long seed[Nensemble];
  unsigned long long *dev_seed;
  for(int i=0;i<Nensemble;i++){
    seed[i]=37*i+53*i*i;
    rand[i]=0.;
  }
  dev_rand= host2dev(Nensemble,rand);
  dev_seed =  host2dev(Nensemble,seed);
  cudaMalloc( (void**)&dev_iniran_state, Nensemble*sizeof(curandState) );
  init_random<<<Nensemble,1>>>(dev_seed,dev_iniran_state);
  UniformRandom<<<Nensemble,1>>>(dev_rand, dev_iniran_state);
  dev2host(rand,Nensemble,dev_rand);
  for(int j=0;j<Nensemble;j++){
    y[0+j*pdim]=rand[j];
  }
  UniformRandom<<<Nensemble,1>>>(dev_rand, dev_iniran_state);
  dev2host(rand,Nensemble,dev_rand);
  for(int j=0;j<Nensemble;j++){
    y[1+j*pdim]=rand[j];
  }
  dev2host(rand_state,Nensemble,dev_iniran_state);
}
/* ----------------------------------------*/
__host__ void diag(double tt, double y[], int Nensemble, FILE* tseries, FILE* diagf){
  int ndim=pdim*Nensemble;
  printf("%lf\t%lf\t%lf\t%lf\t%lf\n",tt,y[0],y[1],y[2],y[3]);
  fprintf(tseries,"%lf\t",tt);
  for (int i=0;i<ndim-1;i++){
    fprintf(tseries,"%lf\t",y[i]);
  }
  fprintf(tseries,"%lf\n",y[ndim-1]);
  fprintf(diagf,"%lf\t",tt);
  for(int i=0; i<Nensemble; i++){
    double xx=y[i*pdim+0];
    double vv=y[i*pdim+1];
    double energy=omega2*xx*xx+vv*vv;
    fprintf(diagf,"%lf\t",energy);
  }
  fprintf(diagf,"\n");
}
/* ----------------------------------------*/




