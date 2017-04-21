#include<math.h>
#include<stdio.h>
#include "CUDA.h"
#include "Random.h"
#include "mycomplex.h"
#include "model.h"
using namespace std;
/* ----------------------------------------*/
/**************************/
  double *yzero;
__device__ double const tau=1.;
__device__ double const one_over_tau=1./tau;
double const iniy_max=10.;
double const iniy_min=-10.;
__device__ inline double kappa(double zz){
  double amp=0.0387;
  double kzero=0.0017;
/*  double zmax=100.;
  double mu=0.1;
  double zabs=abs(zz);
  if (zabs > zmax){
  return mu*zmax*zmax;
  }
  else{
  return mu*zz*zz;
  } */
  return amp*sin(zz)*sin(zz)+kzero;
}
/* ----------------------------------------*/
__global__ void inirand_evolve(unsigned long long seed[]);
/* ----------------------------------------*/
/*__device__ double telegraph(double nu, double tt, int local_index, curandState mystate){
  double poisson_mean=nu*tt;
  double pr=Poisson(poisson_mean,&mystate);
  int ppower=(int) fmod(pr,2.);
  double tele_ran = powf(-1,ppower);
  return tele_ran;
}*/
/* ----------------------------------------*/
__device__ void eval_rhs(double rhs[],double tt,double yy[],int lindex){
  double vv=yy[lindex+1];
  double zdot = vv;
/* --The stochastic part of the equation is added outside the usual integrator - */
  double vdot = -one_over_tau*vv;
/* ---------------------------------------------------- */
  rhs[0]=zdot;
  rhs[1]=vdot;
}
/* ----------------------------------------*/
__device__ void stochastic(double yy[],curandState global_state[], double tlocal,
         double deltat,int lindex)
{
  double pi = 4.*atan(1.);
  double zz=yy[lindex];
  //double r = fmod(zz,pi);
  double mean=0;
  //double sigma=1.;
  double sigma=kappa(zz);
  int tid=lindex/pdim;
  curandState local_state=global_state[tid];
  double uu = Gaussian(mean,sigma,&local_state); 
  global_state[tid] = local_state;
  yy[lindex+1]=yy[lindex+1]+one_over_tau*uu*sqrt(deltat);
}
/*---------------
__global__ void inirand_evolve(unsigned long long seed[], dev_global_state[]){
  int tid = blockIdx.x;
  unsigned long long local_seed = seed[tid];
  curandState local_state;
  local_state = dev_global_state[tid];
  curand_init(local_seed,tid,0, &local_state);
  dev_global_state[tid] = local_state;
}*/
/* ----------------------------------------*/
void iniconf(double y[],int Nensemble, curandState rand_state[]){
  curandState *dev_iniran_state;
  double rand[Nensemble],rand2[Nensemble];
  double *dev_rand;
  unsigned long long seed[Nensemble];
  unsigned long long *dev_seed;
  for(int i=0;i<Nensemble;i++){
    seed[i]=37*i+53*i*i;
    rand[i]=0.;
    rand2[i]=0.;
  }
  dev_rand= host2dev(Nensemble,rand);
  dev_seed =  host2dev(Nensemble,seed);
  cudaMalloc( (void**)&dev_iniran_state, Nensemble*sizeof(curandState) );
  init_random<<<Nensemble,1>>>(dev_seed,dev_iniran_state);
  UniformRandom<<<Nensemble,1>>>(dev_rand, dev_iniran_state);
  dev2host(rand,Nensemble,dev_rand);
  UniformRandom<<<Nensemble,1>>>(dev_rand, dev_iniran_state);
  dev2host(rand2,Nensemble,dev_rand);
  for(int j=0;j<Nensemble;j++){
// Uniformly distributed initial position between iniy_min to iniy_max
      y[0+j*pdim]=iniy_min+rand[j]*(iniy_max-iniy_min);
// and random initial velocity
      y[1+j*pdim]=rand2[j];
      printf("y0,y1,%lf,%lf\n",y[0],y[1]);
  }
  /* copy the state of the random no. generator to host */
  dev2host(rand_state,Nensemble,dev_iniran_state);
  //  inirand_evolve<<<Nensemble,1>>>(dev_seed, dev_rand_state);
}
/* ----------------------------------------*/
__host__ void diag(double tt, double y[], int Nensemble, FILE* tseries, FILE* diagf){
  int ndim=pdim*Nensemble;
  if (tt == 0.) {
     yzero=(double*)malloc(ndim*sizeof(double));
     for (int i=0;i<ndim;i++){
       yzero[i]=y[i];
     }
  }
  //printf("%lf\t%lf\t%lf\t%lf\t%lf\n",tt,y[0],y[1],y[2],y[3]);
  fprintf(tseries,"%lf\t",tt);
  for (int i=0;i<ndim-1;i++){
    fprintf(tseries,"%lf\t",y[i]);
  }
  fprintf(tseries,"%lf\n",y[ndim-1]);
  double meanz=0.;
  double meanv=0.;
  double dzrms=0;
  for(int i=0; i<Nensemble; i++){
    int lindex=pdim*i;
    double zz=y[lindex];
    double dz=y[lindex]-yzero[lindex];
    double vv=y[lindex+1];
    meanz= zz+meanz ;
    dzrms= dz*dz+dzrms ;
    meanv= vv+meanv ;
  }
  meanz=meanz/Nensemble;
  meanv=meanv/Nensemble;
  dzrms=sqrt(dzrms)/Nensemble;
  printf("%lf\t%lf\t%lf\t%lf\n",tt,dzrms,meanz,meanv);
  fprintf(diagf,"%lf\t%lf\t%lf\t%lf\n",tt,dzrms,meanz,meanv);
}
/* ----------------------------------------*/
