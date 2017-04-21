#include<math.h>
#include "CUDA.h"
#include "Random.h"
#include "mycomplex.h"
#include "model.h"
using namespace std;
/* ----------------------------------------*/
/**************************/
__device__ double const kk=1.;
__device__ double ksqr=kk*kk;
__device__ double const alpha0=1.;
__device__ double const tauc=0.1;
__device__ double nu = 1./tauc;
__device__ double const Shear=-4.;
__device__ double const eta=1.;
/* ------------------------------ */
double Biamp = 1.e-5;
struct Bfield{
  /* State variables */
  complex x; // x component of magnetic field
  complex y; // y component of magnetic field
};
/* ----------------------------------------*/
__global__ void inirand_evolve(unsigned long long seed[]);
__host__ __device__ void array2b(Bfield *B,  double y[]);
__host__ __device__ void b2array(double y[], Bfield *B);
/* ----------------------------------------*/
__device__ double telegraph(double nu, double tt, int local_index, curandState mystate){
  double poisson_mean=nu*tt;
  double pr=Poisson(poisson_mean,&mystate);
  int ppower=(int) fmod(pr,2.);
  double tele_ran = powf(-1,ppower);
  return tele_ran;
}
/* ----------------------------------------*/
__device__ void eval_rhs(double rhs[],double tt,double yy[],int lindex){
  Bfield B,dtB;
  complex I=complex(0.,1.);
//
  array2b(&B, &yy[lindex]);
/* --The stochastic part of the equation is added outside the usual integrator - */
  dtB.x= complex(0.,0.) - B.x*(eta*ksqr);
  dtB.y= B.x*Shear -B.y*(eta*ksqr);
// ----------------------------------
  b2array(&rhs[0],&dtB);
//
}
/* ----------------------------------------*/
__device__ void stochastic(double yy[],curandState rstate, double tlocal,
         double deltat,int lindex)
{
  Bfield B,dtB;
  complex I=complex(0.,1.);
  array2b(&B, &yy[lindex]);
  double alpha=alpha0*telegraph(nu,tlocal,lindex,rstate);
  B.x= B.x+I*B.y*(kk*alpha)*sqrt(deltat);
  /* B.x=complex(alpha,0.);
  B.y=complex(alpha,0.); */
  b2array(&yy[lindex],&B);
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
  double rand[Nensemble];
  double *dev_rand;
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
  for(int i=0;i<pdim;i++){
    UniformRandom<<<Nensemble,1>>>(dev_rand, dev_iniran_state);
    dev2host(rand,Nensemble,dev_rand);
    for(int j=0;j<Nensemble;j++){
      y[i+j*pdim]=rand[j]*Biamp;
    }
  }
  /* copy the state of the random no. generator to host */
  dev2host(rand_state,Nensemble,dev_iniran_state);
  //  inirand_evolve<<<Nensemble,1>>>(dev_seed, dev_rand_state);
}
/* ----------------------------------------*/
__host__ void diag(double tt, double y[], int Nensemble, FILE* tseries, FILE* diagf){
  Bfield B;
  int ndim=pdim*Nensemble;
  //printf("%lf\t%lf\t%lf\t%lf\t%lf\n",tt,y[0],y[1],y[2],y[3]);
  fprintf(tseries,"%lf\t",tt);
  for (int i=0;i<ndim-1;i++){
    fprintf(tseries,"%lf\t",y[i]);
  }
  fprintf(tseries,"%lf\n",y[ndim-1]);
  fprintf(diagf,"%lf\t",tt);
  double meanBxreal=0.;
  for(int i=0; i<Nensemble; i++){
    int lindex=pdim*(Nensemble-1);
    array2b(&B, &y[lindex]);
    meanBxreal= B.x.real+meanBxreal ;
    printf("%lf\t%lf\n",tt,meanBxreal);
    fprintf(diagf,"%lf\t",meanBxreal);
  }
  fprintf(diagf,"\n");
}
/* ----------------------------------------*/
__host__ __device__ void array2b(Bfield *B,  double y[]){
  /* real and imaginary part of Bx */
  B->x=complex(y[0],y[1]);
  /* real and imaginary part of By */
  B->y=complex(y[2],y[3]);
}
/* ----------------------------------------*/
__host__ __device__ void b2array(double y[], Bfield *B){
  y[0]=B->x.real;
  y[1]=B->x.imag;
  y[2]=B->y.real;
  y[3]=B->y.imag;
}

