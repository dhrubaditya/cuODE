#ifndef FILE_ODE_SEEN
#define FILE_ODE_SEEN
#include <curand.h>
#include <curand_kernel.h>
/*------------------------------------*/
//extern __global__ void evolve(double time, double dt, double tnext, double yy[], curandState dev_rand_state[]);
extern __global__ void evolve(double time, double dt, double tnext, double yy[], curandState dev_rand_state[] );
extern __device__ void euler(double yy[], double tt, double deltat, int lindex);
extern __device__ void rnkt2(double yy[], double tt, double deltat, int lindex);
extern __device__ void rnkt4(double yy[], double tt,double deltat, int lindex);
/*------------------------------------*/
#endif /* !ODE_SEEN */
