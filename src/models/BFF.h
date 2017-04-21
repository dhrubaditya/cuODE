#ifndef FILE_BFF_SEEN
#define FILE_BFF_SEEN
#include<stdio.h>
/*---------------------------------------*/
#define pdim 2
extern __device__ void eval_rhs(double rhs[],double tt,double yy[],int lindex);
extern __device__ void stochastic(double yy[],curandState global_state[], double tlocal,
         double deltat,int lindex);
void iniconf(double y[],int Nensemble, curandState rand_state[]);
extern __host__ void diag(double tt, double y[], int Nensemble, FILE* tseries, FILE* diagf);
/* ----------------------------------------*/
#endif /* !FILE_BFF_SEEN */
