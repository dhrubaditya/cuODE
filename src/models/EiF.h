#ifndef FILE_EiF_SEEN
#define FILE_EiF_SEEN
/*---------------------------------------*/
#include <iostream>
#include<fstream>
#define pdim 2
extern __device__ void eval_rhs(double rhs[],double tt,double yy[],int lindex);
extern __device__ void stochastic(double yy[],curandState rstate[], double tlocal,
         double deltat,int lindex);
extern __host__ void iniconf(double y[],int Nensemble, curandState rand_state[]);
extern __host__ void diag(double tt, double y[], int Nensemble, FILE* tseries, FILE* diagf);
/* ----------------------------------------*/
#endif /* !FILE_EiF_SEEN */
/* ----------------------------------------*/
