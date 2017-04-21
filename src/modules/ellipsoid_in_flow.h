/*---------------------------------------*/
#include"RigidBody.h"

#include <iostream>
#include<fstream>
#define pdim 2
extern __device__ void eval_rhs(double rhs[],double tt,double yy[],int lindex);
extern __device__ void stochastic(double yy[],curandState rstate[], double tlocal,
         double deltat,int lindex);
extern __host__ void iniconf(double y[],int Nensemble, curandState rand_state[]);
extern __host__ void diag(double tt, double y[], int Nensemble, FILE* tseries, FILE* diagf);




unsigned int const Nensemble=1;
// Properties of Flow
double const kk=0.1;
double const rho = 1.;
double const uzero=0.01;
int Utype = 1; //ABC flow
// Initial rotation of the body about the space-fixed axis 
int q_ini = 2; // rotated about x axis by psi_zero
double psi_zero=pi/4.;
// Initial angular momentum 
int ell_ini = 2; // angular *velocity* along z axis
double omega_zero=0.1; 
// Time step
double const TMAX=1.;
double dt=0.01;
int idiag=10;
/*---------------------------------------*/
