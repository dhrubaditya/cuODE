#include <math.h>
#include "ode.h"
#include "model.h"
#include "Random.h"
using namespace std;
/* ------------------------ */
__global__ void evolve(double time, double dt, double tnext, double yy[], curandState dev_rand_state[]){
  int tid = blockIdx.x;
  int local_index= tid*pdim;
  double tlocal=time;
  double deltat=dt;
  while(tlocal < tnext){
    stochastic(yy,dev_rand_state,tlocal,deltat,local_index);
    //euler(yy,  tlocal, deltat, local_index); 
    //rnkt2(yy,  tlocal, deltat, local_index);
    rnkt4(yy, tlocal, deltat, local_index);
    tlocal=tlocal+deltat;
  }
}
/*********************************/
__device__ void euler(double yy[], double tt,double deltat, int lindex){
  double k1[pdim];
  eval_rhs(k1,tt,yy,lindex);
  for(int idim=0;idim<pdim;idim++){
    yy[lindex+idim]=yy[lindex+idim]+deltat*k1[idim];
  }
}
/*********************************/
__device__ void rnkt2(double yy[], double tt,double deltat, int lindex){
  double temp[pdim],k1[pdim];
  eval_rhs(k1,tt,yy,lindex);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k1[idim]*deltat/2.;
  }
  eval_rhs(k1,tt+(deltat/2.),temp,0);
  for(int idim=0;idim<pdim;idim++){
    yy[lindex+idim]=yy[lindex+idim]+deltat*k1[idim];
  }
}
/*********************************/
__device__ void rnkt4(double yy[], double tt,double deltat, int lindex){
  double  temp[pdim],k1[pdim],k2[pdim],k3[pdim],k4[pdim];
  eval_rhs(k1,tt,yy,lindex);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k1[idim]*deltat/2.;
  }
  eval_rhs(k2,tt+(deltat/2.),temp,0);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k2[idim]*deltat/2.;
  }
  eval_rhs(k3,tt+(deltat/2.),temp,0);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k3[idim]*deltat;
  }
  eval_rhs(k4,tt+deltat,temp,0);
  for(int idim=0;idim<pdim;idim++){
    yy[lindex+idim]=yy[lindex+idim]+deltat*(  (k1[idim]/6.) + (k2[idim]/3.) + (k3[idim]/3.) + (k4[idim]/6.) );
  }
}
/*********************************/
