#include <iostream>
#include<fstream>
#include "CUDA.h"
#include "Random.h"
#include "modules/quaternion.h"
#include "modules/matrix.h"
#include "modules/RigidBody.h"
#include "model.h"
using namespace std;
/**************************/
__device__ double const omega2=1.;
__device__ unsigned int const ndim=pdim*Nensemble;
vec3 inipos(int posini);
__device__ vec3 inivel(int velini);
__device__ quaternion iniq(int qini);
vec3 iniell(int ellini, quaternion qq_zero);
vec3 fluid_velocity(vec3 );
__device__ void U_Omega_GradU(vec3& UU, vec3& OO, double GradU[3][3], vec3 xx);
__device__ vec3 DragForce(vec3 UU,RigidBody *rb);
__device__ vec3 JefferyTorque(vec3 OO,double Sij[3][3],RigidBody *rb);
__host__ void wparam_sim(double time);
/* ----------------------------------------*/
__device__ void stochastic(double yy[],curandState dev_rand_state[], double tlocal,
         double deltat,int lindex){}
/* ----------------------------------------*/
__device__ void eval_rhs(double rhs[],double tt,double yy[],int lindex){
  /* we solve:
        Equations of a rigid body advected by known flow
 */
  double xx=yy[lindex];
  double vv=yy[lindex+1];


  RigidBody rb;
  double GradU[3][3],Sij[3][3];
  vec3 UU,OO,vv,Force,Torque;
  quaternion omegaq,dq_dt;
  //
  for(int ibody=0;ibody<Nensemble;ibody++){
    int irb=pdim*ibody;
    array2rb(&rb, &y[irb]);
    U_Omega_GradU(UU, OO, GradU, rb.xx);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Sij[i][j]= 0.5*(GradU[i][j]+GradU[j][i]);
      }
    }
    /* vec3 UXO = cross(UU,OO);
    cout<<"U, O, UXO"<<"\n";g
    PVec3(UU);
    PVec3g(OO);
    PVec3(UXO);
    cout<<"-----------------------"<<"\n"; */
//----- Calculate velocity -------
    vv = rb.pp*(1./Mass);
//---Calculate Force---------
    Force=DragForce(UU,&rb);
    //PVec3(Force);
// -- Now for zero torque -----------
    Torque=JefferyTorque(OO,Sij,&rb);
/* ---- Calculate the evolution eqn for the quaternion 
   dq_dt = (1/2) omega * q  */
    omegaq = quaternion(0.,rb.omega);
    dq_dt = omegaq * rb.qq *(1./2.);

// ----------------------------------
    rhs[irb+0]=rb.pp.x/Mass;
    rhs[irb+1]=rb.pp.y/Mass;
    rhs[irb+2]=rb.pp.z/Mass;
    //
    rhs[irb+3]=Force.x;
    rhs[irb+4]=Force.y;
    rhs[irb+5]=Force.z;
    //
    rhs[irb+6]=dq_dt.w;
    rhs[irb+7]=dq_dt.u.x;
    rhs[irb+8]=dq_dt.u.y;
    rhs[irb+9]=dq_dt.u.z;
    //
    rhs[irb+10]=Torque.x;
    rhs[irb+11]=Torque.y;
    rhs[irb+12]=Torque.z;
    }







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




