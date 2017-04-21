#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "ode.h"
#include "input.h"
#include "model.h"
#include "CUDA.h"
#include "Random.h"
#include "SpecialFunctions.h"
using namespace std;
/**************************/
/* ----------------------------------------*/
int main(){
  //--------------------------------------//
  unsigned int const ndim=Nensemble*pdim;
  curandState rand_state[Nensemble];
  curandState *dev_rand_state;
  double y[ndim];
/* ---- to debug the noise --------*/
  double noise[Nensemble];
  double *dev_noise;
  for(int i=0;i<Nensemble;i++){
    noise[i]=0.;
  }
 dev_noise=host2dev(Nensemble,noise); 
/*----------------------------------*/
  double time=0.;
  double tnext;
  /* device variables */
  double *dev_y;
  /* end device variables */
  iniconf(y, Nensemble, rand_state);
/*---------------------------------*/
  dev_y=host2dev(ndim,y);
  dev_rand_state=host2dev(Nensemble,rand_state);
  double dtnext=TMAX/Ndiag;
  tnext=dtnext;
  /*----------------------------------------------------------------*/
  FILE *tseries=fopen("tseries.out","w");
  FILE *diagf=fopen("diag.dat","w");
  while (time < TMAX){
    diag(time,y,Nensemble,tseries,diagf);
    evolve<<<Nensemble,1>>>(time,dt,tnext,dev_y,dev_rand_state);
    dev2host(y,ndim,dev_y); 
    dev2host(rand_state,Nensemble,dev_rand_state);
    dev2host(noise,Nensemble,dev_noise);
    time=time+dtnext;
    tnext=time+dtnext;
    printf("Next diagnostic at %f\n",tnext);
  }
  fclose(tseries);
  fclose(diagf);
  FILE *fin_state=fopen("yfin.out","w");
  fprintf(fin_state,"%lf\n",time);
  for (int i=0;i<ndim;i++){
    fprintf(fin_state,"%lf\n",y[i]);
  }
  fclose(fin_state);
//----------------------------
}
/* ----------------------------------------*/
