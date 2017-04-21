#ifndef FILE_RigidBody_SEEN
#define FILE_RigidBody_SEEN
/*---------------------------------------*/
#include<math.h>
#include"3vec.h"
#include"quaternion.h"
#include"matrix.h"
/*---------------------------------------*/
unsigned int const ddim=3;
unsigned int const pdim=4*ddim+1;
double const pi = 4.*atan(1.);
// Properties of the rigid body
/* Input parameters */
double const ap = 1.;
double const lambda = 6.;
double const rhop = 10.;
double const mu = 0.01;
/* Derived Parameters */
double const Mass= ( 4.*pi*rhop* pow(ap,3)  )/(3.*lambda*lambda);
double const taup0 =  (2.*rhop*ap*ap)/(9.*mu*lambda*lambda);
double const I0 = ( 4*pi*rhop*pow(ap,5) )/(15.*lambda*lambda) ;
double const beta=0.;
// Moment of inertial matrix and its inverse 
// PLEASE MAKE SURE THAT Ibody is inverse of Ibodyinv
//double const Ibody[3][3]={{1,0,0},{0,1,0},{0,0,1}}; // sphere
//double const Ibodyinv[3][3]={{1,0,0},{0,1,0},{0,0,1}}; //sphere
//double const Ibody[3][3]={{0.1,0,0},{0,0.1,0},{0,0,0.2}}; // disk in x-y plane
//double const Ibodyinv[3][3]={{10.,0,0},{0,10.,0},{0,0,5.}}; // disk in x-y plane
double const Iperp=I0*( 1. + 1./(lambda*lambda) ) ;
double const Iparl= I0/(lambda*lambda);
double const Ixx = Iperp;
double const Iyy = Iperp;
double const Izz = Iparl;
double const Ibody[3][3]={{Ixx,0,0},{0,Iyy,0},{0,0,Izz}}; // ellipsoid
double const Ibodyinv[3][3]={{1./Ixx,0,0},{0,1./Iyy,0},{0,0,1./Izz}}; // ellipsoid
// The Drag tensor in principal axis frame 
double const tau = lambda*lambda -1.;
double const Kperp=(16.)/(6.*lambda)*( pow(tau,3)  )/(  (2.*tau*tau-1.)*log(lambda+tau) + lambda*tau  );
double const Kparl= (8.)/(6.*lambda)* ( pow(tau,3) )/( (2.*tau*tau + 1.)*log(lambda+tau) -lambda*tau );
double const KK[3][3]={{Kperp,0,0},{0.,Kperp,0.},{0.,0.,Kparl}}; // spheroid
double const Jeff_Rone = ( pow(ap,3)/pow(lambda,3) )*(   mu*16.*pi*pow(tau,3)  )/(3.*(lambda*tau - log(lambda+tau) )  );
double const Jeff_Rtwo = ( pow(ap,3)/pow(lambda,3) )*(  mu*16.*pi*pow(tau,3)*(1.+lambda*lambda) )/(3.*( (2.*tau*tau+1.)*log(lambda+tau) -lambda*tau ) );
double const Jeff_D = (lambda*lambda-1.)/(lambda*lambda+1.);
/* ----------------------------------------*/
struct RigidBody{
  /* State variables */
  vec3 xx; // position
  quaternion qq; // quaternion describing its rotated state 
  vec3 pp, Ell ; // momentum  and angular momentum
  /* Derived quantities */
  double RR[3][3]; // Rotation matrix corresponding to the quaternion;
  double Iinv[3][3]; // inverse of moment of inertial matrix in space fixed frame 
  vec3 omega; // angular velocity
};
/* ----------------------------------------*/
void wparam_rigid_body(){
  cout<< "----------*********----------------\n";
  cout<< "------Parameters of Simulation--------------\n"; 
  cout<<"Major Axis ap="<<ap<<"\n";
  cout<<"Prolate Ellipsoid  lambda( > 1) ="<<lambda<<"\n";
  cout<<"Particle density="<<rhop<<"\n";
  cout<<"Dynamic viscosity mu="<<mu<<"\n";
/* Derived Parameters */
  cout<<"Mass="<<Mass<<"\n";
  cout<<"taup0="<<taup0<<"\n";
  cout<<"I0="<<I0<<"\n";
  cout<<"beta(not used in equations)="<<beta<<"\n";
  cout<<"Iperp="<<Iperp<<"\n";
  cout<<"Iparl="<<Iparl<<"\n";
  cout<<"-- Moment-of-Inertia matrix--"<<"\n";
  PMat(Ibody);
  cout<<"-- Inverse of Moment-of-Inertia matrix--"<<"\n";
  PMat(Ibodyinv);
  cout<<"tau="<<tau<<"\n";
  cout<<"Kperp="<<Kperp<<"\n";
  cout<<"Kparl="<<Kparl<<"\n";
  cout<<"-- Drag matrix--"<<"\n";
  PMat(KK);
  cout<<"----Jeffery's Torque--"<<"\n";
  cout<<"Jeff_Rone="<<Jeff_Rone<<"\n";
  cout<<"Jeff_Rtwo="<<Jeff_Rtwo<<"\n";
  cout<<"Jeff_D="<<Jeff_D<<"\n";
  cout<< "-------------------------------\n"; 
}
/* ----------------------------------------*/
void array2rb(RigidBody *rb,  double y[pdim]){
  double  temp[3][3], Ibinv[3][3];
  /* 0,1,2 are position */
  rb->xx.x=y[0];
  rb->xx.y=y[1];
  rb->xx.z=y[2];
  /* 3,4,5 are momentum */
  rb->pp.x=y[3];
  rb->pp.y=y[4];
  rb->pp.z=y[5];
  /* 6,7,8,9 are quaternion */
  rb->qq.w=y[6];
  rb->qq.u.x=y[7];
  rb->qq.u.y=y[8];
  rb->qq.u.z=y[9];
  /* now normalize the quaternion */
  double normq=norm(rb->qq);
  if (normq!=0){
    rb->qq = rb->qq*(1./normq);
  }
  /* 10, 11, 12 are the angular momentum */
  rb->Ell.x=y[10];
  rb->Ell.y=y[11];
  rb->Ell.z=y[12];
  /* Calculate the rotation matrix equivalent to qq */
  Rot_from_Q(rb->RR, rb->qq);
  // This rotation matrix is the transpose of the matrix A used
  // in the paper, Phys. of Fluids. 20 093302 (2008) by 
  //Mortensent et al (includes Boersma) 
  /* The inverse of moment-of-inertia in body fixed frame */
  for(int j=0;j<3;j++){
    for(int k=0;k<3;k++){
      Ibinv[j][k]=Ibodyinv[j][k];
    }
  }
  RMRt(rb->Iinv,Ibinv,rb->RR);
  rb->omega=MatVec(rb->Iinv, rb->Ell);
  //cout << "omega and ell\n";
  //PVec3(rb->omega);
  //PVec3(rb->Ell);
  //cout << dot(rb->omega,rb->Ell)<<"\n";
  //PVec3(rb->omega);
}
/* ----------------------------------------*/
void rb2array(double y[pdim], RigidBody *rb){
  /* 0,1,2 are position */
  y[0]=rb->xx.x;
  y[1]=rb->xx.y;
  y[2]=rb->xx.z;
  /* 3,4,5 are momentum */
  y[3]=rb->pp.x;
  y[4]=rb->pp.y;
  y[5]=rb->pp.z;
  /* 6,7,8,9 are quaternion */
  y[6]=rb->qq.w;
  y[7]=rb->qq.u.x;
  y[8]=rb->qq.u.y;
  y[9]=rb->qq.u.z;
  /* 10, 11, 12 are the angular momentum */
  y[10]=rb->Ell.x;
  y[11]=rb->Ell.y;
  y[12]=rb->Ell.z;
}
/* ----------------------------------------*/
#endif /* !FILE_RigidBody_SEEN */
