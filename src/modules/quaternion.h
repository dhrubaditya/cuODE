#ifndef FILE_quaternion_SEEN
#define FILE_quaternion_SEEN
/*---------------------------------------*/
#include <iostream>
#include <math.h>
#include "3vec.h"
using namespace std;

class quaternion{
public:
  double w;
  vec3 u;
  quaternion();
  quaternion(int,int,int,int);
  quaternion(float,float,float,float);
  quaternion(double,double,double,double);
  quaternion(double, vec3);
  // definining operators 
  quaternion operator*(quaternion);
  quaternion operator*(double);
  quaternion operator*(vec3);
  quaternion operator+(quaternion);
  quaternion operator-(quaternion);
  quaternion operator/(quaternion){};
  #define qzero quaternion(0.0,0.0,0.0,0.0);
};

quaternion::quaternion(){
  w = 0.0;
  u.x=0.;
  u.y=0.;
  u.z=0.;
}

quaternion::quaternion(int a, int b, int c, int d){
  w = (double) a;
  u.x = (double) b;
  u.y = (double) c;
  u.z = (double) d;
}

quaternion::quaternion(float a, float b, float c, float d){
  w = (double) a;
  u.x = (double) b;
  u.y = (double) c;
  u.z = (double) d;
}

quaternion::quaternion(double a, double b, double c, double d){
  w =  a;
  u.x=b;
  u.y=c;
  u.z=d;
}

quaternion::quaternion(double a, vec3 v){
  w =  a;
  u=v;
}

quaternion quaternion::operator+(quaternion param){
  quaternion temp;
  temp.w = w+param.w;
  temp.u = u+param.u;
  return(temp);
}

quaternion quaternion::operator-(quaternion param){
  quaternion temp;
  temp.w = w-param.w;
  temp.u = u-param.u;
  return(temp);
}

quaternion quaternion::operator*(double param){
  quaternion temp;
  temp.w = w*param;
  temp.u = u*param;
  return(temp);
}

quaternion quaternion::operator*(quaternion param){
  quaternion temp;

  temp.w = param.w*w - dot(u,param.u);
  temp.u = param.u*w + u*param.w + cross(u,param.u);
  return(temp);
}

quaternion quaternion::operator*(vec3 param){
  quaternion temp;
  temp=quaternion(w,u)*quaternion(0,param);
  return(temp);
}

quaternion conjg(quaternion a){
  quaternion temp;
  temp.w = a.w;
  temp.u = vec3(0.0,0.0,0.0)-a.u;
  return(temp);
}

double norm(quaternion a){
  return( sqrt( a.w*a.w+dot(a.u,a.u)  ) );}
/*---------------------------------------*/
void PQ(quaternion a){
  cout<<a.w<<"\t"<<a.u.x<<"\t"<<a.u.y<<"\t"<<a.u.z<<"\n";
}
/*---------------------------------------*/
void Rot_from_Q(double RR[3][3], quaternion q){
  RR[0][0] = 1.- 2.* ( pow(q.u.y,2)  + pow(q.u.z,2) ) ;
  RR[0][1] = 2.* (q.u.x*q.u.y - q.w*q.u.z)   ;
  RR[0][2] = 2.* (q.u.x*q.u.z + q.w*q.u.y)   ;
  RR[1][0] = 2.* (q.u.x*q.u.y + q.w*q.u.z)   ;
  RR[1][1] = 1.- 2.* ( pow(q.u.x,2)  + pow(q.u.z,2) ) ;
  RR[1][2] = 2.* (q.u.y*q.u.z - q.w*q.u.x)   ;
  RR[2][0] = 2.* (q.u.x*q.u.z - q.w*q.u.y)   ;
  RR[2][1] = 2.* (q.u.y*q.u.z + q.w*q.u.x)   ;
  RR[2][2] = 1.- 2.* ( pow(q.u.x,2)  + pow(q.u.y,2) ) ;
}
/*---------------------------------------*/
#endif /* !FILE_quaternion_SEEN */
