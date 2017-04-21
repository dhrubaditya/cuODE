#ifndef FILE_mycomplex_SEEN
#define FILE_mycomplex_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
using namespace std;
/*-------------------------*/
class complex{
public:
  double real,imag;
  __host__ __device__ __inline__ complex();
  __host__ __device__ __inline__ complex(int,int);
  __host__ __device__ __inline__ complex(float,float);
  __host__ __device__ __inline__ complex(double,double);
  // definining operators 
  __host__ __device__ __inline__ complex operator+(complex);
  __host__ __device__ __inline__ complex operator-(complex);
  __host__ __device__ __inline__ complex operator*(complex);
  __host__ __device__ __inline__ complex operator*(double);
};

__host__ __device__ __inline__ complex::complex(){
  real = 0.0;
  imag = 0.0;
}
__host__ __device__ __inline__ complex::complex(int a, int b){
  real = (double) a;
  imag = (double) b;
}
__host__ __device__ __inline__ complex::complex(float a, float b){
  real = (double) a;
  imag = (double) b;
}
__host__ __device__ __inline__ complex::complex(double a, double b){
  real =  a;
  imag =  b;
}
__host__ __device__ __inline__ complex complex::operator+(complex param){
  complex temp;
  temp.real = real+param.real;
  temp.imag = imag+param.imag;
  return(temp);
}
__host__ __device__ __inline__ complex complex::operator-(complex param){
  complex temp;
  temp.real = real-param.real;
  temp.imag = imag-param.imag;
  return(temp);
}
__host__ __device__ __inline__ complex complex::operator*(complex param){
  complex temp;
  temp.real=real*param.real-imag*param.imag;
  temp.imag=real*param.imag+imag*param.real;
  return(temp);
}
__host__ __device__ __inline__ complex complex::operator*(double param){
  complex temp;
  temp.real=real*param;
  temp.imag=imag*param;
  return(temp);
}
__host__ __device__ __inline__ complex conjg(complex a){
  complex temp;
  temp.real=a.real;
  temp.imag=-a.imag;
  return(temp);
}
__host__ __device__ __inline__ double norm(complex a){
  return( sqrt( a.real*a.real+a.imag*a.imag) );}
/*---------------------------------------*/
#endif /* !FILE_mycomplex_SEEN */

