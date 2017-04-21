#ifndef FILE_matrix_SEEN
#define FILE_matrix_SEEN
/*---------------------------------------*/
#include <iostream>
#include <math.h>
#include "3vec.h"
using namespace std;
double const UnitM[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
/*---------------------------------------*/
void InvertDiagMat(double temp[3][3], double M[3][3]){
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      temp[i][j]=M[i][j];
    }
  }
  temp[0][0]=1./M[0][0];
  temp[1][1]=1./M[1][1];
  temp[2][2]=1./M[2][2];
}
/*---------------------------------------*/
void MatTrans(double temp[3][3], double M[3][3]){
  //
  temp[0][0]=M[0][0];
  temp[1][1]=M[1][1];
  temp[2][2]=M[2][2];
  //
  temp[0][1]=M[1][0];
  temp[0][2]=M[2][0];
  //
  temp[1][0]=M[0][1];
  temp[1][2]=M[2][1];
  //
  temp[2][0]=M[0][2];
  temp[2][1]=M[1][2];
}
/*---------------------------------------*/
vec3 MatVec(double M[3][3], vec3 B){
  vec3 A;
  double AdotB[3];
  for(int i=0; i< 3; i++){
    A = vec3(M[i][0],M[i][1],M[i][2]);
    AdotB[i] = dot(A,B);
  }
  return vec3(AdotB[0],AdotB[1],AdotB[2]);
}
/*---------------------------------------*/
void MatMul(double temp[3][3], double A[3][3], double B[3][3]){
  vec3 Cij,Bj;
  for(int i=0; i<3; i++){
    Bj = vec3(B[0][i],B[1][i],B[2][i]);
    Cij = MatVec(A, Bj);
    temp[0][i]=Cij.x;
    temp[1][i]=Cij.y;
    temp[2][i]=Cij.z;
  }
}
/*---------------------------------------*/
void RMRt(double Mrot[3][3], double M[3][3], double R[3][3]){
  double Rt[3][3], temp[3][3];
  MatTrans(Rt,R);
  MatMul(temp,M,Rt);
  MatMul(Mrot,R,temp);
}
/*---------------------------------------*/
void RtMR(double Mrot[3][3], double M[3][3], double R[3][3]){
  double Rt[3][3], temp[3][3];
  MatTrans(Rt,R);
  MatMul(temp,M,R);
  MatMul(Mrot,Rt,temp);
}
/*---------------------------------------*/
void PMat(const double M[3][3]){
  cout<<M[0][0]<<"\t"<<M[0][1]<<"\t"<<M[0][2]<<"\n";
  cout<<M[1][0]<<"\t"<<M[1][1]<<"\t"<<M[1][2]<<"\n";
  cout<<M[2][0]<<"\t"<<M[2][1]<<"\t"<<M[2][2]<<"\n";
}
/*---------------------------------------*/
#endif /* !FILE_matrix_SEEN */

