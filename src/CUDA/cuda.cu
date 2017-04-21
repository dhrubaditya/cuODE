#include "CUDA.h"
/* ------------------------ */
double *allocdev_double(int len){
  double *dev_array;
  cudaMalloc( (void**)&dev_array, len*sizeof(double) );
  return dev_array;
}
int *allocdev_int(int len){
  int *dev_array;
  cudaMalloc( (void**)&dev_array, len*sizeof(int) );
  return dev_array;
}
curandState *host2dev(int len, curandState host_array[]){
  curandState *dev_array;
  cudaMalloc( (void**)&dev_array, len*sizeof(curandState) );
  cudaMemcpy(dev_array, host_array, len*sizeof(curandState), cudaMemcpyHostToDevice);
  return dev_array;
}
double *host2dev(int len, double host_array[]){
  double *dev_array;
  cudaMalloc( (void**)&dev_array, len*sizeof(double) );
  cudaMemcpy(dev_array, host_array, len*sizeof(double), cudaMemcpyHostToDevice);
  return dev_array;
}
int *host2dev(int len, int host_array[]){
  int *dev_array;
  cudaMalloc( (void**)&dev_array, len*sizeof(int) );
  cudaMemcpy(dev_array, host_array, len*sizeof(int), cudaMemcpyHostToDevice);
  return dev_array;
}
unsigned long long  *host2dev(int len, unsigned long long host_array[]){
  unsigned long long *dev_array;
  cudaMalloc( (void**)&dev_array, len*sizeof(unsigned long long) );
  cudaMemcpy(dev_array, host_array, len*sizeof(unsigned long long), cudaMemcpyHostToDevice);
  return dev_array;
}
void  h2d(double *dev_array, int len, double host_array[]){
  cudaMemcpy(dev_array, host_array, len*sizeof(double), cudaMemcpyHostToDevice);
}
void  h2d(int *dev_array, int len, int host_array[]){
  cudaMemcpy(dev_array, host_array, len*sizeof(int), cudaMemcpyHostToDevice);
}
void dev2host(curandState host_array[], int len, curandState *dev_array){
  cudaMemcpy(host_array, dev_array, len*sizeof(curandState),cudaMemcpyDeviceToHost);
}
void dev2host(double host_array[], int len, double *dev_array){
  cudaMemcpy(host_array, dev_array, len*sizeof(double),cudaMemcpyDeviceToHost);
}
void dev2host(int host_array[], int len, int *dev_array){
  cudaMemcpy(host_array, dev_array, len*sizeof(int),cudaMemcpyDeviceToHost);
}
void dev2host(unsigned long long host_array[], int len, unsigned long long *dev_array){
  cudaMemcpy(host_array, dev_array, len*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
}
/*---------------------------------------*/
