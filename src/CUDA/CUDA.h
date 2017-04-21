#ifndef FILE_CUDA_SEEN
#define FILE_CUDA_SEEN
#include <curand.h>
#include <curand_kernel.h>
/* ------------------------ */
extern curandState *host2dev(int len, curandState host_array[]);
extern double *host2dev(int len, double host_array[]);
extern int *host2dev(int len, int host_array[]);
extern unsigned long long  *host2dev(int len, unsigned long long host_array[]);
extern void dev2host(curandState host_array[], int len, curandState *dev_array);
extern void dev2host(double host_array[], int len, double *dev_array);
extern void dev2host(int host_array[], int len, int *dev_array);
extern void dev2host(unsigned long long host_array[], int len, unsigned long long *dev_array);
/* ------------------------ */
#endif /* !FILE_CUDA_SEEN */
