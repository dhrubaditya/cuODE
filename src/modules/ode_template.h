/* ------- Making the code work for both single and double precision -- */
#if defined(DOUBLE) || defined(OUBLE) /* So -DOUBLE works */
typedef double R;
#else
typedef float R;
#endif
/* ----------- templates ------------ */
template <typename Y>
Y euler(Type y, R time, R dt);
void eval_rhs(double time,double y[],double rhs[]);
void rnkt4(unsigned int ndim, double *y, double time,double dt);
void euler(unsigned int ndim, double *y, double time,double dt);
