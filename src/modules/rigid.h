/* ----------------------------------------------- */
unsigned int const Nensemble=1;
double const kk=4.;
int Utype = 1; //ABC flow
// Initial rotation of the body about the space-fixed axis 
int q_ini = 2; // rotated about x axis by psi_zero
double psi_zero=pi/4.;
// Initial angular momentum 
int ell_ini = 2; // angular *velocity* along z axis
double omega_zero=0.1; 
// Time step
double const TMAX=1000.;
double dt=0.01;
/* ----------------------------------------*/
