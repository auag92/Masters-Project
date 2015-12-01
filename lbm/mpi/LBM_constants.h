#define dx       1
#define dy       dx
#define dt       1.0
#define Mx       100
#define My       Mx
#define M2       Mx*My
#define Q        9   // no. of nodes in the LBM model
#define nu       1.0 // vsicosity
#define Re       u_lid*Mx/nu
#define inv_tau  (1.*dt)/(3.*nu*dt + 0.5)//(deltax*Cs)/(3.*nu+0.5*Cs*deltax)
#define omega    0.9 //dt*inv_tau
#define Rho_init 1.0
// #define ldc
#ifdef ldc
  #define u_lid   0.5
#endif
#define pipeflow
#ifdef pipeflow
  #define rho_in      1.05
  #define rho_out     1.0
  #define inv_rho_in  1./rho_in
  #define inv_rho_out 1./rho_out
#endif
#define tol     10e-6
#define tsteps  3000 // no. of iterations
#define savet   10   //file saving steps
#define ftag    6

//--------------------------------------------------------------//
// #define model_1  //compressible flow
#define model_2  //incompressible flow
// #define mdoel_3  //modified incompressible flow
//--------------------------------------------------------------//
