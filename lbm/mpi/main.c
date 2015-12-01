#include "stdio.h"
#include "stdlib.h"
#include "math.h"
// #include "mpi.h"
//-----------------------------------------------------------//
// #include "constants.h"
// #include "variables.h"
// #include "functions.h"
// #include "MPI_def.h"
//------------------------------------------------------------//
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
  #define rho_in      1.01
  #define rho_out     1.0
  #define inv_rho_in  1./rho_in
  #define inv_rho_out 1./rho_out
#endif
#define tol     10e-6
#define tsteps  100000 // no. of iterations
#define savet   1000   //file saving steps
#define ftag    1

//--------------------------------------------------------------//
// #define model_1  //compressible flow
#define model_2  //incompressible flow
// #define mdoel_3  //modified incompressible flow
//--------------------------------------------------------------//
double *f, *f_str;
double feq[Q];
double *u, *v, *rho;
double w[Q]    = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
double e[Q][2] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0};
int t;
//--------------------------------------------------------------//
void    allocate_memory();
void    init();

void    calculate_rho();
void    calculate_velocities();
void    Fk_eq( double u, double v, double r );
void    collision_step();
void    streaming_step();

void    boundary();
void    boundary_pipeflow();

void    write2file1( int t, int m, double *c, double *d, char *fname);  // for vector array - velocity
void    write2file( int t, int m, double *c, char *fname);              // for scalar array - phi
void    free_memory();
//-------------------------------------------------------------//
void main(int argc, char *argv[]){
  int i,j,k;
  char fname1[100], fname2[100];
  int  m = sprintf(fname1,"velocities");
  int  n = sprintf(fname2,"rho");
  allocate_memory();
  init();
  // boundary();
  write2file1(0,Mx,u,v,fname1);
  for (t = 1; t <= tsteps; t++) {
    collision_step();
    streaming_step();
    boundary_pipeflow();
    calculate_rho();
    calculate_velocities();
    if(t%savet == 0){
      write2file1(t,Mx,u,v,fname1);
      // write2file(t,Mx,rho,fname2);
    }
    printf("iteration no. %d \n",t);
  }
  free_memory();
}

void allocate_memory() {
  //velocity matrices
  u      = (double *)malloc(Mx*My*sizeof(double));
  v      = (double *)malloc(Mx*My*sizeof(double));
  rho    = (double *)malloc(Mx*My*sizeof(double));
  // probability distribution functions
  f     = (double *)malloc(Q*Mx*My*sizeof(double));
  f_str = (double *)malloc(Q*Mx*My*sizeof(double));
}
void free_memory() {
  free(u);
  free(v);
  free(rho);
  free(f);
  free(f_str);
}
void calculate_rho(){
  int i,j,k,indx, indx_f;
  double tmp = 0.0;
  for(i = 1; i < Mx-1 ; i++) {
    for(j = 1; j < My-1 ; j++) {
      tmp  = 0.0;
      indx = i*Mx + j;
      for(k = 0; k < Q; k++ ){
        indx_f = k*M2 + indx;
        tmp += f[indx_f];
      }
      rho[indx] = tmp;
    }
  }
}
void calculate_velocities() {
  int i, j, k, indx_f, indx;
  double tmp1, tmp2, inv_rho;
  for(i = 1; i < Mx-1; i++) {
    for(j = 1; j < My-1; j++) {
      indx = i*Mx + j;
      tmp1 = 0.0;
      tmp2 = 0.0;
      for(k = 0; k < Q; k++ ){
        indx_f  = k*M2 + indx;
        tmp1   += e[k][0]*f[indx_f];
        tmp2   += e[k][1]*f[indx_f];
      }
    #ifdef model_1
      u[indx] = tmp1/rho[indx];
      v[indx] = tmp2/rho[indx];
    #endif
    #ifdef model_2
      u[indx] = tmp1;
      v[indx] = tmp2;
    #endif
    }
  }
}
void Fk_eq( double u, double v, double r ) {
  int k;
  double e_dot_u, u_dot_u;
  u_dot_u = u*u+v*v;
  for (k = 0; k < Q; k++) {
    e_dot_u = e[k][0]*u+e[k][1]*v;
  #ifdef model_1
    feq[k] = w[k]*r*(1. + 3.0*(e_dot_u) - 1.5*(u_dot_u) + 4.5*(e_dot_u*e_dot_u));
  #endif
  #ifdef model_2
    feq[k] = w[k]*(r + 3.0*(e_dot_u) - 1.5*(u_dot_u) + 4.5*(e_dot_u*e_dot_u));
  #endif
  }
}
void collision_step() {
  int i, j, k, indx_f, indx;
  for(i = 1; i < Mx-1 ; i++) {
    for(j = 1; j < My-1 ; j++) {
      indx   = i*My + j;
      Fk_eq( u[indx], v[indx], rho[indx] );
      for( k = 0; k < Q; k++ ){
        indx_f        = k*M2 + indx;
        f_str[indx_f] = omega*feq[k] + (1.0 - omega)*f[indx_f];
      }
    }
  }
}
void streaming_step() {
  int i, j, k, indx, indx_f, indx_next;
  for(i = 1; i < Mx - 1; i++) {
    for(j = 1; j < My - 1; j++) {
      indx = i*My + j;
      for(k = 0; k < Q; k++ ){
        indx_f       = k*M2 + indx;
        indx_next    = indx_f + e[k][1]*My + e[k][0];
        f[indx_next] = f_str[indx_f];
      }
    }
  }
}
void init() {
  int i, j, k, indx, indx_f;
  for(i = 0; i < Mx; i++) {
    for(j = 0; j < My; j++) {
      indx = i*My+j;
      rho[indx] = Rho_init;
      u[indx]   = 0.0;
      v[indx]   = 0.0;
      Fk_eq(u[indx],v[indx],rho[indx]);
      for(k = 0; k < Q; k++ ){
        indx_f = k*M2 + indx;
        f[indx_f]     = feq[k];
        f_str[indx_f] = 0.0;
      }
    }
  }
}
void    boundary_pipeflow(){
  int i, j, indx_left, indx_right, indx_top, indx_bot;
  int indx_next;
  double u_wall;
  double two_by_three = 0.6666666666667;
  double one_by_six   = 0.1666666666667;
  for(i = 2; i < Mx-2; i++) {
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  // left wall - inlet - eastside
    indx_left    = My*i + 1;
    u_wall       = rho_in - (f[0*M2+indx_left]+f[2*M2+indx_left]+f[4*M2+indx_left]+2.*(f[3*M2+indx_left]+f[7*M2+indx_left]+f[6*M2+indx_left]));
    u[indx_left] = u_wall;
    v[indx_left] = 0.0;
    rho[indx_left]    = rho_in;
    f[1*M2+indx_left] = f[3*M2+indx_left] + two_by_three*u_wall;//
    f[5*M2+indx_left] = f[7*M2+indx_left] - 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;
    f[8*M2+indx_left] = f[6*M2+indx_left] + 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;//
  }
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  for(i = 2; i < Mx-2; i++) {
  // right wall
    indx_right   = My*i + My - 2;
    u_wall       = -1.0*rho_out + (f[0*M2+indx_right]+f[2*M2+indx_right]+f[4*M2+indx_right]+2.*(f[1*M2+indx_right]+f[5*M2+indx_right]+f[8*M2+indx_right]));

    u[indx_right] = u_wall;
    v[indx_right] = 0.0;
    rho[indx_right]    = rho_out;

    f[3*M2+indx_right] = f[1*M2+indx_right] - two_by_three*u_wall;
    f[7*M2+indx_right] = f[5*M2+indx_right] + 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;
    f[6*M2+indx_right] = f[8*M2+indx_right] - 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;//
  }
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  for (i = 2; i < Mx - 2; i++) {
  // bottom wall
    indx_bot          = i + My;

    u[indx_bot] = 0.0;
    v[indx_bot] = 0.0;
    rho[indx_bot]     = f[0*M2+indx_bot]+f[1*M2+indx_bot]+f[3*M2+indx_bot]+2.0*(f[4*M2+indx_bot]+f[7*M2+indx_bot]+f[8*M2+indx_bot]);

    f[2*M2+indx_bot]  = f[4*M2+indx_bot];
    f[5*M2+indx_bot]  = f[7*M2+indx_bot] - 0.5*(f[1*M2+indx_bot] - f[3*M2+indx_bot]);
    f[6*M2+indx_bot]  = f[8*M2+indx_bot] + 0.5*(f[1*M2+indx_bot] - f[3*M2+indx_bot]);
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  // top wall
    indx_top         = i + (Mx-2)*My;

    u[indx_top] = 0.0;
    v[indx_top] = 0.0;
    rho[indx_top]     = f[0*M2+indx_top]+f[1*M2+indx_top]+f[3*M2+indx_top]+2.0*(f[2*M2+indx_top]+f[5*M2+indx_top]+f[6*M2+indx_top]);


    f[4*M2+indx_top] = f[2*M2+indx_top];
    f[7*M2+indx_top] = f[5*M2+indx_top] + 0.5*(f[1*M2+indx_top] - f[3*M2+indx_top]);
    f[8*M2+indx_top] = f[6*M2+indx_top] - 0.5*(f[1*M2+indx_top] - f[3*M2+indx_top]);
  }
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  // top left corner
  indx_top         = 1 + (Mx-2)*My;
  u[indx_top] = 0.0;
  v[indx_top] = 0.0;
  f[1*M2+indx_top] = f[3*M2+indx_top];
  f[4*M2+indx_top] = f[2*M2+indx_top];
  f[8*M2+indx_top] = f[6*M2+indx_top];
  rho[indx_top]     = f[0*M2+indx_top]+f[1*M2+indx_top]+f[3*M2+indx_top]+2.0*(f[2*M2+indx_top]+f[5*M2+indx_top]+f[6*M2+indx_top]);
  f[5*M2+indx_top] = 0.5*(rho_in - (f[0*M2+indx_top]+f[1*M2+indx_top]+f[2*M2+indx_top]+f[3*M2+indx_top]+f[4*M2+indx_top]+f[6*M2+indx_top]+f[8*M2+indx_top]));
  // f[5*M2+indx_top] = 0.5*(1. -f[0*M2+indx_top] - f[3*M2+indx_top] - f[2*M2+indx_top] - f[6*M2+indx_top] - f[4*M2+indx_top] - f[8*M2+indx_top]);
  f[7*M2+indx_top] = f[5*M2+indx_top];
  // bottom left corner
  indx_bot          = 1 + My;
  f[1*M2+indx_bot]  = f[3*M2+indx_bot];
  f[2*M2+indx_bot]  = f[4*M2+indx_bot];
  f[5*M2+indx_bot]  = f[7*M2+indx_bot];
  f[6*M2+indx_bot]  = 0.5*(rho_in - (f[0*M2+indx_bot]+f[1*M2+indx_bot]+f[2*M2+indx_bot]+f[3*M2+indx_bot]+f[4*M2+indx_bot]+f[5*M2+indx_bot]+f[7*M2+indx_bot]));
  // f[6*M2+indx_bot]  = 0.5*(1. - (f[0*M2+indx_bot]+f[1*M2+indx_bot]+f[2*M2+indx_bot]+f[3*M2+indx_bot]+f[4*M2+indx_bot]+f[5*M2+indx_bot]+f[7*M2+indx_bot]));
  f[8*M2+indx_bot]  = f[6*M2+indx_bot];
  //top right corner
  indx_top         = My - 2 + (Mx-2)*My;
  u[indx_top] = 0.0;
  v[indx_top] = 0.0;
  f[3*M2+indx_top] = f[1*M2+indx_top];
  f[4*M2+indx_top] = f[2*M2+indx_top];
  f[7*M2+indx_top] = f[5*M2+indx_top];
  rho[indx_top]     = f[0*M2+indx_top]+f[1*M2+indx_top]+f[3*M2+indx_top]+2.0*(f[2*M2+indx_top]+f[5*M2+indx_top]+f[6*M2+indx_top]);
  f[6*M2+indx_top] = 0.5*(rho[indx_top] - (f[0*M2+indx_top]+f[1*M2+indx_top]+f[2*M2+indx_top]+f[3*M2+indx_top]+f[4*M2+indx_top]+f[5*M2+indx_top]+f[7*M2+indx_top]));
  f[8*M2+indx_top] = f[6*M2+indx_top];
  // Bottom right corner
  indx_bot          = My - 2 + My;
  f[3*M2+indx_bot]  = f[1*M2+indx_bot];
  f[2*M2+indx_bot]  = f[4*M2+indx_bot];
  f[6*M2+indx_bot]  = f[8*M2+indx_bot];
  f[5*M2+indx_bot]  = 0.5*(rho_out - (f[0*M2+indx_bot]+f[1*M2+indx_bot]+f[2*M2+indx_bot]+f[3*M2+indx_bot]+f[4*M2+indx_bot]+f[6*M2+indx_bot]+f[8*M2+indx_bot]));
  f[7*M2+indx_bot]  = f[5*M2+indx_bot];

//*****************************************************************************************//
}
void write2file1( int t, int m, double *c, double *d, char *fname) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles%d/%s_%d.dat",ftag, fname, t);
  fp = fopen(filename,"w");
  for ( i = 0; i < m; i++)
  {
    for ( j = 0; j < m; j++)
    {
      z= i*m + j;
      fprintf(fp,"%d %d %le %le\n",j,i,c[z],d[z]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
void write2file( int t, int m, double *c, char *fname) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles%d/%s_%d.dat",ftag,fname, t);
  fp = fopen(filename,"w");
  for ( i = 0; i < m; i++)
  {
    for ( j=0; j < m; j++)
    {
      z= i*m + j;
      fprintf(fp,"%d %d %le\n",j,i,c[z]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
#ifdef ldc
void    boundary(){
  int i, j, indx_left, indx_right, indx_top, indx_bot;
  int indx_next;
  double u_wall = u_lid, v_wall = 0.0, rho_n;
  for(i = 1; i < Mx-2; i++) {
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  // left wall
    indx_left         = My*i + 1;
    f[1*M2+indx_left] = f[3*M2+indx_left];//

    f[5*M2+indx_left] = f[7*M2+indx_left] - 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]);

    f[8*M2+indx_left] = f[6*M2+indx_left] + 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]);//
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  // right wall
    indx_right         = My*i + My - 2;
    f[3*M2+indx_right] = f[1*M2+indx_right];

    f[7*M2+indx_right] = f[5*M2+indx_right] + 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]);

    f[6*M2+indx_right] = f[8*M2+indx_right] - 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]);//
  /*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  }
    // bottom wall
  for (i = 1; i < Mx - 1; i++) {
    indx_bot          = i + My;

    f[2*M2+indx_bot]  = f[4*M2+indx_bot];

    f[5*M2+indx_bot]  = f[7*M2+indx_bot] - 0.5*(f[1*M2+indx_bot] - f[3*M2+indx_bot]);

    f[6*M2+indx_bot]  = f[8*M2+indx_bot] + 0.5*(f[1*M2+indx_bot] - f[3*M2+indx_bot]);
  }
  for ( i = 2; i < Mx - 2; i++) {
  // top wall
    indx_top = i + (Mx-2)*My;
    rho_n    = (f[0*M2+indx_top]+f[1*M2+indx_top]+f[3*M2+indx_top]+2*(f[2*M2+indx_top]+f[5*M2+indx_top]+f[6*M2+indx_top]));
    rho_n = 1;
    f[4*M2+indx_top] = f[2*M2+indx_top];
    f[7*M2+indx_top] = f[5*M2+indx_top] + 0.5*(f[1*M2+indx_top] - f[3*M2+indx_top]) - rho_n*u_lid*0.5;
    f[8*M2+indx_top] = f[6*M2+indx_top] - 0.5*(f[1*M2+indx_top] - f[3*M2+indx_top]) + rho_n*u_lid*0.5;
  }
  // top left corner
  indx_top = 1 + (Mx-2)*My;
  f[1*M2+indx_top] = f[3*M2+indx_top] + (2.*u_lid/3);
  f[4*M2+indx_top] = f[2*M2+indx_top];
  f[8*M2+indx_top] = f[6*M2+indx_top] + (u_lid)/6.;
  f[5*M2+indx_top] = -0.5*f[0*M2+indx_top] - f[3*M2+indx_top] - f[2*M2+indx_top] - f[6*M2+indx_top] + 0.5*(1-2.*u_lid/3);
  f[7*M2+indx_top] = -0.5*f[0*M2+indx_top] - f[3*M2+indx_top] - f[2*M2+indx_top] - f[6*M2+indx_top] + 0.5*(1.-u_lid) ;

  //top right corner
  indx_top = My-2 + (Mx-2)*My;
  f[3*M2+indx_top] = f[1*M2+indx_top];
  f[2*M2+indx_top] = f[4*M2+indx_top];
  f[7*M2+indx_top] = f[5*M2+indx_top];
  f[8*M2+indx_top] = -0.5*f[0*M2+indx_top] - f[1*M2+indx_top] - f[5*M2+indx_top] - f[7*M2+indx_top] + 0.5*(1.);
  f[6*M2+indx_top] = -0.5*f[0*M2+indx_top] - f[1*M2+indx_top] - f[5*M2+indx_top] - f[7*M2+indx_top] + 0.5*(1.) ;
}
#endif
