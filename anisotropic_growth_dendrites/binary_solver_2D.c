#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#define MESHX 100
#define MESHX2 (MESHX*MESHX)
#define deltat (0.02)
#define deltax (1)
#define deltax2 (deltax*deltax)

#define RADIUS2 100
#define K (0.5)             /*Partition Coefficient*/
#define G (1.0)             /*Surface Energy*/
#define M (1.0)             /*Mobility*/
#define E (6.0)             /*Relaxation factor - dimensions of length [m]*/
#define tau (1.0)
#define Dab (0.06)

#define ntimesteps (100000)
#define saveT (100)
#define deltaMu (0.4)
#define Mu (1.0)

// #define DIRICHLET
#define CONSTANT_X_0 (0.1)
#define CONSTANT_X_1 (0.5)
#define NEUMANN

#define STEP_PROFILE

//#define isotropy
#define anisotropy
double    phi_new[MESHX2], phi_old[MESHX2];
double    mu_new[MESHX2], mu_old[MESHX2];
double    lap_phi[MESHX2], lap_mu[MESHX2];
double    conc[MESHX2];
double    dphi_now[4*MESHX], dphi_next[4*MESHX];
double    inv_deltax2 = (1.0/deltax2);
double    inv_deltax = (1.0/deltax);

void      initialize();
void      boundary(double *c);
void      write2file ( long t);
void      update();
void      laplacian(double *f, double *lap);
void      concentration();
void      fnupdate();
void      grad_phi(long i, double *d_phi);
double    div_phi(long i);
void anisotropic_solverloop();

void main() {
  long t;
  initialize();
  boundary(phi_old);
  boundary(mu_old);

  for (t=0; t < ntimesteps; t++) {

    laplacian(phi_old, lap_phi);
    laplacian(mu_old, lap_mu);
    anisotropic_solverloop();
    boundary(phi_new);
    boundary(mu_new);
    //concentration();
    update();
    if((t%saveT) ==0) {
      write2file(t);
    }
  }
  printf("Yay! We are done.\n");
}
void update() {
  long i, j, z;
  for (i=0; i < MESHX; i++) {
    for (j=0; j < MESHX; j++) {
      z= i*MESHX + j;
      phi_old[z]=phi_new[z];
      mu_old[z]=mu_new[z];
    }
  }
}
void laplacian(double *f, double *lap) {
  long i,j,z;

  for (i=1; i< MESHX -1; i++)
  {
    for (j=1; j< MESHX -1; j++)
    {
      z= i*MESHX + j;
      lap[z] = (f[z-1] + f[z+1] + f[z+MESHX] + f[z-MESHX]  -4.0*f[z])*inv_deltax2;
    }
  }
 }
void concentration(){
  double p,u,h;
  int i, j, z;
  for ( i = 0; i < MESHX; i++)
  {
    for ( j = 0; j < MESHX; j++){
      z= i*MESHX + j;
      p = phi_new[z];
      u = mu_new[z];
      h = p*p*(3-2*p);
      conc[z] = u*(1-h) + (h)*K*u;
    }
  }
}
void initialize() {
  long i,j,z;
  double r;

  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      r= (i-MESHX/2)*(i-MESHX/2) + (j-MESHX/2)*(j-MESHX/2);
      z= i*MESHX + j;
      if(r < RADIUS2){
        phi_old[z] = 1.0;
      }
      else{
        phi_old[z] = 0.0;
      }
      mu_old[z] = Mu - deltaMu;

    }
  }
}
void boundary(double *c) {
  int i ,y ,z;

  for (i=0; i<MESHX -1; i++)
  {
    y= i*MESHX;
    z= i*MESHX + MESHX-1;

    //left - right
    c[y]        = c[y+1];
    c[z]  = c[z-1];

    // up - down
    c[i]        = c[MESHX + i];
    c[MESHX2 - MESHX + i]  = c[MESHX2 - 2*MESHX + i];
  }
}
void grad_phi(long i, double *d_phi){

 	long j,z;

  for(j=1; j<MESHX-1; j++){
    z = i * MESHX + j;
    if (i == MESHX -1 ){
      d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
      d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
    }
    else {
      d_phi[j] = (phi_old[z] - phi_old[z-1])*inv_deltax;
      d_phi[MESHX+j] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
      d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
      d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
     }
  }
  if (i != MESHX-1) {
    z = i * MESHX + j;
    d_phi[j] = (phi_old[z] - phi_old[z-1])*inv_deltax;
    d_phi[MESHX+j] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
  }
}

double dqdx( double phi_x, double phi_y) {

	double a, phi_x2, phi_x4, phi_y2, phi_y4, inv_phi;
	long z;
  double ans = 0;
  double part1, part2, part3, part4;
  phi_x2 = phi_x *phi_x;
  phi_y2 = phi_y *phi_y;
  phi_y4 = phi_y2 *phi_y2;
  phi_x4 = phi_x2 *phi_x2;

  if ((phi_x2> 1e-15) && (phi_y2> 1e-15)){
    inv_phi = 1/(phi_x2+phi_y2);
    part1 = (1-Dab*(3-4*(phi_x4+phi_y4)*inv_phi*inv_phi));
    part2 = 2*G*E*part1*part1*phi_x;
    part3 = 32*G*E*Dab*(phi_x2+phi_y2)*(part1);
    part4 = phi_x2*phi_x*inv_phi*inv_phi - phi_x*(phi_x4+phi_y4)*inv_phi*inv_phi*inv_phi;
    ans = part2 + part3*part4;
  }

  return ans;
}

double div_phi(long i){

  double ans;
	double x_next, x_now;
	double y_next, y_now;

  x_now = dqdx(dphi_now[i], dphi_now[i+MESHX]);
  x_next = dqdx(dphi_now[i+1], dphi_now[i+1+MESHX]);
  y_now = dqdx(dphi_now[i+2*MESHX], dphi_now[i+3*MESHX]);
  y_next = dqdx(dphi_next[i+2*MESHX], dphi_next[i+3*MESHX]);
	ans = ((x_next - x_now) + ( y_next - y_now))*inv_deltax;

  return ans;
}
void fnupdate()
{
  long i;

  for( i=0; i<MESHX; i++) {
    dphi_now[i] = dphi_next[i];
    dphi_now[MESHX+i] = dphi_next[MESHX+i];
    dphi_now[2*MESHX+i] = dphi_next[2*MESHX+i];
    dphi_now[3*MESHX+i] = dphi_next[3*MESHX+i];
  }
}
void write2file (long t) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles/phi_%ld.dat",t);
  fp = fopen(filename,"w");
  if (fp) {
    for ( i = 0; i < MESHX; i++)
    {
      for ( j=0; j < MESHX; j++)
      {
        z= i*MESHX + j;
        fprintf(fp,"%d %d %le\n",i,j, phi_new[z]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  else {
    printf("#Error:%s could not be opened to write",filename);
    exit(1);
    }
  }
void anisotropic_solverloop(){

  long      i, j, z;
  double    p,dp,du;
  double    drv_frce, alln_chn;
  double    Gamma;

  grad_phi(1, dphi_now);
  for (i=1; i < (MESHX-1); i++) {
    grad_phi(i+1, dphi_next);
    for (j=1; j < (MESHX-1); j++){

      z= i*MESHX + j;
      p = phi_old[z];
      //Gamma = 2*G*lap_phi[z];
      Gamma = div_phi(j);
      drv_frce = (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p);
      alln_chn = E*Gamma - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p);
      dp = deltat*(alln_chn + drv_frce)/(tau*E);
      phi_new[z] = p + dp;
      du = deltat*M*lap_mu[z] - (K-1)*mu_old[z]*6*p*(1-p)*dp;
      mu_new[z] = mu_old[z]  + du/(1+(K-1)*p*p*(3-2*p));
    }
    fnupdate();
  }
}
