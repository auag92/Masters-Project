#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "variables.h"
//#include "mg_solver.c"
#define ntimesteps 1
#define phi_tol 0.5

//#define LDC
#define PipeFlow

void update(double *old, double *now,int M);
void laplacian(double *f, double *lap, int M);
void fluid_initialize();
void boundary_fluid();
void computeH(double *u, double *v,double *hx, double *hy);
void RHS_fn(double *hx, double *hy, double *fn, int M);
void V_update(int m);
void V_str(int m);
void LHS_fn();
void Gauss_siedel(double *P, double *fn, double *a_x, double *a_y);
void boundary_pressure();

fluid_solver(){

  int i,j,z,x;
  double dp_dx,dp_dy,du_dt,dv_dt;
  double phi,phi_1,phi_2,phi_3,phi_4;

  boundary_fluid();
  boundary_pressure();
  computeH(u_old,v_old,Hx,Hy);
  RHS_fn(Hx,Hy,rhs_fn,MESHX);
  LHS_fn();
  //multigrid(P, rhs_fn, a_x, a_y);
  Gauss_siedel(P, rhs_fn, a_x, a_y);
  for(i=1; i<MESHX-1; i++){
    for(j=1; j<MESHX-1; j++){
	     z = i*MESHX +j;
	     x = (i-1)*(pmesh)+ (j-1);

	     phi = phi_old[z];

	     dp_dx = 0.5*inv_deltax*(P[x+1]+P[x+pmesh+1]-P[x+pmesh]-P[x]);
	     dp_dy = 0.5*inv_deltax*(P[x+pmesh]+P[x+pmesh+1]-P[x+1]-P[x]);

	     du_dt = Hx[z] - dp_dx*(1.0-phi);// - Mu*2.757*phi*phi*u_str[z]*inv_deltax2;
	     dv_dt = Hy[z] - dp_dy*(1.0-phi);// - Mu*2.757*phi*phi*v_str[z]*inv_deltax2;

	     u_now[z] = u_str[z] + deltat * du_dt;
	     v_now[z] = v_str[z] + deltat * dv_dt;
    }
  }
  V_update(MESHX);
  boundary_fluid();
}

void LHS_fn(){
  int i,j,x,m,x1,x2,y,y1,y2,z,z1,z2,c,d1,d2,t;

  t = pmesh-1;
  for (i=0;i<t;i++){
    for (j=0;j<t-1;j++){
      phi_indx = (i+1)*MESHX + j+1;
      a_indx = i*t +j;
      a_x[x] = 1.0 - 0.5*(phi_old[phi_indx] + phi_old[phi_indx+MESHX]);
    }
  }
  t = pmesh-2;
  for (i=0;i<t;i++){
    for (j=0;j<t+1;j++){
      phi_indx = (i+1)*MESHX + j+1;
      a_indx = i*t + j;
      a_y[y] = 1.0 - 0.5*(phi_old[z2] + phi_old[z2+1]);
    }
  }
}


void update(double *old, double *now,int M) {
  long i, j, z;
  for (i=0; i < M; i++) {
    for (j=0; j < M; j++){
      z= i*M + j;
      old[z]=now[z];
      old[z]=now[z];
    }
  }
}
void laplacian(double *f, double *lap, int M) {
  long i,j,z;

  for (i=1; i< M -1; i++)
  {
    for (j=1; j< M -1; j++)
    {
      z= i*M + j;
      lap[z] = (f[z-1] + f[z+1] + f[z+M] + f[z-M] -4.0*f[z])*inv_deltax2;
    }
  }
}
void boundary_fluid() {
  int i;
  for (i=1; i<MESHX-1; i++)
  {
#ifdef PipeFlow
    v_old[i] = Vs; // south
    v_old[i*MESHX] = Ve; // east
    v_old[i*MESHX+MESHX-1] = Vw; // west
    v_old[MESHX2-MESHX+i] = Vn; // north

    //u_old[i] = Us; // south
    u_old[i*MESHX] = Ue; // east
    u_old[i*MESHX+MESHX-1] = Uw; // west
    //u_old[MESHX2-MESHX+i] = Un;  // north
#endif
#ifdef LDC
    v_old[i] =0.0; // south
    v_old[i*MESHX] = 0.0; // east
    v_old[i*MESHX+MESHX-1] = 0.0; // west
    v_old[MESHX2-MESHX+i] = 0.0; // north

    u_old[i] = 2; // south
    u_old[i*MESHX] = 0; // east
    u_old[i*MESHX+MESHX-1] = 0; // west
    u_old[MESHX2-MESHX+i] = 0;  // north
#endif
  }
}
void computeH(double *u, double *v,double *hx, double *hy){

  int i,j,z,x;
  double du_dx,du_dy,dv_dx,dv_dy,du_dt,dv_dt;
  double lap_u[MESHX2], lap_v[MESHX2];
  V_str(MESHX);
  laplacian(u_str,lap_u, MESHX);
  laplacian(v_str,lap_v, MESHX);
  for(i=1; i<MESHX-1; i++){
    for(j=1; j<MESHX-1; j++){
      z = i*MESHX +j;
      du_dx = 0.5*inv_deltax*(u_str[z+1] - u_str[z-1]);
      du_dy = 0.5*inv_deltax*(u_str[z+MESHX] - u_str[z-MESHX]);
      dv_dx = 0.5*inv_deltax*(v_str[z+1] - v_str[z-1]);
      dv_dy = 0.5*inv_deltax*(v_str[z+MESHX] - v_str[z-MESHX]);
      hx[z] = inv_Re*lap_u[z] - (u_old[z]*du_dx + v_old[z]*du_dy);
      hy[z] = inv_Re*lap_v[z] - (u_old[z]*dv_dx + v_old[z]*dv_dy);
    }
  }
}
void RHS_fn(double *hx, double *hy, double *fn, int M){
  int i,j,x,z;
  double dhx_dx,dhy_dy,phi, dp_dx, dp_dy, dphi_dx,dphi_dy;
  for (i=1; i<M-2; i++){
    for (j=1; j<M-2; j++){
      z = i*M + j;
      x = i*(M-1) + j;
      dhx_dx = 0.5*inv_deltax*(hx[z+1]-hx[z]+hx[z+1+M]-hx[z+M]);
      dhy_dy = 0.5*inv_deltax*(hy[z+1+M]-hy[z+1]+hy[z+M]-hy[z]);
      fn[x] = dhx_dx+dhy_dy;
    }
  }
}
void V_update(int m){
  int i,j,z;
  double rem;
  for (i=0; i<m; i++){
    for(j=0; j<m; j++){
      z = i*m+j;
      if(phi_old[z] < phi_tol){
	      rem = 1.0/(1.0-phi_old[z]);
	      v_old[z]=v_now[z]*rem;
	      u_old[z]=u_now[z]*rem;
      }
      else{
        v_old[z] = 0.0;
        u_old[z] = 0.0;
      }
    }
  }
}
void V_str(int m){
  int i,j,z;
  for (i=0; i<m; i++){
    for(j=0; j<m; j++){
      z = i*m+j;
      v_str[z]=v_old[z]*(1-phi_old[z]);
      u_str[z]=u_old[z]*(1-phi_old[z]);
    }
  }
}
void Gauss_siedel(double *P, double *fn, double *a_x, double *a_y) {
  long i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck;
  double P_old;
  double error;
  double tol=1.0e-6;

  for(;;) {
    for(i=1; i < pmesh-1; i++) {
      for(j=1;j < pmesh-1;j++) {
        indx         = i*pmesh      + j;
        indx_rght    = indx         + 1;
        indx_frnt    = indx         + MESHY;
        indx_lft     = indx         - 1;
        indx_bck     = indx         - MESHY;
        indx_ax      = i*(pmesh-1)  +j;
        indx_ay      = i*(pmesh-2)   +j;

        if (((i+j)%2) == 0) {
          P_old     = P[indx];

          P[indx]  = a_x[indx_ax-1]*P[indx_lft] + a_x[indx_ax]*P[indx_rght]
          + a_y[indx_ay-(pmesh-2)]*P[indx_bck] + a_y[indx_ay]*P[indx_frnt] + deltax2*fn[indx];
          P[indx] /= a_x[indx_ax-1]+a_x[indx_ax]+a_y[indx_ay-(pmesh-2)]+a_y[indx_ay];
        }
      }
    }
    for(i=1; i < pmesh-1; i++) {
      for(j=1;j < pmesh-1;j++) {
        indx         = i*pmesh      + j;
        indx_rght    = indx         + 1;
        indx_frnt    = indx         + MESHY;
        indx_lft     = indx         - 1;
        indx_bck     = indx         - MESHY;
        indx_ax      = i*(pmesh-1)  +j;
        indx_ay      = i*(pmesh-2)   +j;

        if (((i+j)%2) != 0) {
          P_old     = P[indx];

          P[indx]  = a_x[indx_ax-1]*P[indx_lft] + a_x[indx_ax]*P[indx_rght]
          + a_y[indx_ay-(pmesh-2)]*P[indx_bck] + a_y[indx_ay]*P[indx_frnt] + deltax2*fn[indx];
          P[indx] /= a_x[indx_ax-1]+a_x[indx_ax]+a_y[indx_ay-(pmesh-2)]+a_y[indx_ay];
        }
      }
    }
    error = compute_error(a_x,a_y,P);
    if (fabs(error) < tol) {
      break;
    }
  }
}
double compute_error(double *a_x, double *a_y, double *P) {
  double error=0.0;
  long indx, indx_rght, indx_frnt, indx_lft, indx_bck, i, j;
  for(i=1; i < pmesh-1; i++) {
    for(j=1;j < pmesh-1; j++) {
      indx  = i*pmesh    + j;
      indx_right         = indx      + 1;
      indx_front         = indx      + MESHY;
      indx_left          = indx      - 1;
      indx_back          = indx      - MESHY;
      indx_ax      = i*(pmesh-1)  +j;
      indx_ay      = i*(pmesh-2)   +j;

      error += fabs(a_x[indx_ax-1]*P[indx_lft] + a_x[indx_ax]*P[indx_rght]
      + a_y[indx_ay-(pmesh-2)]*P[indx_bck] + a_y[indx_ay]*P[indx_frnt] +
      a_x[indx_ax-1]+a_x[indx_ax]+a_y[indx_ay-(pmesh-2)]+a_y[indx_ay] - deltax2*fn[indx]);
    }
  }
  return(error);
}
void boundary_pressure(){

  int i ,y ,z;
  int m2 = m*m;
  for (i=0; i< pmesh; i++)
  {
    y= i*pmesh;
    z= i*pmesh + pmesh-1;

    #ifdef neumann_pressure
    //left - right
    P[y]        = P[y+1];
    P[z]        = P[z-1];
    // up - down
    P[i]          = P[MESHX + i];
    P[m2-m+i]     = P[m2-2*m+i];
    #endif
    #ifdef dirichlet
    //left - right
    P[y]        = p_left;
    P[z]        = p_right;
    // up - down
    P[i]        = p_up;
    P[pmesh*pmesh-pmesh+i]  = p_down;
    #endif
  }
}
