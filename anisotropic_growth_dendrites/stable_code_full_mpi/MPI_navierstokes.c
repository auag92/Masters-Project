#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "MPI_def.h"

void allocate_memory_fluid();
void velocities_initialize();
void boundary_velocities();
void boundary_pressure();
void ns_solver(int start, int end, double *phi);

void NS_solver() {
  allocate_memory_fluid();
  if (taskid == MASTER) {
    velocities_initialize();
    boundary_velocities();
    fluid_distribute();
  }else {
    fluid_distribute();
    computeH(u_old,v_old,Hx,Hy);
    RHS_fn(Hx,Hy,rhs_fn,MESHX);
    LHS_fn(MESHX, pmesh, rows);
    boundary_pressure(taskid);
    gs_mpi();
    ns_solver(start, end, phi_old);
    V_update(MESHX);
  }
}

void initialize_velocities() {
  long i,j,z, r1,r2;

  for (i=0; i< MESHX; i++)
  {
    for (j=0; j< MESHX; j++)
    {

      z = i*MESHX+j;
      v_old[z] = 0.0;
      u_old[z] = 0.0;
    }
  }
}
void initialize_pressure(double *P,int taskid, int Px) {
  int i, j;
  if(taskid == 1) {
    for (i=0; i < rows; i++){
      for (j=0; j < Px; j++){
        z= i*Px + j;
        P[z] = 0.0;
        rhs_fn[z] = 0.0;
        if(j==0) {
          P[z] = p_left;
        }else if(j == (Px-1)){
          P[z] = p_right;
        }
        if(i == 0) {
          P[z] = p_down;
        }
      }
    }
  }else if( taskid == numworkers ){
    for (i=0; i < rows; i++){
      for (j=0; j < Px; j++){
        z= i*Px + j;
        P[z] = 0.0;
        rhs_fn[z] = 0.0;
        if(j==0) {
          P[z] = p_left;
        }else if(j == (Px-1)){
          P[z] = p_right;
        }
        if(i == rows-1) {
          P[z] = p_up;
        }
      }
    }
  }else {
    for (i = 0; i < rows+1; i++){
      for (j = 0; j < Px; j++){
        z= i*Px + j;
        P[z] = 0.0;
        rhs_fn[z] = 0.0;
        if(j==0) {
          P[z] = p_left;
        }else if(j == (Px-1)){
          P[z] = p_right;
        }
      }
    }
  }
}
void boundary_velocities() {
  int i;
  for (i=0; i<MESHX; i++)
  {
#ifdef PipeFlow
    v_old[i] = Vs; // south
    v_old[i*MESHX] = Ve; // east
    v_old[i*MESHX+MESHX-1] = Vw; // west
    v_old[MESHX2-MESHX+i] = Vn; // north

    u_old[i] = Us; // south
    //u_old[i*MESHX] = Ue; // east
    //u_old[i*MESHX+MESHX-1] = Uw; // west
    u_old[MESHX2-MESHX+i] = Un;  // north
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
void LHS_fn(int Mx, int Px, int r){
  int i,j,x,m,x1,x2,y,y1,y2,z,z1,z2,c,d1,d2,t;
  int phi_indx,a_indx;

  t = pmesh-1;
  for (i=0;i<t-1;i++){
    for (j=0;j<t;j++){
      phi_indx = (i+1)*MESHX + j+1;
      a_indx = i*t +j;
      a_x[a_indx] = 1.0 - 0.5*(phi_old[phi_indx] + phi_old[phi_indx+MESHX]);
    }
  }
  t = pmesh-2;
  for (i=0;i<t+1;i++){
    for (j=0;j<t;j++){
      phi_indx = (i+1)*MESHX + j+1;
      a_indx = i*t + j;
      a_y[a_indx] = 1.0 - 0.5*(phi_old[phi_indx] + phi_old[phi_indx+1]);
    }
  }
}
void laplacian(double *f, double *lap, int Mx, int start, int end) {
  long i,j,z;

  for (i = start; i <= end; i++)
  {
    for (j = 1; j < Mx - 1; j++)
    {
      z = i*Mx + j;
      lap[z] = (f[z-1] + f[z+1] + f[z+Mx] + f[z-Mx] -4.0*f[z])*inv_deltax2;
    }
  }
}
void computeH(double *u, double *v,double *hx, double *hy){

  int i,j,z,x;
  double du_dx,du_dy,dv_dx,dv_dy,du_dt,dv_dt;

  V_str(MESHX);
  laplacian(u,lap_u, MESHX, start, end);
  laplacian(v,lap_v, MESHX, start, end);
  for(i=start; i<= end; i++){
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
void ns_solver(int start, int end, double *phi) {

  int     i,j,z,x;
  double  dp_dx,dp_dy,du_dt,dv_dt;
  double  phi;

  for(i=start; i<=end; i++){
    for(j=1; j<MESHX-1; j++){
      z = i*MESHX +j;
      x = (i-1)*(pmesh)+ (j-1);

      phi = phi[z];

      dp_dx = 0.5*inv_deltax*(P[x+1]+P[x+pmesh+1]-P[x+pmesh]-P[x]);
      dp_dy = 0.5*inv_deltax*(P[x+pmesh]+P[x+pmesh+1]-P[x+1]-P[x]);

      du_dt = Hx[z] - dp_dx*(1.0-phi);
      dv_dt = Hy[z] - dp_dy*(1.0-phi);

      u_now[z] = u_str[z] + deltat * du_dt;
      v_now[z] = v_str[z] + deltat * dv_dt;
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
void V_str(int m, int start, int end){
  int i,j,z;
  for (i=start; i<=end; i++){
    for(j=0; j<m; j++){
      z = i*m+j;
      v_str[z]=v_old[z]*(1-phi_old[z]);
      u_str[z]=u_old[z]*(1-phi_old[z]);
    }
  }
}
