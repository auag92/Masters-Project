void allocate_memory_FluidVariables(int Mx, int Px);
void initialize_velocities(int Mx);
void initialize_pressure(double *P,int taskid, int Px);
void computeH(double *u, double *v,double *hx, double *hy);
void RHS_fn(double *hx, double *hy, double *fn, int Mx, int Px);
void LHS_fn(int Mx, int Px);
void V_str(int Mx);
void update(double *old, double *now,int Mx);
void laplacian(double *f, double *lap, int M);
void update_velocities(int Mx);
void gauss_siedel();
void red_solver(double *P, double *fn, double *a_x, double *a_y);
void black_solver(double *P, double *fn, double *a_x, double *a_y);
double compute_error_mpi( double *P, double *fn, double *a_x, double *a_y);

void allocate_memory_FluidVariables(int Mx, int Px) {
  if( taskid == MASTER ) {
    v_old       =   (double *)malloc(Mx*Mx*sizeof(double));
    u_old       =   (double *)malloc(Mx*Mx*sizeof(double));
  }else {
    if((taskid ==1) || (taskid == numworkers)) {
      v_old =   (double *)malloc((rows+1)*Mx*sizeof(double));
      v_now =   (double *)malloc((rows+1)*Mx*sizeof(double));
      u_old =   (double *)malloc((rows+1)*Mx*sizeof(double));
      u_now =   (double *)malloc((rows+1)*Mx*sizeof(double));
      v_str =   (double *)malloc((rows+1)*Mx*sizeof(double));
      u_str =   (double *)malloc((rows+1)*Mx*sizeof(double));
    } else {
      v_old =   (double *)malloc((rows+2)*Mx*sizeof(double));
      v_now =   (double *)malloc((rows+2)*Mx*sizeof(double));
      u_old =   (double *)malloc((rows+2)*Mx*sizeof(double));
      u_now =   (double *)malloc((rows+2)*Mx*sizeof(double));
      v_str =   (double *)malloc((rows+2)*Mx*sizeof(double));
      u_str =   (double *)malloc((rows+2)*Mx*sizeof(double));

    }
    P       = (double *)malloc((rows+2)*Px*sizeof(double));
    rhs_fn  = (double *)malloc((rows+2)*Px*sizeof(double));
    a_x     = (double *)malloc(rows*Px*sizeof(double));
    a_y     = (double *)malloc(rows*Px*sizeof(double));
    lap_u   = (double *)malloc(rows*Mx*sizeof(double));
    lap_v   = (double *)malloc(rows*Mx*sizeof(double));
    Hx      = (double *)malloc(rows*Mx*sizeof(double));
    Hy      = (double *)malloc(rows*Mx*sizeof(double));
  }
}
void initialize_velocities(int Mx) {
  long i,j,z;
  for (i=0; i< Mx; i++)
  {
    for (j=0; j< Mx; j++)
    {
      z = i*Mx+j;
      v_old[z] = 0.0;
      u_old[z] = 0.0;
    }
  }
  for (i=0; i<Mx; i++)
  {
#ifdef PipeFlow
    v_old[i] = Vs; // south
    v_old[i*Mx] = Ve; // east
    v_old[i*Mx+Mx-1] = Vw; // west
    v_old[Mx*Mx-Mx+i] = Vn; // north

    u_old[i] = Us; // south
    //u_old[i*Mx] = Ue; // east
    //u_old[i*Mx+Mx-1] = Uw; // west
    u_old[Mx*Mx-Mx+i] = Un;  // north
#endif
#ifdef LDC
    v_old[i] =0.0; // south
    v_old[i*Mx] = 0.0; // east
    v_old[i*Mx+Mx-1] = 0.0; // west
    v_old[Mx*Mx-Mx+i] = 0.0; // north

    u_old[i] = 2; // south
    u_old[i*Mx] = 0; // east
    u_old[i*Mx+Mx-1] = 0; // west
    u_old[Mx*Mx-Mx+i] = 0;  // north
#endif
  }
}
void initialize_pressure(double *P,int taskid, int Px) {
  int i, j, z;
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
void computeH(double *u, double *v,double *hx, double *hy){
  int i,j,z,x;
  double du_dx, du_dy, dv_dx, dv_dy, du_dt, dv_dt;
  V_str(MESHX);
  laplacian(u,lap_u, MESHX);
  laplacian(v,lap_v, MESHX);
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
void V_str(int Mx){
  int i,j,z;
  for ( i = start; i <= end; i++ ){
    for( j = 0; j < Mx; j++ ){
      z = i*Mx+j;
      v_str[z]=v_old[z]*(1-phi_old[z]);
      u_str[z]=u_old[z]*(1-phi_old[z]);
    }
  }
}
void RHS_fn(double *hx, double *hy, double *fn, int Mx, int Px){
  int i,j,x,z;
  double dhx_dx,dhy_dy,phi, dp_dx, dp_dy, dphi_dx,dphi_dy;
  for ( i = start; i <= end; i++){
    for ( j = 1; j < Mx-2; j++){
      z = i*(Mx) + j;
      x = i*(Px) + j;
      dhx_dx = 0.5*inv_deltax*(hx[z+1]-hx[z]+hx[z+1+Mx]-hx[z+Mx]);
      dhy_dy = 0.5*inv_deltax*(hy[z+1+Mx]-hy[z+1]+hy[z+Mx]-hy[z]);
      fn[x] = dhx_dx + dhy_dy;
    }
  }
}
void LHS_fn(int Mx, int Px){
  int i,j,x,m,x1,x2,y,y1,y2,z,z1,z2,c,d1,d2,t;
  int phi_indx,a_indx;

  t = Px-1;
  for (i = 0; i <= end; i++){
    for (j = 0; j < t; j++){
      phi_indx = (i+1)*Mx + j+1;
      a_indx = i*t +j;
      a_x[a_indx] = 1.0 - 0.5*(phi_old[phi_indx] + phi_old[phi_indx+Mx]);
    }
  }
  t = Px-2;
  for (i = 0; i <= (end+1); i++){
    for (j = 0; j < t; j++){
      phi_indx = (i+1)*MESHX + j+1;
      a_indx = i*t + j;
      a_y[a_indx] = 1.0 - 0.5*(phi_old[phi_indx] + phi_old[phi_indx+1]);
    }
  }
}
void update(double *old, double *now,int Mx) {
  long i, j, z;
  for (i=start; i <= end; i++) {
    for (j=0; j < Mx; j++){
      z= i*Mx + j;
      old[z] = now[z];
      old[z] = now[z];
    }
  }
}
// void laplacian(double *f, double *lap, int Mx) {
//   long i,j,z;
//
//   for (i = start; i <= end; i++)
//   {
//     for (j = 1; j < Mx-1; j++)
//     {
//       z = i*Mx + j;
//       lap[z] = (f[z-1] + f[z+1] + f[z+Mx] + f[z-Mx] -4.0*f[z])*inv_deltax2;
//     }
//   }
// }
void update_velocities( int Mx){
  int i,j,z;
  double rem;
  for (i = start; i <= end; i++){
    for(j = 1; j < Mx-1; j++){
      z = i*Mx+j;
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
void ns_solver(int start, int end, double *Phi) {

  int     i,j,z,x;
  double  dp_dx,dp_dy,du_dt,dv_dt;
  double  phi;

  for(i=start; i<=end; i++){
    for(j=1; j<MESHX-1; j++){
      z = i*MESHX +j;
      x = (i-1)*(pmesh)+ (j-1);

      phi = Phi[z];

      dp_dx = 0.5*inv_deltax*(P[x+1]+P[x+pmesh+1]-P[x+pmesh]-P[x]);
      dp_dy = 0.5*inv_deltax*(P[x+pmesh]+P[x+pmesh+1]-P[x+1]-P[x]);

      du_dt = Hx[z] - dp_dx*(1.0-phi);
      dv_dt = Hy[z] - dp_dy*(1.0-phi);

      u_now[z] = u_str[z] + deltat * du_dt;
      v_now[z] = v_str[z] + deltat * dv_dt;
    }
  }
}
void gauss_siedel() {
  double error, err;
  int flag;
  if (taskid == MASTER) {
    int iter = 0;
    for (;;) {
      error=0.0;
      err = 0.0;
      for(rank = 1; rank <= numworkers; rank ++){
        MPI_Recv (&err,     1,     MPI_DOUBLE,    rank,   ERROR,  MPI_COMM_WORLD, &status);
        error  += err;
      }
      iter ++;
      printf( " %d, %le \n", iter, error );
      if ( error < gs_tol ) {
        flag = 1;
      } else {
        flag = 0;
      }
      for( rank = 1; rank <= numworkers; rank++ ) {
        msgtype = BREAK;
        dest = rank;
        MPI_Send (&flag,     1,     MPI_INT,    dest,   msgtype,  MPI_COMM_WORLD);
      }
      if (flag == 1) {
        break;
      }
    }
  }else{
    for (; ;) {
      mpiexchange(taskid, P, pmesh);
      red_solver(P, rhs_fn, a_x, a_y);
      mpiexchange(taskid, P, pmesh);
      black_solver(P, rhs_fn, a_x, a_y);
      error = compute_error_mpi(P, rhs_fn, a_x, a_y);

      MPI_Send(&error,     1,     MPI_DOUBLE,    MASTER,   ERROR,  MPI_COMM_WORLD);
      MPI_Recv(&flag,      1,     MPI_INT,       MASTER,   BREAK,  MPI_COMM_WORLD,  &status);

      if (flag == 1) {
        break;
      }
    }
  }
}
void red_solver(double *P, double *fn, double *a_x, double *a_y){
  int i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck, indx_ax, indx_ay;
  for(i=start; i <= end; i++) {
    for(j=1;j < pmesh-1;j++) {
      indx         = i*pmesh      + j;
      indx_rght    = indx         + 1;
      indx_frnt    = indx         + pmesh;
      indx_lft     = indx         - 1;
      indx_bck     = indx         - pmesh;
      indx_ax      = (i-1)*(pmesh-1)  +(j-1);
      indx_ay      = (i-1)*(pmesh-2)  +(j-1);
      if (((i+j)%2) == 0) {
        #ifdef const_coeff
          P[indx]  = -1.0*(P[indx_lft] + P[indx_rght]
          + P[indx_bck] + P[indx_frnt]) + deltax*deltax*fn[indx];
          P[indx] /= -1.0*4.0;
        #endif
        #ifdef var_coeff
          P[indx]  = -1.0*(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
          + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]) + deltax*deltax*fn[indx];
          P[indx] /= -1.0*(a_x[indx_ax+1]+a_x[indx_ax]+a_y[indx_ay+(pmesh-2)]+a_y[indx_ay]);
        #endif
      }
    }
  }
}
void black_solver(double *P, double *fn, double *a_x, double *a_y){
  int i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck, indx_ax, indx_ay;
  for(i=start; i <= end; i++) {
    for(j=1;j < pmesh-1;j++) {
      indx         = i*pmesh      + j;
      indx_rght    = indx         + 1;
      indx_frnt    = indx         + pmesh;
      indx_lft     = indx         - 1;
      indx_bck     = indx         - pmesh;
      indx_ax      = (i-1)*(pmesh-1)  + (j-1);
      indx_ay      = (i-1)*(pmesh-2)  + (j-1);
      if (((i+j)%2) != 0) {
        #ifdef const_coeff
          P[indx]  = -1.0*(P[indx_lft] + P[indx_rght]
          + P[indx_bck] + P[indx_frnt]) + deltax*deltax*fn[indx];
          P[indx] /= -1.0*4.0;
        #endif
        #ifdef var_coeff
          P[indx]  = -1.0*(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
          + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]) + deltax*deltax*fn[indx];
          P[indx] /= -1.0*(a_x[indx_ax+1]+a_x[indx_ax]+a_y[indx_ay+(pmesh-2)]+a_y[indx_ay]);
        #endif
      }
    }
  }
}
double compute_error_mpi( double *P, double *fn, double *a_x, double *a_y) {
  double error=0.0;
  int indx, indx_rght, indx_frnt, indx_lft, indx_bck, i, j,indx_ax,indx_ay;
  for(i=start; i <= end; i++) {
    for(j=1;j < pmesh-1; j++) {
      indx         = i*pmesh    + j;
      indx_frnt    = indx      + pmesh;
      indx_rght    = indx      + 1;
      indx_lft     = indx      - 1;
      indx_bck     = indx      - pmesh;
      indx_ax      = (i-1)*(pmesh-1)  + j-1;
      indx_ay      = (i-1)*(pmesh-2)  + j-1;
      #ifdef const_coeff
        error += fabs(P[indx_lft] + P[indx_rght]
        + P[indx_bck] + P[indx_frnt] - (4.0)*P[indx] - deltax*deltax*fn[indx]);
      #endif
      #ifdef var_coeff
        error += fabs(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
        + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]
        - (a_x[indx_ax]+a_x[indx_ax+1]+a_y[indx_ay]+a_y[indx_ay+(pmesh-2)])*P[indx]
        - deltax*deltax*fn[indx]);
      #endif
    }
  }
  return(error);
}
