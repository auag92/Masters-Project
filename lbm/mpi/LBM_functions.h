//--------------------------------------------------------------//
void    allocate_memory();
void    init();
void    assign_weights();
void    assign_eks();
void    calculate_rho();
void    calculate_velocities();
void    Fk_eq( double u, double v, double r );
void    collision_step();
void    streaming_step();
void    boundary();
void    boundary_pipeflow();
double  calc_rho(int indx);
void    write2file1( int t, int m, double *c, double *d, char *fname);  // for vector array - velocity
void    write2file( int t, int m, double *c, char *fname);              // for scalar array - phi
void    free_memory();
//-------------------------------------------------------------//
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
void assign_weights() {
  w[0] = 4.0/9.0;
  w[1] = 1.0/9.0;
  w[2] = 1.0/9.0;
  w[3] = 1.0/9.0;
  w[4] = 1.0/9.0;
  w[5] = 1.0/36.0;
  w[6] = 1.0/36.0;
  w[7] = 1.0/36.0;
  w[8] = 1.0/36.0;
}
void assign_eks() {
  e[0][0] =  0.0;
  e[0][1] =  0.0;
  e[1][0] =  1.0;
  e[1][1] =  0.0;
  e[2][0] =  0.0;
  e[2][1] =  1.0;
  e[3][0] = -1.0;
  e[3][1] =  0.0;
  e[4][0] =  0.0;
  e[4][1] = -1.0;
  e[5][0] =  1.0;
  e[5][1] =  1.0;
  e[6][0] = -1.0;
  e[6][1] =  1.0;
  e[7][0] = -1.0;
  e[7][1] = -1.0;
  e[8][0] =  1.0;
  e[8][1] = -1.0;
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

  if(taskid == 1){
    for (i = 2; i < My-2; i++) {
    // bottom wall - south side - no slip
      indx_bot          = i + My;
      u[indx_bot] = 0.0;
      v[indx_bot] = 0.0;
      rho[indx_bot]     = f[0*M2+indx_bot]+f[1*M2+indx_bot]+f[3*M2+indx_bot]+2.0*(f[4*M2+indx_bot]+f[7*M2+indx_bot]+f[8*M2+indx_bot]);
      f[2*M2+indx_bot]  = f[4*M2+indx_bot];
      f[5*M2+indx_bot]  = f[7*M2+indx_bot] - 0.5*(f[1*M2+indx_bot] - f[3*M2+indx_bot]);
      f[6*M2+indx_bot]  = f[8*M2+indx_bot] + 0.5*(f[1*M2+indx_bot] - f[3*M2+indx_bot]);
    }
    for(i = start+1; i <= end; i++) {
    // left wall - inlet - eastside
      indx_left    = My*i + 1;
      u_wall       = rho_in - (f[0*M2+indx_left]+f[2*M2+indx_left]+f[4*M2+indx_left]+2.*(f[3*M2+indx_left]+f[7*M2+indx_left]+f[6*M2+indx_left]));
      u[indx_left] = u_wall;
      v[indx_left] = 0.0;
      rho[indx_left]    = rho_in;
      f[1*M2+indx_left] = f[3*M2+indx_left] + two_by_three*u_wall;//
      f[5*M2+indx_left] = f[7*M2+indx_left] - 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;
      f[8*M2+indx_left] = f[6*M2+indx_left] + 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;//

    // right wall - outlet - westside
      indx_right   = My*i + My - 2;
      u_wall       = -1.0*rho_out + (f[0*M2+indx_right]+f[2*M2+indx_right]+f[4*M2+indx_right]+2.*(f[1*M2+indx_right]+f[5*M2+indx_right]+f[8*M2+indx_right]));
      u[indx_right] = u_wall;
      v[indx_right] = 0.0;
      rho[indx_right]    = rho_out;
      f[3*M2+indx_right] = f[1*M2+indx_right] - two_by_three*u_wall;
      f[7*M2+indx_right] = f[5*M2+indx_right] + 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;
      f[6*M2+indx_right] = f[8*M2+indx_right] - 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;
    }
    // bottom left corner
    indx_bot          = 1 + My;
    f[1*M2+indx_bot]  = f[3*M2+indx_bot];
    f[2*M2+indx_bot]  = f[4*M2+indx_bot];
    f[5*M2+indx_bot]  = f[7*M2+indx_bot];
    f[6*M2+indx_bot]  = 0.5*(rho_in - (f[0*M2+indx_bot]+f[1*M2+indx_bot]+f[2*M2+indx_bot]+f[3*M2+indx_bot]+f[4*M2+indx_bot]+f[5*M2+indx_bot]+f[7*M2+indx_bot]));
    f[8*M2+indx_bot]  = f[6*M2+indx_bot];
    // Bottom right corner
    indx_bot          = My - 2 + My;
    f[3*M2+indx_bot]  = f[1*M2+indx_bot];
    f[2*M2+indx_bot]  = f[4*M2+indx_bot];
    f[6*M2+indx_bot]  = f[8*M2+indx_bot];
    f[5*M2+indx_bot]  = 0.5*(rho_out - (f[0*M2+indx_bot]+f[1*M2+indx_bot]+f[2*M2+indx_bot]+f[3*M2+indx_bot]+f[4*M2+indx_bot]+f[6*M2+indx_bot]+f[8*M2+indx_bot]));
    f[7*M2+indx_bot]  = f[5*M2+indx_bot];
  }
  else if(taskid == numworkers){
    for (i = 1; i < My - 2; i++) {
    // top wall - north side - no slip
      indx_top         = i + (Mx-2)*My;

      u[indx_top] = 0.0;
      v[indx_top] = 0.0;
      rho[indx_top]     = f[0*M2+indx_top]+f[1*M2+indx_top]+f[3*M2+indx_top]+2.0*(f[2*M2+indx_top]+f[5*M2+indx_top]+f[6*M2+indx_top]);

      f[4*M2+indx_top] = f[2*M2+indx_top];
      f[7*M2+indx_top] = f[5*M2+indx_top] + 0.5*(f[1*M2+indx_top] - f[3*M2+indx_top]);
      f[8*M2+indx_top] = f[6*M2+indx_top] - 0.5*(f[1*M2+indx_top] - f[3*M2+indx_top]);
    }
    for(i = start; i <= (end-1); i++) {
    // left wall - inlet - eastside
      indx_left    = My*i + 1;
      u_wall       = rho_in - (f[0*M2+indx_left]+f[2*M2+indx_left]+f[4*M2+indx_left]+2.*(f[3*M2+indx_left]+f[7*M2+indx_left]+f[6*M2+indx_left]));
      u[indx_left] = u_wall;
      v[indx_left] = 0.0;
      rho[indx_left]    = rho_in;
      f[1*M2+indx_left] = f[3*M2+indx_left] + two_by_three*u_wall;//
      f[5*M2+indx_left] = f[7*M2+indx_left] - 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;
      f[8*M2+indx_left] = f[6*M2+indx_left] + 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;//

    // right wall - outlet - westside
      indx_right   = My*i + My - 2;
      u_wall       = -1.0*rho_out + (f[0*M2+indx_right]+f[2*M2+indx_right]+f[4*M2+indx_right]+2.*(f[1*M2+indx_right]+f[5*M2+indx_right]+f[8*M2+indx_right]));
      u[indx_right] = u_wall;
      v[indx_right] = 0.0;
      rho[indx_right]    = rho_out;
      f[3*M2+indx_right] = f[1*M2+indx_right] - two_by_three*u_wall;
      f[7*M2+indx_right] = f[5*M2+indx_right] + 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;
      f[6*M2+indx_right] = f[8*M2+indx_right] - 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;
    }
    // top left corner
    indx_top         = 1 + (Mx-2)*My;
    u[indx_top] = 0.0;
    v[indx_top] = 0.0;
    f[1*M2+indx_top] = f[3*M2+indx_top];
    f[4*M2+indx_top] = f[2*M2+indx_top];
    f[8*M2+indx_top] = f[6*M2+indx_top];
    rho[indx_top]     = f[0*M2+indx_top]+f[1*M2+indx_top]+f[3*M2+indx_top]+2.0*(f[2*M2+indx_top]+f[5*M2+indx_top]+f[6*M2+indx_top]);
    f[5*M2+indx_top] = 0.5*(rho_in - (f[0*M2+indx_top]+f[1*M2+indx_top]+f[2*M2+indx_top]+f[3*M2+indx_top]+f[4*M2+indx_top]+f[6*M2+indx_top]+f[8*M2+indx_top]));
    f[7*M2+indx_top] = f[5*M2+indx_top];

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
  }else{
    for(i = start; i <= end; i++) {
    // left wall - inlet - eastside
      indx_left    = My*i + 1;
      u_wall       = rho_in - (f[0*M2+indx_left]+f[2*M2+indx_left]+f[4*M2+indx_left]+2.*(f[3*M2+indx_left]+f[7*M2+indx_left]+f[6*M2+indx_left]));
      u[indx_left] = u_wall;
      v[indx_left] = 0.0;
      rho[indx_left]    = rho_in;
      f[1*M2+indx_left] = f[3*M2+indx_left] + two_by_three*u_wall;//
      f[5*M2+indx_left] = f[7*M2+indx_left] - 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;
      f[8*M2+indx_left] = f[6*M2+indx_left] + 0.5*(f[2*M2+indx_left] - f[4*M2+indx_left]) + one_by_six*u_wall;//

    // right wall - outlet - westside
      indx_right   = My*i + My - 2;
      u_wall       = -1.0*rho_out + (f[0*M2+indx_right]+f[2*M2+indx_right]+f[4*M2+indx_right]+2.*(f[1*M2+indx_right]+f[5*M2+indx_right]+f[8*M2+indx_right]));
      u[indx_right] = u_wall;
      v[indx_right] = 0.0;
      rho[indx_right]    = rho_out;
      f[3*M2+indx_right] = f[1*M2+indx_right] - two_by_three*u_wall;
      f[7*M2+indx_right] = f[5*M2+indx_right] + 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;
      f[6*M2+indx_right] = f[8*M2+indx_right] - 0.5*(f[2*M2+indx_right] - f[4*M2+indx_right]) - one_by_six*u_wall;
    }
  }
}
double calc_rho(int indx){
  int i, k, j, indx_f;
  double tmp = 0.0;
    for(k = 0; k < Q; k++ ) {
      indx_f = k*M2 + indx;
      tmp   += f[indx_f];
    }
    return tmp;
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
