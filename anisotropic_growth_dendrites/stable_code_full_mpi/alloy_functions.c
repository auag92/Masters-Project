void    allocate_memory_phasefields(int Mx);
void    initialize_phasefields();
void    mpi_distribute(int Mx);
void    boundary_neuman_mpi(int taskid, double *c, int Mx);
void    laplacian(double *f, double *lap, int Mx);
void    solverloop();
void    phi_update();
void    grad_phi(int i, double *d_phi);
double  dqdx( double phi_x, double phi_y);
double  div_phi(int i);
void    fnupdate();
void    write2file_phi ( int t, int m, double *c);
void    free_memory();


void allocate_memory_phasefields(int Mx) {
  if (taskid == MASTER) {
    phi_old     =   (double *)malloc(Mx*Mx*sizeof(double));
    mu_old      =   (double *)malloc(Mx*Mx*sizeof(double));
    averow    =   Mx/numworkers;
    extra     =   Mx%numworkers;
    for ( rank=1; rank <= (numworkers); rank++) {
      rows      =   (rank <= extra) ? averow+1 : averow;
      dest      =   rank;
      MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
    }
  }
  else{
    source =  MASTER;
    MPI_Recv(&rows, 1, MPI_INT, source, BEGIN, MPI_COMM_WORLD, &status);
    if((taskid ==1) || (taskid == numworkers)) {
      phi_old       =   (double *)malloc((rows+1)*Mx*sizeof(double));
      phi_new       =   (double *)malloc((rows+1)*Mx*sizeof(double));
      mu_old        =   (double *)malloc((rows+1)*Mx*sizeof(double));
      mu_new        =   (double *)malloc((rows+1)*Mx*sizeof(double));
    } else {
      phi_old       =   (double *)malloc((rows+2)*Mx*sizeof(double));
      phi_new       =   (double *)malloc((rows+2)*Mx*sizeof(double));
      mu_old        =   (double *)malloc((rows+2)*Mx*sizeof(double));
      mu_new        =   (double *)malloc((rows+2)*Mx*sizeof(double));
    }
    lap_phi        =   (double *)malloc((rows)*Mx*sizeof(double));
    lap_mu         =   (double *)malloc((rows)*Mx*sizeof(double));
    dphi_now       =   (double *)malloc(Mx*4*sizeof(double));
    dphi_next      =   (double *)malloc(Mx*4*sizeof(double));
  }
}
void free_memory(){
  if(taskid == MASTER){
    free(phi_old);
    free(mu_old);
  }else {
    free(phi_old);
    free(mu_old);
    free(mu_new);
    free(phi_new);
    free(lap_phi);
    free(lap_mu);
    free(dphi_now);
    free(dphi_next);
  }
}
void initialize_phasefields() {
  long i,j,z;
  double r;
#ifdef Centre
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      r= (i-MESHX*0.5)*(i-MESHX*0.5) + (j-MESHX*0.5)*(j-MESHX*0.5);
      z= i*MESHX + j;
      if(r < radius2){
	      phi_old[z] = 1.0;
      }
      else{
	      phi_old[z] = 0.0;
      }
      mu_old[z] = Mu - deltaMu;
    }
  }
#endif
#ifdef Corner
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      r= (i)*(i) + (j)*(j);
      z= i*MESHX + j;
      if(r < radius2){
	      phi_old[z] = 1.0;
      }
      else{
	      phi_old[z] = 0.0;
      }
      mu_old[z] = Mu - deltaMu;
    }
  }
#endif
#ifdef Nothing
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      z= i*MESHX + j;
      phi_old[z] = 0.0;
      mu_old[z] = 0.0;
    }
  }
#endif
}
void boundary_neuman_mpi(int taskid, double *c, int Mx){
  int i, indx;
  int indx_up, indx_dwn, indx_lft, indx_rght;
  if ( (taskid == 1) || (taskid == numworkers) ) {
    if ( taskid == 1 ){
      for (i = 0; i < Mx; i++ ) {
        c[i]        = c[Mx + i];
      }
    }
    else if (taskid == numworkers) {
      for (i = 0; i < Mx; i++ ) {
        indx_up    = (rows-1)*Mx + i;
        indx       = (rows)*Mx   + i;
        c[indx]    = c[indx_up];
      }
    }
  }
  for (i=start; i <= end; i++){
    indx_rght     = i*Mx;
    indx_lft      = i*Mx + Mx - 1;
    c[indx_lft]   = c[indx_lft + 1];
    c[indx_rght]  = c[indx_rght - 1];
  }
}
void laplacian(double *f, double *lap, int Mx) {
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
void phi_update() {
  long i, j, z;
  for (i=start; i <= end; i++) {
    for (j = 0; j < MESHX; j++){
      z= i*MESHX + j;
      phi_old[z]=phi_new[z];
      mu_old[z]=mu_new[z];
    }
  }
}
void solverloop(){

  int       i, j, z;
  double    p,dp_dt,dmu_dt;
  double    drv_frce, alln_chn;
  double    Gamma, kai;
  double    dc_dx, dc_dy, V_gradC = 0.0;

  #ifdef ANISO
  grad_phi(1, dphi_now);
  #endif

  for (i=start; i <= end; i++) {

    #ifdef ANISO
    grad_phi(i+1, dphi_next);
    #endif

    for (j=1; j < (MESHX-1); j++){

      z =   i*MESHX + j;
      p =   phi_old[z];

      #ifdef ISO
      Gamma = 2*G*lap_phi[z];
      #endif

      #ifdef ANISO
      Gamma         =     div_phi(j);
      #endif

      drv_frce      =     (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p);
      alln_chn      =     E*Gamma - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p);
      dp_dt         =     (alln_chn + drv_frce)/(tau*E);

      phi_new[z]    =     p + deltat*dp_dt;

      #ifdef FLUID
      dc_dx         =     (conc[z+1]-conc[z-1])*0.5*inv_deltax;
      dc_dy         =     (conc[z+MESHX]-conc[z-MESHX])*0.5*inv_deltax;
      V_gradC       =     u_old[z]*dc_dx + v_old[z]*dc_dy;
      #endif

      dmu_dt        =     Mob*lap_mu[z] - V_gradC - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
      kai           =     1+(K-1)*p*p*(3-2*p);
      mu_new[z]     =     mu_old[z]  + deltat*dmu_dt/kai;
    }
    #ifdef ANISO
    fnupdate();
    #endif
  }
}
void grad_phi(int i, double *d_phi){

 	int j,z;

  for(j=1; j<MESHX-1; j++){
    z = i * MESHX + j;
    if (  i == MESHX -1 ){
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
  if (  i != MESHX-1  ) {
    z = i * MESHX + j;
    d_phi[j]           = (phi_old[z] - phi_old[z-1])*inv_deltax;
    d_phi[MESHX+j]     = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
  }
}
double dqdx( double phi_x, double phi_y) {

	double   a, phi_x2, phi_x4, phi_y2, phi_y4, inv_phi;
	int      z;
  double   ans = 0;
  double   part1, part2, part3, part4;
  phi_x2    =   phi_x *phi_x;
  phi_y2    =   phi_y *phi_y;
  phi_y4    =   phi_y2 *phi_y2;
  phi_x4    =   phi_x2 *phi_x2;

  if ((phi_x2> 1e-15) && (phi_y2> 1e-15)){
    inv_phi    =    1/(phi_x2+phi_y2);
    part1      =    (1-Dab*(3-4*(phi_x4+phi_y4)*inv_phi*inv_phi));
    part2      =    2*G*E*part1*part1*phi_x;
    part3      =    32*G*E*Dab*(phi_x2+phi_y2)*(part1);
    part4      =    phi_x2*phi_x*inv_phi*inv_phi - phi_x*(phi_x4+phi_y4)*inv_phi*inv_phi*inv_phi;

    ans        =    part2 + part3*part4;
  }

  return ans;
}
double div_phi(int i){

  double    ans;
	double    x_next,    x_now;
	double    y_next,    y_now;

  x_now     = dqdx(dphi_now[i], dphi_now[i+MESHX]);
  x_next    = dqdx(dphi_now[i+1], dphi_now[i+1+MESHX]);
  y_now     = dqdx(dphi_now[i+2*MESHX], dphi_now[i+3*MESHX]);
  y_next    = dqdx(dphi_next[i+2*MESHX], dphi_next[i+3*MESHX]);
	ans       = ((x_next - x_now) + ( y_next - y_now))*inv_deltax;

  return ans;
}
void fnupdate()
{
  int i;

  for( i=0; i < MESHX; i++ ) {
    dphi_now[i]           =   dphi_next[i];
    dphi_now[MESHX+i]     =   dphi_next[MESHX+i];
    dphi_now[2*MESHX+i]   =   dphi_next[2*MESHX+i];
    dphi_now[3*MESHX+i]   =   dphi_next[3*MESHX+i];
  }
}
