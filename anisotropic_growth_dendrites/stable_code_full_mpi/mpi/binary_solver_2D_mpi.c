void      phi_update();
void      phi_initialize();
void      neuman_boundary(double *c, int m);
void      concentration(double *phi, double *mu, double *c, int m );
void      write2file_phi ( int t, int m,double *phi);
void      allocate_memory();
void      free_memory();
void      isotropic_solverloop();
void      anisotropic_solverloop();
void      grad_phi(int i, double *d_phi);
double    dqdx( double phi_x, double phi_y);
double    div_phi(int i);
void      fnupdate();
/******************************************************************************/
void allocate_memory() {
  phi_old     =   (double *)malloc(MESHX*MESHX*sizeof(double));
  phi_new     =   (double *)malloc(MESHX*MESHX*sizeof(double));
  mu_old      =   (double *)malloc(MESHX*MESHX*sizeof(double));
  mu_new      =   (double *)malloc(MESHX*MESHX*sizeof(double));
  lap_phi     =   (double *)malloc(MESHX*MESHX*sizeof(double));
  lap_mu      =   (double *)malloc(MESHX*MESHX*sizeof(double));
  conc        =   (double *)malloc(MESHX*MESHX*sizeof(double));
  dphi_now    =   (double *)malloc(MESHX*4*sizeof(double));
  dphi_next   =   (double *)malloc(MESHX*4*sizeof(double));
  P           =   (double *)malloc(pmesh*pmesh*sizeof(double));
  a_x         =   (double *)malloc((pmesh-2)*(pmesh-1)*sizeof(double));
  a_y         =   (double *)malloc((pmesh-1)*(pmesh-2)*sizeof(double));
  rhs_fn      =   (double *)malloc(pmesh*pmesh*sizeof(double));
  v_old       =   (double *)malloc(MESHX*MESHX*sizeof(double));
  u_old       =   (double *)malloc(MESHX*MESHX*sizeof(double));
  v_now       =   (double *)malloc(MESHX*MESHX*sizeof(double));
  u_now       =   (double *)malloc(MESHX*MESHX*sizeof(double));
  v_str       =   (double *)malloc(MESHX*MESHX*sizeof(double));
  u_str       =   (double *)malloc(MESHX*MESHX*sizeof(double));
  Hx          =   (double *)malloc(MESHX*MESHX*sizeof(double));
  Hy          =   (double *)malloc(MESHX*MESHX*sizeof(double));
}
void free_memory(){
  free(phi_old);
  free(phi_new);
  free(mu_old) ;
  free(mu_new);
  free(lap_phi);
  free(lap_mu);
  free(conc);
  free(P);
  free(v_old);
  free(u_old);
  free(v_now);
  free(u_now);
  free(v_str);
  free(u_str);
  free(a_x);
  free(a_y);
  free(rhs_fn);
  free(Hx);
  free(Hy);
  //===========
  free(dphi_now);
  free(dphi_next);
}
void phi_update() {
  long i, j, z;
  for (i=0; i < MESHX; i++) {
    for (j=0; j < MESHX; j++){
      z= i*MESHX + j;
      phi_old[z]=phi_new[z];
      mu_old[z]=mu_new[z];
    }
  }
}
void phi_initialize() {
  long i,j,z;
  double r;
#ifdef Centre
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      r= (i-pmesh*0.5)*(i-pmesh*0.5) + (j-pmesh*0.5)*(j-pmesh*0.5);
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
void neuman_boundary(double *c, int m) {
  int i ,y ,z;
  int m2 = m*m;
  for (i=0; i<m -1; i++)
  {
    y= i*m;
    z= i*m + m-1;
    //left - right
    c[y]        = c[y+1];
    c[z]  = c[z-1];
    // up - down
    c[i]        = c[MESHX + i];
    c[m2-m+i]  = c[m2-2*m+i];
  }
}
void concentration(double *phi, double *mu, double *c, int m ){
  double p,u,h;
  int i, j, z;
  for ( i = 0; i < m; i++)
  {
    for ( j = 0; j < m; j++){
      z= i*m + j;
      p = phi[z];
      u = mu[z];
      h = p*p*(3-2*p);
      c[z] = u*(1-h) + (h)*K*u;
    }
  }
}

void anisotropic_solverloop(){

  int       i, j, z;
  double    p,dp_dt,dmu_dt;
  double    drv_frce, alln_chn;
  double    Gamma, kai;
  double    dc_dx, dc_dy, V_gradC;

  #ifdef ANISO
  grad_phi(1, dphi_now);
  #endif

  for (i=1; i < (MESHX-1); i++) {

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

      dc_dx         =     (conc[z+1]-conc[z-1])*0.5*inv_deltax;
      dc_dy         =     (conc[z+MESHX]-conc[z-MESHX])*0.5*inv_deltax;
      V_gradC       =     u_old[z]*dc_dx + v_old[z]*dc_dy;

      dmu_dt        =     Mob*lap_mu[z] - V_gradC - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
      // dmu_dt        =     Mob*lap_mu[z] - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
      kai           =     1+(K-1)*p*p*(3-2*p);
      mu_new[z]     =     mu_old[z]  + deltat*dmu_dt/kai;
    }
    fnupdate();
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
