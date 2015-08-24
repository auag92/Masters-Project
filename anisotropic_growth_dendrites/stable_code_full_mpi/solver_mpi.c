void solverloop(int start, int end){

  int       i, j;
  int       indx, indx_up, indx_lft;
  int       indx_rght, indx_up, indx_dwn;
  double    phi,dphi_dt,dmu_dt;
  double    drv_frce, alln_chn;
  double    Gamma, kai;
  double    dc_dx, dc_dy, V_gradC;

  #ifdef ANISO
    grad_phi(1, dphi_now);
  #endif

  for (i=start; i < end; i++) {

    #ifdef ANISO
      grad_phi(i+1, dphi_next);
    #endif

    for (j=1; j < (MESHX-1); j++){

      indx          =   i*MESHX + j;
      indx_lft      =   indx - 1;
      indx_rght     =   indx + 1;
      indx_up       =   indx - MESHX;
      indx_dwn      =   indx + MESHX;

      phi           =   phi_old[indx];

      #ifdef ISO
        Gamma       = 2*G*lap_phi[indx];
      #endif
      #ifdef ANISO
        Gamma       =     div_phi(j);
      #endif

      drv_frce      =     (mu_old[indx] - Mu)*(K-1)*(mu_old[indx])*6*phi*(1-phi);
      alln_chn      =     E*Gamma - (G/E)*18.0*(phi)*(1.0-phi)*(1.0-2.0*phi);
      dp_dt         =     (alln_chn + drv_frce)/(tau*E);

      phi_new[z]    =     phi + deltat*dphi_dt;

      #ifdef ANISO
        dc_dx       =     (conc[indx_rght] - conc[indx_lft])*0.5*inv_deltax;
        dc_dy       =     (conc[indx_dwn]  - conc[indx_up] )*0.5*inv_deltax;
        V_gradC     =     u_old[indx]*dc_dx + v_old[indx]*dc_dy;
      #endif
      #ifdef ISO
        V_gradC     =     0.0;
      #endif
      dmu_dt        =     Mob*lap_mu[indx] - V_gradC - (K-1)*mu_old[indx]*6*phi*(1-phi)*dphi_dt;
      kai           =     1+(K-1)*phi*phi*(3-2*phi);
      mu_new[indx]  =     mu_old[indx]  + deltat*dmu_dt/kai;
    }
    #ifdef ANISO
      fnupdate();
    #endif
  }
}
void boundary_pressure_mpi(int taskid){
  int i ,y ,z;
  int indx_up, indx_dwn, indx_lft, indx_rght;
  if ( (taskid == 1) || (taskid == numworkers) ) {
    for (i = 0; i < pmesh; i++ ) {
      if ( taskid == 1 ){
        indx_up       = i;
        P[indx_up]    = p_up;
      }
      else if (taskid == numworkers) {
        indx_dwn      = end + i;
        P[indx_dwn]   = p_down;
      }
    }
  }
  else{
    for (i=start; i <= end; i++){
      indx_rght     = i*pmesh;
      indx_lft      = i*pmesh + pmesh - 1;
      P[indx_lft]   = p_left;
      P[indx_rght]  = p_right;
    }
  }
}
void neuman_boundary(double *c, int m) {
  int i ,y ,z;
  int m2 = m*m;
  for (i=0; i<m -1; i++)
  {
    y= i*m;
    z= i*m + m-1;
    //left - right
    c[y]       = c[y+1];
    c[z]       = c[z-1];
    // up - down
    c[i]       = c[MESHX + i];
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
