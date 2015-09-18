void master_allocate_memory(int Mx, int My) {
  if (taskid == MASTER) {
    phi_old     =   (double *)malloc(Mx*My*sizeof(double));
    mu_old      =   (double *)malloc(Mx*My*sizeof(double));
  }
}
void worker_allocate_memory(int Mx, int My) {
  if ( taskid == MASTER ) {
    averow    =   My/numworkers;
    extra     =   My%numworkers;
    for ( rank=1; rank <= (numworkers); rank++) {
      rows      =   (rank <= extra) ? averow+1 : averow;
      dest      =   rank;
      MPI_Send(&rows,                 1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
    }
  }
  if(taskid != MASTER) {
    source =  MASTER;
    MPI_Recv(&rows,          1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
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
  }
}
mpi_distribute(int My){
  if ( taskid == MASTER ) {
    averow    =   My/numworkers;
    extra     =   My%numworkers;
    offset    =   0;
    for ( rank=1; rank <= (numworkers); rank++) {
      rows         =   (rank <= extra) ? averow+1 : averow;
      left_node    =   rank - 1;
      right_node   =   rank + 1;

      if ( rank == 1 ) {
        left_node  = NONE;
      }
      if ( rank == (numworkers) ) {
        right_node = NONE;
      }

      dest = rank;

      MPI_Send(&offset,               1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&rows,                 1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&left_node,            1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&right_node,           1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&phi_old[offset*My],      rows*My,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&mu_old[offset*My],       rows*My,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      offset = offset + rows;
    }
  }else{
    source =  MASTER;
    MPI_Recv(&offset,        1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
    MPI_Recv(&rows,          1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
    MPI_Recv(&left_node,     1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
    MPI_Recv(&right_node,    1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);

    start = 1;
    if((taskid ==1) || (taskid == numworkers)) {
      if(taskid == 1) {
        MPI_Recv(&phi_old[0],   rows*My,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&mu_old[0],    rows*My,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      }
      else {
        MPI_Recv(&phi_old[My],  rows*My,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&mu_old[My],   rows*My,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      }
      end = rows-1;
    } else {
      MPI_Recv(&phi_old[My],    rows*My,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&mu_old[My],     rows*My,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      end = rows;
    }
  }
}
void phasefield_initialize() {
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
void boundary_mpi(int taskid, double *c){
  int i ,y ,z;
  int indx_up, indx_dwn, indx_lft, indx_rght;
  if ( (taskid == 1) || (taskid == numworkers) ) {
    for (i = 0; i < MESHX; i++ ) {
      if ( taskid == 1 ){
        c[i]        = c[MESHX + i];
      }
      else if (taskid == numworkers) {
        indx    = rows*MESHX - MESHX   + i;
        indx_up = rwos*MESHX - 2*MESHX + i;
        c[indx]  = c[indx_up];
      }
    }
  }
  else{
    for (i=start; i <= end; i++){
      indx_rght     = i*MESHX;
      indx_lft      = i*MESHX + MESHX - 1;
      c[indx_lft]   = c[indx_lft + 1];
      c[indx_rght]  = c[indx_rght - 1];
    }
  }
}

void solverloop(int start, int end){

  int       i, j;
  int       indx, indx_up, indx_lft;
  int       indx_rght, indx_up, indx_dwn;
  double    phi,dphi_dt,dmu_dt;
  double    drv_frce, alln_chn;
  double    Gamma, kai;
  double    dc_dx, dc_dy, V_gradC = 0.0;

  // #ifdef ANISO
  //   grad_phi(1, dphi_now);
  // #endif

  for (i=start; i <= end; i++) {

    // #ifdef ANISO
    //   grad_phi(i+1, dphi_next);
    // #endif

    for (j=1; j < (MESHX-1); j++){

      indx          =   i*MESHX + j;
      indx_lft      =   indx - 1;
      indx_rght     =   indx + 1;
      indx_up       =   indx - MESHX;
      indx_dwn      =   indx + MESHX;

      phi           =   phi_old[indx];

      // #ifdef ISO
        Gamma       = 2*G*lap_phi[indx];
      // #endif
      // #ifdef ANISO
      //   Gamma       =     div_phi(j);
      // #endif

      drv_frce      =     (mu_old[indx] - Mu)*(K-1)*(mu_old[indx])*6*phi*(1-phi);
      alln_chn      =     E*Gamma - (G/E)*18.0*(phi)*(1.0-phi)*(1.0-2.0*phi);
      dp_dt         =     (alln_chn + drv_frce)/(tau*E);

      phi_new[z]    =     phi + deltat*dphi_dt;

      // #ifdef ANISO
      //   dc_dx       =     (conc[indx_rght] - conc[indx_lft])*0.5*inv_deltax;
      //   dc_dy       =     (conc[indx_dwn]  - conc[indx_up] )*0.5*inv_deltax;
      //   V_gradC     =     u_old[indx]*dc_dx + v_old[indx]*dc_dy;
      // #endif
      // #ifdef ISO
        // V_gradC     =     0.0;
      // #endif
      dmu_dt        =     Mob*lap_mu[indx] - V_gradC - (K-1)*mu_old[indx]*6*phi*(1-phi)*dphi_dt;
      kai           =     1+(K-1)*phi*phi*(3-2*phi);
      mu_new[indx]  =     mu_old[indx]  + deltat*dmu_dt/kai;
    }
    // #ifdef ANISO
    //   fnupdate();
    // #endif
  }
}
void mpiexchange(int taskid, double *P, int Mx) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      MPI_Send(&P[end*Mx],      Mx, MPI_DOUBLE,  right_node, LTAG,      MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*Mx],  Mx, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
    }
    MPI_Send(&P[start*Mx],      Mx, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    source  = left_node;
    msgtype = LTAG;
    MPI_Recv(&P[0],             Mx, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG;
       MPI_Recv(&P[0],          Mx, MPI_DOUBLE,   source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&P[start*Mx],   Mx, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*Mx],  Mx, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&P[(end)*Mx],    Mx, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
    }
  }
}
void sendtomaster(int taskid, double *c) {
  dest = MASTER;

  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&c[0],     rows*MESHX, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&c[MESHX], rows*MESHX, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  }
}
void receivefrmworker() {
  int rank;
  for (rank=1; rank <= numworkers; rank++) {
    source = rank;
    MPI_Recv(&offset,             1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,               1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,          1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node,         1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&phi_old[offset*MESHX],    rows*MESHX,    MPI_DOUBLE,    source,   WRITE,  MPI_COMM_WORLD, &status);
  }
}
