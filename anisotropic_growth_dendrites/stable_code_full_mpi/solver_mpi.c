void solverloop(){

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

void sendtomaster(int taskid) {
  dest = MASTER;

  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&P[0],     rows*pmesh, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&P[pmesh], rows*pmesh, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
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
    MPI_Recv(&P[offset*pmesh],    rows*pmesh,    MPI_DOUBLE,    source,   WRITE,  MPI_COMM_WORLD, &status);
  }
}
void mpiexchange(int taskid) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      MPI_Send(&P[end*pmesh],       pmesh, MPI_DOUBLE,  right_node, LTAG,      MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*pmesh],  pmesh, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
    }
    MPI_Send(&P[start*pmesh],      pmesh, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    source  = left_node;
    msgtype = LTAG;
    MPI_Recv(&P[0],                pmesh, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG;
       MPI_Recv(&P[0],             pmesh, MPI_DOUBLE,   source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&P[start*pmesh],   pmesh, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*pmesh],  pmesh, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&P[(end)*pmesh],    pmesh, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
    }
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
void allocate_memory(){
  if ( taskid == MASTER ) {
    averow    =   pmesh/numworkers;
    extra     =   pmesh%numworkers;

    rank      =    1;
    rows      =   (rank <= extra) ? averow+1 : averow;
    dest      =   rank;
    rows++;
    MPI_Send(&rows,       1,      MPI_INT,      dest,     BEGIN,     MPI_COMM_WORLD);
    rows--;

    for ( rank=2; rank <= (numworkers); rank++) {
      rows         =   (rank <= extra) ? averow+1 : averow;
      dest = rank;
      MPI_Send(&rows,      1,     MPI_INT,       dest,   BEGIN,  MPI_COMM_WORLD);
    }
  }
  if(taskid != MASTER) {
    source =  MASTER;
    MPI_Recv(&rows,          1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
    if((taskid ==1) || (taskid == numworkers)) {
      phi_old     =   (double *)malloc(rows*MESHX*sizeof(double));
      phi_new     =   (double *)malloc(rows*MESHX*sizeof(double));
      mu_old      =   (double *)malloc(rows*MESHX*sizeof(double));
      mu_new      =   (double *)malloc(rows*MESHX*sizeof(double));
      lap_phi     =   (double *)malloc(rows*MESHX*sizeof(double));
      lap_mu      =   (double *)malloc(rows*MESHX*sizeof(double));
      conc        =   (double *)malloc(rows*MESHX*sizeof(double));
    } else {
      phi_old     =   (double *)malloc(rows*MESHX*sizeof(double));
      phi_new     =   (double *)malloc(rows*MESHX*sizeof(double));
      mu_old      =   (double *)malloc(rows*MESHX*sizeof(double));
      mu_new      =   (double *)malloc(rows*MESHX*sizeof(double));
      lap_phi     =   (double *)malloc(rows*MESHX*sizeof(double));
      lap_mu      =   (double *)malloc(rows*MESHX*sizeof(double));
      conc        =   (double *)malloc(rows*MESHX*sizeof(double));
    }
    dphi_now    =   (double *)malloc(MESHX*4*sizeof(double));
    dphi_next   =   (double *)malloc(MESHX*4*sizeof(double));
  }
}
