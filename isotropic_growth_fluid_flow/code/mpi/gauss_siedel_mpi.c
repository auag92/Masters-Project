#define MASTER 0
#define NONE 0
#define BEGIN 999
#define LTAG  777
#define RTAG  666
#define WRITE 555

gs_mpi(double *P, double *fn, double *a_x, double *a_y){
  double error;
  double tol=1.0e-6;
  int iter  = 0;
  if ( taskid == MASTER ){
    averow    =   MESHX/numworkers;
    extra     =   MESHX%numworkers;
    offset    =   0;
    for (rank=1; rank <= (numworkers); rank++){
      rows        =    (rank <= extra) ? averow+1 : averow;

      left_node    =   rank - 1;
      right_node   =   rank + 1;

      if ( rank == 1 ){
        left_node  = NONE;
      }
      if ( rank == (numworkers) ){
        right_node = NONE;
      }

      dest = rank;
      MPI_Send(&offset,         1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&rows,           1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&left_node,      1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&right_node,     1,                   MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&P[offset],      rows*pmesh,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&fn[offset],     rows*pmesh,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&a_x[offset],    (rows-2)*(pmesh-1),  MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&a_y[offset],    (rows-1)*(pmesh-2),  MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);


      offset = offset + rows;
    }
    for (t=1; t < ntimesteps; t++) {
      if (t%saveT == 0) {
        receivefrmworker();
        check_error();
      }
    }
  }
  if(taskid != MASTER){

    MPI_Recv(&offset,         1,                   MPI_INT,         source,   BEGIN,  MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,           1,                   MPI_INT,         source,   BEGIN,  MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,      1,                   MPI_INT,         source,   BEGIN,  MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node,     1,                   MPI_INT,         source,   BEGIN,  MPI_COMM_WORLD, &status);

    P            =   (double *)malloc(rows*pmesh*sizeof(double));
    fn           =   (double *)malloc(rows*pmesh*sizeof(double));
    a_x          =   (double *)malloc(rows*pmesh*sizeof(double));
    a_y          =   (double *)malloc(rows*pmesh*sizeof(double));

    if(taskid ==1){
      MPI_Recv(&P[0],      rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&fn[0],     rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&a_x[0],    (rows-2)*(pmesh-1),  MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&a_y[0],    (rows-1)*(pmesh-2),  MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
    }
    else if(taskid == numworkers){
      MPI_Recv(&P[1],      rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&fn[1],     rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&a_x[1],    (rows-2)*(pmesh-1),  MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&a_y[1],    (rows-1)*(pmesh-2),  MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&P[0],      rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&fn[0],     rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&a_x[0],    (rows-2)*(pmesh-1),  MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&a_y[0],    (rows-1)*(pmesh-2),  MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
    }

    boundary_pressure(taskid);
    red_solver();
    mpiexchange();
    black_solver();
    sendtomaster();
  }
}
void red_solver(){
  int i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck, indx_ax, indx_ay;
  for(i=1; i < rows-1; i++) {
    for(j=1;j < pmesh-1;j++) {
      indx         = i*pmesh      + j;
      indx_rght    = indx         + 1;
      indx_frnt    = indx         + pmesh;
      indx_lft     = indx         - 1;
      indx_bck     = indx         - pmesh;
      indx_ax      = (i-1)*(pmesh-1)  +(j-1);
      indx_ay      = (i-1)*(pmesh-2)  +(j-1);

      if (((i+j)%2) == 0) {
        P[indx]  = -1.0*(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
        + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]) + deltax*deltax*fn[indx];
        P[indx] /= -1.0*(a_x[indx_ax+1]+a_x[indx_ax]+a_y[indx_ay+(pmesh-2)]+a_y[indx_ay]);
      }
    }
  }
}
void black_solver(int rows){
  int i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck, indx_ax, indx_ay;
  for(i=1; i < rows-1; i++) {
    for(j=1;j < pmesh-1;j++) {
      indx         = i*pmesh      + j;
      indx_rght    = indx         + 1;
      indx_frnt    = indx         + pmesh;
      indx_lft     = indx         - 1;
      indx_bck     = indx         - pmesh;
      indx_ax      = (i-1)*(pmesh-1)  + (j-1);
      indx_ay      = (i-1)*(pmesh-2)  + (j-1);

      if (((i+j)%2) != 0) {
        P[indx]  =-1.0*(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
        + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]) + deltax*deltax*fn[indx];
        P[indx] /= -1.0*(a_x[indx_ax]+a_x[indx_ax+1]+a_y[indx_ay]+a_y[indx_ay+(pmesh-2)]);
      }
    }
  }
}
double compute_error(double *a_x, double *a_y, double *P, double *fn) {
  double error=0.0;
  int indx, indx_rght, indx_frnt, indx_lft, indx_bck, i, j,indx_ax,indx_ay;
  for(i=1; i < pmesh-1; i++) {
    for(j=1;j < pmesh-1; j++) {
      indx  = i*pmesh    + j;
      indx_rght         = indx      + 1;
      indx_frnt         = indx      + pmesh;
      indx_lft          = indx      - 1;
      indx_bck          = indx      - pmesh;
      indx_ax      = (i-1)*(pmesh-1)  + j-1;
      indx_ay      = (i-1)*(pmesh-2)  + j-1;

      error += fabs(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
      + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]
      - (a_x[indx_ax]+a_x[indx_ax+1]+a_y[indx_ay]+a_y[indx_ay+(pmesh-2)])*P[indx]
      - deltax*deltax*fn[indx]);
    }
  }
  return(error);
}
void sendtomaster(int taskid) {
  dest = MASTER;

  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&P[0], rows*pmesh, MPI_DOUBLE,         dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&P[1], rows*pmesh, MPI_DOUBLE,         dest, WRITE, MPI_COMM_WORLD);
  }
}
void receivefrmworker() {
  int rank;
  for (rank=1; rank <= numworkers; rank++) {
    source = rank;
    MPI_Recv(&offset,       1,       MPI_INT,         source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,         1,       MPI_INT,         source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,    1,       MPI_INT,         source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node,   1,       MPI_INT,         source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&P[offset],    rows*pmesh,    MPI_DOUBLE,      source, WRITE, MPI_COMM_WORLD, &status);
  }
}
void mpiexchange(int taskid) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      MPI_Send(&P[end],    pmesh, MPI_DOUBLE,  right_node, LTAG,     MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[end+1],  pmesh, MPI_DOUBLE,  source,     msgtype,  MPI_COMM_WORLD, &status);
    }
    MPI_Send(&P[start],    pmesh, MPI_DOUBLE,  left_node,  RTAG,     MPI_COMM_WORLD);
    source  = left_node;
    msgtype = LTAG;
    MPI_Recv(&P[0],        pmesh, MPI_DOUBLE,  source,     msgtype,  MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG;
       MPI_Recv(&P[0],     pmesh, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&P[start], pmesh, MPI_DOUBLE, left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[end+1],  pmesh, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&P[end],    pmesh, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
    }
  }
}
void boundary_pressure(){

  int i ,y ,z;

  for (i=0; i< pmesh; i++)
  {
    y= i*pmesh;
    z= i*pmesh + pmesh-1;

    #ifdef neumann_pressure
    //left - right
    P[y]        = P[y+1];
    P[z]        = P[z-1];
    // up - down
    P[i]          = P[pmesh + i];
    P[pmesh2-pmesh+i]     = P[pmesh2-2*pmesh+i];
    #endif
    #ifdef dirichlet_pressure
    //left - right
    P[y]        = p_left;
    P[z]        = p_right;
    // up - down
    P[i]        = p_up;
    P[pmesh2-pmesh+i]  = p_down;
    #endif
  }
}
