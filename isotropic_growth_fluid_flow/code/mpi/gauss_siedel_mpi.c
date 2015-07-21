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

      MPI_Send(&offset,       1,           MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&rows,         1,           MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&left_node,    1,           MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&right_node,   1,           MPI_INT,         dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&P[offset],      rows*pmesh,  MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&fn[offset],     rows*pmesh,  MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
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
    boundary_pressure(taskid);
    red_solver();
    black_solver();
  }
}
void red_solver(int rows){
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
