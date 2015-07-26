#define MASTER 0
#define NONE 0
#define BEGIN 999
#define LTAG  777
#define RTAG  666
#define WRITE 555
#define ERROR 888
#define BREAK 111
//------------------------------------------------
int numtasks, numworkers, taskid, rank, dest;
int averow, extra, offset;
int left_node, right_node;
int start, end;
int source, msgtype;
int rows;
MPI_Status status;
int flag = 0;
double error, err;
int iter = 0;
//-------------------------------------------------
void    boundary_pressure_mpi(int taskid);
void    red_solver(double *P, double *fn, int start, int end);
void    black_solver(double *P, double *fn, int start, int end);
double  compute_error_mpi(double *P, double *fn, int start, int end);
void    mpiexchange(int taskid);
void    sendtomaster(int taskid);
void    receivefrmworker();
// gs_mpi( double *P, double *fn ,double *a_x, double *a_y ){
void gs_mpi() {
  if ( taskid == MASTER ) {
    averow    =   pmesh/numworkers;
    extra     =   pmesh%numworkers;
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
      MPI_Send(&P[offset*pmesh],      rows*pmesh,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&rhs_fn[offset*pmesh], rows*pmesh,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);

      offset = offset + rows;
    }
    iter = 0;
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
      if ( error < gs_tol ) {
        receivefrmworker();
        break;
      }
    }
    printf("I have finished my iterations\n");
  }
  if(taskid != MASTER) {
    // printf("hello world from worker no. %d\n", taskid);
    source =  MASTER;
    MPI_Recv(&offset,        1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
    MPI_Recv(&rows,          1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
    MPI_Recv(&left_node,     1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);
    MPI_Recv(&right_node,    1,      MPI_INT,     source,    BEGIN,   MPI_COMM_WORLD,  &status);

    start = 1;
    if((taskid ==1) || (taskid == numworkers)) {
      P            =   (double *)malloc((rows+1)*pmesh*sizeof(double));
      rhs_fn       =   (double *)malloc((rows+1)*pmesh*sizeof(double));
      if(taskid == 1) {
        MPI_Recv(&P[0],      rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&rhs_fn[0],     rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      }
      else {
        MPI_Recv(&P[pmesh],      rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&rhs_fn[pmesh],     rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      }
      end = rows-1;
    } else {
      P            =   (double *)malloc((rows+2)*pmesh*sizeof(double));
      rhs_fn           =   (double *)malloc((rows+2)*pmesh*sizeof(double));
      MPI_Recv(&P[pmesh],          rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&rhs_fn[pmesh],     rows*pmesh,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      end = rows;
    }

    for (; ;) {
      //boundary_pressure_mpi(taskid);
      mpiexchange(taskid);
      // printf("I have recieved my share of work from master: %d\n", taskid);
      red_solver(P, rhs_fn, start, end);
      mpiexchange(taskid);
      black_solver(P, rhs_fn, start, end);
      error = compute_error_mpi(P, rhs_fn, start, end);

      MPI_Send(&error,     1,     MPI_DOUBLE,    MASTER,   ERROR,  MPI_COMM_WORLD);
      MPI_Recv(&flag,      1,     MPI_INT,       MASTER,   BREAK,  MPI_COMM_WORLD,  &status);

      if (flag == 1) {
        sendtomaster(taskid);
        break;
      }
    }
    free(P);
    free(rhs_fn);
  }
}
void red_solver(double *P, double *fn, int start, int end){
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
        P[indx]  = -1.0*(P[indx_lft] + P[indx_rght]
        + P[indx_bck] + P[indx_frnt]) + deltax*deltax*fn[indx];
        P[indx] /= -1.0*4.0;
      }
      // if (((i+j)%2) == 0) {
      //   P[indx]  = -1.0*(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
      //   + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]) + deltax*deltax*fn[indx];
      //   P[indx] /= -1.0*(a_x[indx_ax+1]+a_x[indx_ax]+a_y[indx_ay+(pmesh-2)]+a_y[indx_ay]);
      // }
    }
  }
}
void black_solver(double *P, double *fn, int start, int end){
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
        P[indx]  = -1.0*(P[indx_lft] + P[indx_rght]
        + P[indx_bck] + P[indx_frnt]) + deltax*deltax*fn[indx];
        P[indx] /= -1.0*4.0;
      }

      // if (((i+j)%2) != 0) {
      //   P[indx]  =-1.0*(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
      //   + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]) + deltax*deltax*fn[indx];
      //   P[indx] /= -1.0*(a_x[indx_ax]+a_x[indx_ax+1]+a_y[indx_ay]+a_y[indx_ay+(pmesh-2)]);
      // }
    }
  }
}
double compute_error_mpi(double *P, double *fn, int start, int end ) {
  double error=0.0;
  int indx, indx_rght, indx_frnt, indx_lft, indx_bck, i, j,indx_ax,indx_ay;
  for(i=start; i <= end; i++) {
    for(j=1;j < pmesh-1; j++) {
      indx  = i*pmesh    + j;
      indx_frnt         = indx      + pmesh;
      indx_rght         = indx      + 1;
      indx_lft          = indx      - 1;
      indx_bck          = indx      - pmesh;
      indx_ax      = (i-1)*(pmesh-1)  + j-1;
      indx_ay      = (i-1)*(pmesh-2)  + j-1;

      error += fabs(P[indx_lft] + P[indx_rght]
      + P[indx_bck] + P[indx_frnt] - (4.0)*P[indx] - deltax*deltax*fn[indx]);
      // error += fabs(a_x[indx_ax]*P[indx_lft] + a_x[indx_ax+1]*P[indx_rght]
      // + a_y[indx_ay]*P[indx_bck] + a_y[indx_ay+(pmesh-2)]*P[indx_frnt]
      // - (a_x[indx_ax]+a_x[indx_ax+1]+a_y[indx_ay]+a_y[indx_ay+(pmesh-2)])*P[indx]
      // - deltax*deltax*fn[indx]);
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
    // MPI_Recv(&P[0],                pmesh, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
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
      // MPI_Send(&P[(end)*pmesh],    pmesh, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
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
