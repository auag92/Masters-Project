#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mpi.h"
//-----------------------------------------------------------//
#include "LBM_constants.h"
// #include "LBM_variables.h"
// #include "LBM_functions.h"
#include "MPI_def.h"
//------------------------------------------------------------//
/*Variable definitions*/
//------------------------------------------------------------//
double *f, *f_str;
double feq[Q];
double *u, *v, *rho;
double w[Q]    = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
double e[Q][2] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0};
int t;
/*Function Definitions*/
void allocate_rows();
void allocate_memory();
void free_memory();

void    calculate_rho();
void    calculate_velocities();
void    Fk_eq( double u, double v, double r );
void    collision_step();
void    streaming_step();
void    boundary_pipeflow();
void    init();

void mpiexchange(int taskid, double *c, int M);
void sendtomaster(int taskid, double *c, int M);
void receivefrmworker(double *c, int M);

void write2file( int t, int m, double *c, char *fname);
void write2file1( int t, int m, double *c, double *d, char *fname);
/*main*/
void main(int argc, char *argv[]) {
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  numworkers = numtasks - 1;

  allocate_rows();
  allocate_memory();

  char fname1[100], fname2[100];
  int  m = sprintf(fname1,"velocities");
  int  n = sprintf(fname2,"rho");

  if (taskid == MASTER){
    printf("Hi, I am the Master node.\n");
    for (t = 1; t <= tsteps; t++) {
      if(t%savet == 0){
        receivefrmworker(u, Mx);
        receivefrmworker(v, Mx);
        printf("%d\n",t);
        write2file1(t,Mx,u,v,fname1);
        // write2file(t,Mx,rho,fname2);
      }
    }
  }else {
    printf("Hi, I am worker node no. %d, with no. of rows = %d. start = %d, end = %d\n", taskid, rows, start, end);
    init();
    for (t = 1; t <= tsteps; t++) {
      collision_step();
      streaming_step();
      // mpiexchange(taskid, f, Q*Mx);
      // boundary_pipeflow();
      calculate_rho();
      calculate_velocities();
      if(t%savet == 0){
        sendtomaster(taskid, u, Mx);
        sendtomaster(taskid, v, Mx);
      }
    }
  }
  free_memory();
  MPI_Finalize();
}

void allocate_rows(){
  if (taskid == MASTER) {
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
  }
}
void allocate_memory() {

  if(taskid ==  MASTER) {
  //velocity matrices
  u      = (double *)malloc(Mx*My*sizeof(double));
  v      = (double *)malloc(Mx*My*sizeof(double));
  rho    = (double *)malloc(Mx*My*sizeof(double));

  }else {
    //velocity matrices
    u      = (double *)malloc((rows+2)*My*sizeof(double));
    v      = (double *)malloc((rows+2)*My*sizeof(double));
    rho    = (double *)malloc((rows+2)*My*sizeof(double));
    start = 1;
    // probability distribution functions
    if(taskid == 1 || taskid == numworkers){
      f     = (double *)malloc(Q*(rows+1)*My*sizeof(double));
      f_str = (double *)malloc(Q*(rows+1)*My*sizeof(double));
      end   = rows-1;
    }else{
      f     = (double *)malloc(Q*(rows+2)*My*sizeof(double));
      f_str = (double *)malloc(Q*(rows+2)*My*sizeof(double));
      end   = rows;
    }
  }
}
void free_memory() {
  free(u);
  free(v);
  free(rho);
  free(f);
  free(f_str);
}
void calculate_rho(){
  int i,j,k,indx, indx_f;
  double tmp = 0.0;
  for(i = start; i <= end ; i++) {
    for(j = 1; j < My-1 ; j++) {
      tmp  = 0.0;
      indx = i*Mx + j;
      for(k = 0; k < Q; k++ ){
        indx_f  = i*Q*My + k*My +j;
        tmp    += f[indx_f];
      }
      rho[indx] = tmp;
    }
  }
}
void calculate_velocities() {
  int i, j, k, indx_f, indx;
  double tmp1, tmp2, inv_rho;
  for(i = start; i <= end; i++) {
    for(j = 1; j < My-1; j++) {
      indx = i*Mx + j;
      tmp1 = 0.0;
      tmp2 = 0.0;
      for(k = 0; k < Q; k++ ){
        indx_f  = i*Q*My + k*My +j;
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
  for(i = start; i <= end ; i++) {
    for(j = 1; j < My-1 ; j++) {
      indx   = i*My + j;
      Fk_eq( u[indx], v[indx], rho[indx] );
      for( k = 0; k < Q; k++ ){
        indx_f        = i*Q*My + k*My +j;
        f_str[indx_f] = omega*feq[k] + (1.0 - omega)*f[indx_f];
      }
    }
  }
}
void streaming_step() {
  int i, j, k, indx, indx_f, indx_next;
  for(i = start; i <= end; i++) {
    for(j = 1; j < My - 1; j++) {
      indx = i*My + j;
      for(k = 0; k < Q; k++ ){
        indx_f       = i*Q*My + k*My +j;
        indx_next    = indx_f + e[k][1]*My + e[k][0];
        f[indx_next] = f_str[indx_f];
      }
    }
  }
}
void init() {
  int i, j, k, indx, indx_f;
  for(i = start; i < end; i++) {
    for(j = 0; j < My; j++) {
      indx = i*My+j;
      rho[indx] = Rho_init;
      u[indx]   = 0.0;
      v[indx]   = 0.0;
      Fk_eq(u[indx],v[indx],rho[indx]);
      for(k = 0; k < Q; k++ ){
        indx_f        = i*Q*My + k*My +j;
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
    for (j = 2; j < My-2; j++) {
    // bottom wall - south side - no slip
      indx_bot = 1*Q*My + j;
      u[indx_bot] = 0.0;
      v[indx_bot] = 0.0;
      rho[indx_bot]     = f[0*My+indx_bot]+f[1*My+indx_bot]+f[3*My+indx_bot]+2.0*(f[4*My+indx_bot]+f[7*My+indx_bot]+f[8*My+indx_bot]);
      f[2*My+indx_bot]  = f[4*My+indx_bot];
      f[5*My+indx_bot]  = f[7*My+indx_bot] - 0.5*(f[1*My+indx_bot] - f[3*My+indx_bot]);
      f[6*My+indx_bot]  = f[8*My+indx_bot] + 0.5*(f[1*My+indx_bot] - f[3*My+indx_bot]);
    }
    for(i = start+1; i <= end; i++) {
    // left wall - inlet - eastside
      indx_left    = Q*My*i + 1;
      u_wall       = rho_in - (f[0*My+indx_left]+f[2*My+indx_left]+f[4*My+indx_left]+2.*(f[3*My+indx_left]+f[7*My+indx_left]+f[6*My+indx_left]));
      u[indx_left] = u_wall;
      v[indx_left] = 0.0;
      rho[indx_left]    = rho_in;
      f[1*My+indx_left] = f[3*My+indx_left] + two_by_three*u_wall;//
      f[5*My+indx_left] = f[7*My+indx_left] - 0.5*(f[2*My+indx_left] - f[4*My+indx_left]) + one_by_six*u_wall;
      f[8*My+indx_left] = f[6*My+indx_left] + 0.5*(f[2*My+indx_left] - f[4*My+indx_left]) + one_by_six*u_wall;//

    // right wall - outlet - westside
      indx_right   = Q*My*i + My - 2;
      u_wall       = -1.0*rho_out + (f[0*My+indx_right]+f[2*My+indx_right]+f[4*My+indx_right]+2.*(f[1*My+indx_right]+f[5*My+indx_right]+f[8*My+indx_right]));
      u[indx_right] = u_wall;
      v[indx_right] = 0.0;
      rho[indx_right]    = rho_out;
      f[3*My+indx_right] = f[1*My+indx_right] - two_by_three*u_wall;
      f[7*My+indx_right] = f[5*My+indx_right] + 0.5*(f[2*My+indx_right] - f[4*My+indx_right]) - one_by_six*u_wall;
      f[6*My+indx_right] = f[8*My+indx_right] - 0.5*(f[2*My+indx_right] - f[4*My+indx_right]) - one_by_six*u_wall;
    }
    // bottom left corner
    indx_bot          = 1 + Q*My;
    f[1*My+indx_bot]  = f[3*My+indx_bot];
    f[2*My+indx_bot]  = f[4*My+indx_bot];
    f[5*My+indx_bot]  = f[7*My+indx_bot];
    f[6*My+indx_bot]  = 0.5*(rho_in - (f[0*My+indx_bot]+f[1*My+indx_bot]+f[2*My+indx_bot]+f[3*My+indx_bot]+f[4*My+indx_bot]+f[5*My+indx_bot]+f[7*My+indx_bot]));
    f[8*My+indx_bot]  = f[6*My+indx_bot];
    // Bottom right corner
    indx_bot          = Q*My - 2 + My;
    f[3*My+indx_bot]  = f[1*My+indx_bot];
    f[2*My+indx_bot]  = f[4*My+indx_bot];
    f[6*My+indx_bot]  = f[8*My+indx_bot];
    f[5*My+indx_bot]  = 0.5*(rho_out - (f[0*My+indx_bot]+f[1*My+indx_bot]+f[2*My+indx_bot]+f[3*My+indx_bot]+f[4*My+indx_bot]+f[6*My+indx_bot]+f[8*My+indx_bot]));
    f[7*My+indx_bot]  = f[5*My+indx_bot];
  }
  else if(taskid == numworkers){
    for (i = 1; i < My - 2; i++) {
    // top wall - north side - no slip
      indx_top         = i + (end-2)*Q*My;

      u[indx_top] = 0.0;
      v[indx_top] = 0.0;
      rho[indx_top]     = f[0*My+indx_top]+f[1*My+indx_top]+f[3*My+indx_top]+2.0*(f[2*My+indx_top]+f[5*My+indx_top]+f[6*My+indx_top]);

      f[4*My+indx_top] = f[2*My+indx_top];
      f[7*My+indx_top] = f[5*My+indx_top] + 0.5*(f[1*My+indx_top] - f[3*My+indx_top]);
      f[8*My+indx_top] = f[6*My+indx_top] - 0.5*(f[1*My+indx_top] - f[3*My+indx_top]);
    }
    for(i = start; i <= (end-1); i++) {
    // left wall - inlet - eastside
      indx_left    = Q*My*i + 1;
      u_wall       = rho_in - (f[0*My+indx_left]+f[2*My+indx_left]+f[4*My+indx_left]+2.*(f[3*My+indx_left]+f[7*My+indx_left]+f[6*My+indx_left]));
      u[indx_left] = u_wall;
      v[indx_left] = 0.0;
      rho[indx_left]    = rho_in;
      f[1*My+indx_left] = f[3*My+indx_left] + two_by_three*u_wall;//
      f[5*My+indx_left] = f[7*My+indx_left] - 0.5*(f[2*My+indx_left] - f[4*My+indx_left]) + one_by_six*u_wall;
      f[8*My+indx_left] = f[6*My+indx_left] + 0.5*(f[2*My+indx_left] - f[4*My+indx_left]) + one_by_six*u_wall;//

    // right wall - outlet - westside
      indx_right   = Q*My*i + My - 2;
      u_wall       = -1.0*rho_out + (f[0*My+indx_right]+f[2*My+indx_right]+f[4*My+indx_right]+2.*(f[1*My+indx_right]+f[5*My+indx_right]+f[8*My+indx_right]));
      u[indx_right] = u_wall;
      v[indx_right] = 0.0;
      rho[indx_right]    = rho_out;
      f[3*My+indx_right] = f[1*My+indx_right] - two_by_three*u_wall;
      f[7*My+indx_right] = f[5*My+indx_right] + 0.5*(f[2*My+indx_right] - f[4*My+indx_right]) - one_by_six*u_wall;
      f[6*My+indx_right] = f[8*My+indx_right] - 0.5*(f[2*My+indx_right] - f[4*My+indx_right]) - one_by_six*u_wall;
    }
    // top left corner
    indx_top         = 1 + Q*(end-2)*My;
    u[indx_top] = 0.0;
    v[indx_top] = 0.0;
    f[1*My+indx_top] = f[3*My+indx_top];
    f[4*My+indx_top] = f[2*My+indx_top];
    f[8*My+indx_top] = f[6*My+indx_top];
    rho[indx_top]     = f[0*My+indx_top]+f[1*My+indx_top]+f[3*My+indx_top]+2.0*(f[2*My+indx_top]+f[5*My+indx_top]+f[6*My+indx_top]);
    f[5*My+indx_top] = 0.5*(rho_in - (f[0*My+indx_top]+f[1*My+indx_top]+f[2*My+indx_top]+f[3*My+indx_top]+f[4*My+indx_top]+f[6*My+indx_top]+f[8*My+indx_top]));
    f[7*My+indx_top] = f[5*My+indx_top];

    //top right corner
    indx_top         = My - 2 + Q*(end-2)*My;
    u[indx_top] = 0.0;
    v[indx_top] = 0.0;
    f[3*My+indx_top] = f[1*My+indx_top];
    f[4*My+indx_top] = f[2*My+indx_top];
    f[7*My+indx_top] = f[5*My+indx_top];
    rho[indx_top]     = f[0*My+indx_top]+f[1*My+indx_top]+f[3*My+indx_top]+2.0*(f[2*My+indx_top]+f[5*My+indx_top]+f[6*My+indx_top]);
    f[6*My+indx_top] = 0.5*(rho[indx_top] - (f[0*My+indx_top]+f[1*My+indx_top]+f[2*My+indx_top]+f[3*My+indx_top]+f[4*My+indx_top]+f[5*My+indx_top]+f[7*My+indx_top]));
    f[8*My+indx_top] = f[6*My+indx_top];
  }else{
    for(i = start; i <= end; i++) {
    // left wall - inlet - eastside
      indx_left    = My*i + 1;
      u_wall       = rho_in - (f[0*My+indx_left]+f[2*My+indx_left]+f[4*My+indx_left]+2.*(f[3*My+indx_left]+f[7*My+indx_left]+f[6*My+indx_left]));
      u[indx_left] = u_wall;
      v[indx_left] = 0.0;
      rho[indx_left]    = rho_in;
      f[1*My+indx_left] = f[3*My+indx_left] + two_by_three*u_wall;//
      f[5*My+indx_left] = f[7*My+indx_left] - 0.5*(f[2*My+indx_left] - f[4*My+indx_left]) + one_by_six*u_wall;
      f[8*My+indx_left] = f[6*My+indx_left] + 0.5*(f[2*My+indx_left] - f[4*My+indx_left]) + one_by_six*u_wall;//

    // right wall - outlet - westside
      indx_right   = Q*My*i + My - 2;
      u_wall       = -1.0*rho_out + (f[0*My+indx_right]+f[2*My+indx_right]+f[4*My+indx_right]+2.*(f[1*My+indx_right]+f[5*My+indx_right]+f[8*My+indx_right]));
      u[indx_right] = u_wall;
      v[indx_right] = 0.0;
      rho[indx_right]    = rho_out;
      f[3*My+indx_right] = f[1*My+indx_right] - two_by_three*u_wall;
      f[7*My+indx_right] = f[5*My+indx_right] + 0.5*(f[2*My+indx_right] - f[4*My+indx_right]) - one_by_six*u_wall;
      f[6*My+indx_right] = f[8*My+indx_right] - 0.5*(f[2*My+indx_right] - f[4*My+indx_right]) - one_by_six*u_wall;
    }
  }
}

void mpiexchange(int taskid, double *c, int M) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Send(&c[end*M],      M, MPI_DOUBLE,  right_node, LTAG,      MPI_COMM_WORLD);
      MPI_Recv(&c[(end+1)*M],  M, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
    }
    source  = left_node;
    msgtype = LTAG;
    MPI_Send(&c[start*M],      M, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    MPI_Recv(&c[0],            M, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG;
       MPI_Recv(&c[0],         M, MPI_DOUBLE,   source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&c[start*M],   M, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Send(&c[(end)*M],    M, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
      MPI_Recv(&c[(end+1)*M],  M, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
    }
  }
}

void sendtomaster(int taskid, double *c, int M) {
  dest = MASTER;

  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&c[0],     rows*M, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&c[M], rows*M, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  }
}
void receivefrmworker(double *c, int M) {
  int rank;
  for (rank=1; rank <= numworkers; rank++) {
    source = rank;
    MPI_Recv(&offset,             1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,               1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,          1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node,         1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&c[offset*M],    rows*M,    MPI_DOUBLE,    source,   WRITE,  MPI_COMM_WORLD, &status);
  }
}

void write2file1( int t, int m, double *c, double *d, char *fname) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles/%s_%d.dat", fname, t);
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
  sprintf(filename,"./datafiles/%s_%d.dat",fname, t);
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
