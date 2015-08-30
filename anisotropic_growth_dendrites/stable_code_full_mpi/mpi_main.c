#include "mpi.h"
#include "time.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "constants.h"
#include "variables.h"
#include "binary_solver_2D_mpi.c"
#include "mpi_variables.h"
int t;

int main(int argc, char *argv[]){

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  numworkers = numtasks - 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  master_allocate_memory(MESHX, MESHX);
  worker_allocate_memory(MESHX, MESHX);
  if (taskid == MASTER) {
    phasefield_initialize();
    mpi_distribute(MESHX);
    for ( t = 0; t < phi_timesteps; t++ ) {
    }
  }else{
    mpi_distribute(MESHX);
    for ( t = 0; t < phi_timesteps; t++ ) {
    }
  }
}

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
      MPI_Send(&phi_old[offset*pmesh],      rows*My,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&mu_old[offset*pmesh],       rows*My,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
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
