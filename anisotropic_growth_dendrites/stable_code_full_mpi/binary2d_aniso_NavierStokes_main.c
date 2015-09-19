#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "constants.h"
#include "variables.h"
#include "MPI_def.h"
#include "fluid_functions.c"
#include "alloy_functions.c"
#include "mpi_functions.c"
#include "mpi_fn.c"

void mpi_distribute(int Mx);
void write2file( int t, int m, double *c, char *fname);
void write2file1( int t, int m, double *c, double *d, char *fname);

void main(int argc, char *argv[]){
  int t;
  char fname1[100], fname2[100];
  int m = sprintf(fname1,"phi");
  int n = sprintf(fname2,"velocities");

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  numworkers = numtasks - 1;
  allocate_memory_phasefields(MESHX);
#ifdef FLUID
  allocate_memory_FluidVariables(MESHX,pmesh);
#endif
  if(taskid == MASTER) {
    initialize_phasefields();
#ifdef FLUID
    initialize_velocities(MESHX);
#endif
    mpi_distribute(MESHX);
    for (t = 0; t <= phi_timesteps; t++){
#ifdef FLUID
      gauss_siedel();
#endif
      if (t%save_phi == 0) {
        receivefrmworker(phi_old);
        receivefrmworker(u_old);
        receivefrmworker(v_old);
        write2file( t, MESHX, phi_old, fname1);
        write2file1( t, MESHX, u_old, v_old, fname2);
      }
    }
  }else {
    mpi_distribute(MESHX);

    for (t = 0; t <= phi_timesteps; t++){
      mpiexchange(taskid, phi_old, MESHX);
      mpiexchange(taskid, mu_old, MESHX);
      boundary_mpi(taskid, phi_old, MESHX);
      boundary_mpi(taskid, mu_old, MESHX);
      laplacian(phi_old, lap_phi, MESHX);
      laplacian(mu_old, lap_mu, MESHX);
#ifdef FLUID
      computeH(u_old,v_old,Hx,Hy);
      RHS_fn(Hx,Hy,rhs_fn,MESHX, pmesh);
      LHS_fn(MESHX, pmesh);
      boundary_pressure(taskid);
      gauss_siedel();
      ns_solver(start, end, phi_old);
      update_velocities(MESHX);
#endif
      solverloop();
      phi_update();
      if (t%save_phi == 0) {
        sendtomaster(taskid, phi_old);
      }
    }
  }
}
void mpi_distribute(int Mx){
  if ( taskid == MASTER ) {
    averow    =   Mx/numworkers;
    extra     =   Mx%numworkers;
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
      MPI_Send(&phi_old[offset*Mx],      rows*Mx,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&mu_old[offset*Mx],       rows*Mx,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&u_old[offset*Mx],       rows*Mx,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
      MPI_Send(&v_old[offset*Mx],       rows*Mx,          MPI_DOUBLE,      dest,   BEGIN,  MPI_COMM_WORLD);
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
        MPI_Recv(&phi_old[0],   rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&mu_old[0],    rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
#ifdef FLUID
        MPI_Recv(&u_old[0],   rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&v_old[0],    rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
#endif
      }
      else {
        MPI_Recv(&phi_old[Mx],  rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&mu_old[Mx],   rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
#ifdef FLUID
        MPI_Recv(&u_old[Mx],  rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&v_old[Mx],   rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
#endif
      }
      end = rows-1;
    } else {
      MPI_Recv(&phi_old[Mx],    rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&mu_old[Mx],     rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
#ifdef FLUID
      MPI_Recv(&u_old[Mx],    rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&v_old[Mx],     rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
#endif
      end = rows;
    }
  }
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
void write2file1( int t, int m, double *c, double *d, char *fname) {
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
      fprintf(fp,"%d %d %le %le\n",j,i,c[z],d[z]);

    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
