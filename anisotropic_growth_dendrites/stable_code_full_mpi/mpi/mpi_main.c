#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "constants.h"
#include "variables.h"
#include "functions.h"
#include "mpi_variables.h"

int main(int argc, char *argv[]){
  int t;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  master_allocate_memory(MESHX, MESHX);
  worker_allocate_memory(MESHX, MESHX);
  numworkers = numtasks - 1;

  if ( taskid == MASTER ) {
    phasefield_initialize();
    mpi_distribute(MESHX);
    for ( t = 0; t < phi_timesteps; t++ ) {
      if((t%save_phi) == 0) {
       receivefrmworker();
       write2file_phi(t, MESHX, phi_old);
      }
    }
  }else{
    mpi_distribute(MESHX);

    for ( t = 0; t < phi_timesteps; t++ ) {
      boundary_mpi(taskid, phi_old);
      boundary_mpi(taskid, mu_old);
      mpiexchange( taskid, phi_old, MESHX );
      mpiexchange( taskid, mu_old, MESHX );
      laplacian(phi_old, lap_phi, MESHX);
      laplacian(mu_old,  lap_mu,  MESHX);
      solverloop(start, end);
      update(phi_old, phi_new, MESHX);
      update(mu_old, mu_new, MESHX);
      if((t%save_phi) == 0) {
       sendtomaster();
      }
    }
  }
}
