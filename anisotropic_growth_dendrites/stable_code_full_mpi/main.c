#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "constants.h"
#include "variables.h"
#include "gauss_siedel_mpi.c"
#include "fluid_solver_mpi.c"
#include "binary_solver_2D_mpi.c"
#include "time.h"

int t;

int main(int argc, char *argv[]){

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  numworkers = numtasks - 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  if(taskid == MASTER){ // This code fragment runs in the master porcesses
    // clock_t start_t, end_t, total_t;
    allocate_memory();
    phi_initialize();
    fluid_initialize();
    gs_allocate();
    // start_t = clock();
    // printf("Starting of the program, start_t = %ld\n", start_t);
    for ( t = 0; t < phi_timesteps; t++ ) {
    #ifdef growth
      neuman_boundary(phi_old, MESHX);
      neuman_boundary(mu_old, MESHX);
      concentration(phi_old, mu_old, conc, MESHX);
      neuman_boundary(conc, MESHX);
      laplacian(phi_old, lap_phi, MESHX);
      laplacian(mu_old,  lap_mu,  MESHX);
      anisotropic_solverloop();
      update(phi_old, phi_new, MESHX);
      update(mu_old, mu_new, MESHX);
      if((t%save_phi) == 0) {
       write2file_phi(t, MESHX,phi_old);
      }
    #endif
    #ifdef FLUID
      if (t>SMOOTH) {
        fluid_solver();
        if((t%save_fluid) ==0) {
             write2file_fluid (t,u_old,v_old,MESHX);
        }
      }
    #endif
      printf("t=%d\n",t);
    }
    free_memory();

    // end_t = clock();
    // printf("End of the big loop, end_t = %ld\n", end_t);
    // total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    // printf("Total time taken by CPU: %f\n", total_t  );
    printf("Exiting of the program...\n");

  } else { // This code fragment runs in the worker processes
    gs_allocate();
    for ( t = 0; t < phi_timesteps; t++ ) {
      if ( t > SMOOTH ) {
        gs_mpi();
      }
    }
    free(P);
    free(rhs_fn);
    free(a_x);
    free(a_y);
  }
  MPI_Finalize();
  return(0);
}
