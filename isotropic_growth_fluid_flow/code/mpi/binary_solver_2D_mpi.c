#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "constants.h"
#include "variables.h"
//------------------------------------------------------------
#define MASTER 0
#define NONE 0
#define BEGIN 999
#define LTAG  777
#define RTAG  666
#define WRITE 555
//-------------------------------------------------------------
int numtasks, numworkers, taskid, rank, dest;
int averow, extra, offset;
int left_node, right_node;
int start, end;
int t;

void phi_initialize();
void solverloop(int Mx, int My);

void main(int argc, char *argv[]){

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  numworkers = numtasks - 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  if(taskid == MASTER){
    phi_old = (double *)malloc(MESHX*MESHX*sizeof(double));
    mu_old = (double *)malloc(MESHX*MESHX*sizeof(double));
    phi_initialize();

    averow  = MESHX/numworkers;
    extra   = MESHX%numworkers;
    offset = 0;

    for (rank=1; rank <= (numworkers); rank++) {
     rows = (rank <= extra) ? averow+1 : averow;
     left_node  = rank - 1;
     right_node = rank + 1;
     if(rank == 1) {
       left_node  = NONE;
     }
     if(rank == (numworkers)) {
       right_node = NONE;
     }
     dest = rank;
     MPI_Send(&offset,          1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
     MPI_Send(&rows,            1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
     MPI_Send(&left_node,       1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
     MPI_Send(&right_node,      1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
     MPI_Send(&phi_old[offset], rows*MESHX, MPI_DOUBLE,      dest, BEGIN, MPI_COMM_WORLD);
     MPI_Send(&mu_old[offset],  rows*MESHX, MPI_DOUBLE,      dest, BEGIN, MPI_COMM_WORLD);
     offset = offset + rows;
   }

   for (t=0; t < phi_timesteps; t++) {
     if (t%saveT == 0) {
       receivefrmworker();
       writetofile(t);
     }
   }
   free(phi_old, mu_old);
  }
  if(taskid != MASTER) {
    source  = MASTER;
    msgtype = BEGIN;
    MPI_Recv(&offset,     1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,       1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);

    start = 1;
    if ((taskid == 1) || (taskid == numworkers)) {
      phi_old = (double *)malloc((rows+2)*MESHX*sizeof(double));
      mu_old = (double *)malloc((rows+2)*MESHX*sizeof(double));

      if (taskid == 1) {
        MPI_Recv(&phi_old[1], rows*MESHX, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&mu_old[1], rows*MESHX, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
      }
      else {
        MPI_Recv(&phi_old[1], rows*MESHX, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&mu_old[1], rows*MESHX, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
      }
      end   = rows-1;
    }
    else {
      phi_old   =   (double *)malloc((rows+2)*MESHX*sizeof(double));
      mu_old    =   (double *)malloc((rows+2)*MESHX*sizeof(double));
      MPI_Recv(&phi_old[1], rows, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&mu_old[1], rows, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
      end = rows;
    }
    for (t=1; t < ntimesteps; t++) {
      mpiexchange(taskid);
      lap_phi   =   (double *)malloc((rows+2)*MESHX*sizeof(double));
      lap_mu    =   (double *)malloc((rows+2)*MESHX*sizeof(double));
      laplacian(phi_old, lap_phi, rows, MESHX);
      laplacian(mu_old,  lap_mu, rows, MESHX);

      solverloop(comp);
      apply_boundary_conditions(taskid);
      if (t%saveT == 0) {
        sendtomaster(taskid);
      }
    }
    free(phi_old, mu_old);
  }
  MPI_Finalize();
}
void solverloop(int Mx, int My){
  int i,j,z;
  double p,dp_dt,dmu_dt, kai;
  for (i=1; i < (My-1); i++) {
    for (j=1; j < (Mx-1); j++){
      z= i*Mx + j;

      p = phi_old[z];
      dp_dt = (G*E*lap_phi[z] - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p) + (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p))/(tau*E);
      phi_new[z] = p + deltat*(dp_dt);

      dmu_dt = Mob*lap_mu[z] - V_gradC - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
      kai = 1+(K-1)*p*p*(3-2*p);

      mu_new[z] = mu_old[z]  + deltat*dmu_dt/kai;
    }
  }
}

void phi_initialize() {
  long i,j,z;
  double r;
#ifdef Centre
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      r= (i-pmesh*0.5)*(i-pmesh*0.5) + (j-pmesh*0.5)*(j-pmesh*0.5);
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
}
void laplacian(double *f, double *lap, int Mx, int My) {
  long i,j,z;

  for (i=1; i< My -1; i++)
  {
    for (j=1; j< Mx -1; j++)
    {
      z= i*M + j;
      lap[z] = (f[z-1] + f[z+1] + f[z+M] + f[z-M] -4.0*f[z])*inv_deltax2;
    }
  }
}
