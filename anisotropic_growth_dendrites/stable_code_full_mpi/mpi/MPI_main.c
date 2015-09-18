#include "mpi.h"
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "constants.h"
#include "variables.h"
//-----------------------------------------------------------------------------//
#define MASTER  0
#define NONE    0
#define BEGIN   999
#define LTAG    777
#define RTAG    666
#define WRITE   555
#define ERROR   888
#define BREAK   111
//----------------------------------------------------------------------------//
int numtasks, numworkers, taskid, rank, dest;
int averow, extra, offset;
int left_node, right_node;
int start, end;
int source, msgtype;
int rows;
MPI_Status status;

void allocate_memory(int Mx);
void phasefield_initialize();
void mpi_distribute(int Mx);
void mpiexchange(int taskid, double *c, int Mx);
void sendtomaster(int taskid, double *c);
void receivefrmworker(double *c)
int main(int argc, char *argv[]) {
  int t;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  numworkers = numtasks - 1;

  allocate_memory(MESHX);
  if (taskid == MASTER) {
    printf("Hello World, from the master\n");
    phasefield_initialize();
    mpi_distribute(MESHX);
    receivefrmworker(phi_old);
  }
  else {
    printf("Hello Wrold, from worker no. %d\n",taskid);
    mpi_distribute(MESHX);
    mpiexchange(taskid, phi_old, MESHX);
    mpiexchange(taskid, mu_old, MESHX);
    sendtomaster(taskid, phi_old);
  }
  MPI_Finalize();
  return 0;
}
void allocate_memory(int Mx) {
  if (taskid == MASTER) {
    phi_old     =   (double *)malloc(Mx*Mx*sizeof(double));
    mu_old      =   (double *)malloc(Mx*Mx*sizeof(double));
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
      }
      else {
        MPI_Recv(&phi_old[Mx],  rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
        MPI_Recv(&mu_old[Mx],   rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      }
      end = rows-1;
    } else {
      MPI_Recv(&phi_old[Mx],    rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      MPI_Recv(&mu_old[Mx],     rows*Mx,          MPI_DOUBLE,      source,   BEGIN,  MPI_COMM_WORLD, &status);
      end = rows;
    }
  }
}
void mpiexchange(int taskid, double *c, int Mx) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      MPI_Send(&c[end*Mx],      Mx, MPI_DOUBLE,  right_node, LTAG,      MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&c[(end+1)*Mx],  Mx, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
    }
    MPI_Send(&c[start*Mx],      Mx, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    source  = left_node;
    msgtype = LTAG;
    MPI_Recv(&c[0],             Mx, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG;
       MPI_Recv(&c[0],          Mx, MPI_DOUBLE,   source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&c[start*Mx],   Mx, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&c[(end+1)*Mx],  Mx, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&c[(end)*Mx],    Mx, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
    }
  }
}
void sendtomaster(int taskid, double *c) {
  dest = MASTER;

  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&c[0],     rows*MESHX, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&c[MESHX], rows*MESHX, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  }
}
void receivefrmworker(double *c) {
  int rank;
  for (rank=1; rank <= numworkers; rank++) {
    source = rank;
    MPI_Recv(&offset,             1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,               1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,          1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node,         1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&c[offset*MESHX],    rows*MESHX,    MPI_DOUBLE,    source,   WRITE,  MPI_COMM_WORLD, &status);
  }
}
