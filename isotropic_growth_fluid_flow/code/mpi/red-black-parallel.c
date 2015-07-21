#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"mpi.h"

#define MESHX 50
#define MESHY 50

#define deltat (0.005)
#define deltax (1.0)
#define deltay (1.0)
#define deltax2 (deltax*deltax)
#define deltay2 (deltay*deltay)

#define CONSTANT_V_left  (0.1)
#define CONSTANT_V_right (0.1)
#define CONSTANT_V_back  (0.0)
#define CONSTANT_V_front (20.0)

#define POTENTIAL_GRADIENT ((CONSTANT_V_front - CONSTANT_V_back)/MESHX)

#define K 1.2
#define mu_0 1.0

#define e 1.0
#define Z 1.0

#define zeta_alpha (5.0)
#define zeta_beta (1.0)
#define CIRCLE
#ifdef CIRCLE
#define centerX (MESHX/2)
#define centerY (MESHY/2)
#define radius  (5.0)
#endif

#define mu_initial (-gamma/radius)

#define MASTER 0
#define NONE 0

#define BEGIN 999
#define LTAG  777
#define RTAG  666
#define WRITE 555
double phi[MESHX*MESHY], deltaphi[MESHX*MESHY] ,deltamu[MESHX*MESHY];
double inv_deltax2 = (1.0/deltax2);
double inv_deltay2 = (1.0/deltay2);
double *comp;
long start, end;

int numtasks, taskid;
int numworkers;
  
int source, msgtype;

int averow, extra, rank, offset, dest;
int left_node, right_node;

long i, t,rows;

FILE *fp;
void initialize(double *phi,double *comp);
void apply_boundary_conditions(double *phi,double *comp);
void writetofile ();
void mpiexchange(int taskid);
void receivefrmworker();
void sendtomaster(int taskid);
double zeta(double phi);
void copyYZ(long copyto, long copyfrom, double *phi);
void copyXZ(long copyto, long copyfrom, double *phi);
void Gauss_siedel(double *phi, double *mu,double *V);
void apply_DIRICHLET_X_0(double *comp);
void apply_DIRICHLET_X_END(double *comp);
void apply_DIRICHLET_Y_0(double *comp);
void apply_DIRICHLET_Y_END(double *comp);
double compute_error(double *A, double *B, double *C, double *E, double *F, double *comp);
MPI_Status status;

void main(int argc, char *argv[]) {
  
  long i, j, index, index_left, index_back;
  int numtasks,taskid,numworker,start,end;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  numworkers = numtasks-1;
if (taskid == MASTER) {
    comp = (double *)malloc((MESHX)*(MESHY)*sizeof(double));

    initialize(phi,comp);
    apply_boundary_conditions(phi,comp);

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

      MPI_Send(&offset,       1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&rows,         1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&left_node,    1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&right_node,   1,    MPI_INT,         dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&comp[offset], rows*MESHY, MPI_DOUBLE,      dest, BEGIN, MPI_COMM_WORLD);

      offset = offset + (rows)*MESHY;
    }

      receivefrmworker();
      writetofile();
    free(comp);
  }

if(taskid != MASTER) {
   long i, j, index, index_right, index_front, index_left, index_back;
   double A[MESHX*MESHY],B[MESHX*MESHY],C[MESHX*MESHY],E[MESHX*MESHY],F[MESHX*MESHY];
   double error;
   double tol=1.0e-6;
    source  = MASTER;
    msgtype = BEGIN;
    MPI_Recv(&offset,     1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,       1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);


    start = 1;
    if ((taskid == 1) || (taskid == numworkers)) {
      comp = (double *)malloc((rows+1)*MESHY*sizeof(double));

      if (taskid == 1) {
        MPI_Recv(&comp[0], rows*MESHY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
      } else {
        MPI_Recv(&comp[1], rows*MESHY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
      }
      end   = rows-1;
    } else {
      comp = (double *)malloc((rows+2)*MESHY*sizeof(double));
      MPI_Recv(&comp[1], rows*MESHY, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
      end = rows;
    }

     for(i=1; i < MESHX-1; i++) {
      for(j=1; j< MESHY-1; j++) {
	index               = i*MESHY    + j;
	index_right         = index      + 1;
	index_front         = index      + MESHY;
	index_left          = index      - 1;
	index_back          = index      - MESHY;
	
	A[index]             = (zeta(phi[index_left])  + zeta(phi[index]));  //left
	B[index]             = (zeta(phi[index_right]) + zeta(phi[index]));  //right
	C[index]             = (zeta(phi[index_back])  + zeta(phi[index]));  //back
        E[index]             = (zeta(phi[index_front]) + zeta(phi[index])); //front
        F[index]             = -(A[index] + B[index] + C[index] + E[index]);
      }
   }

    for (;;) {

      for(i=1; i < MESHX-1; i++) {
      for(j=1;j < MESHY-1;j++) {
	index               = i*MESHY    + j;
	index_right         = index      + 1;
	index_front         = index      + MESHY;
	index_left          = index      - 1;
	index_back          = index      - MESHY;

	if (((i+j)%2) == 0) {
      mpiexchange(taskid);

	  comp[index]  = A[index]*comp[index_left] + B[index]*comp[index_right]
		    + C[index]*comp[index_back] + E[index]*comp[index_front];
	  comp[index] /= (-F[index]);
	}
      }
    }
    for(i=1; i < MESHX-1; i++) {
      for(j=1;j < MESHY-1; j++) {
	index               = i*MESHY    + j;
	index_right         = index      + 1;
	index_front         = index      + MESHY;
	index_left          = index      - 1;
	index_back          = index      - MESHY;

	if (((i+j)%2) != 0) {
      mpiexchange(taskid);


	  comp[index]  = A[index]*comp[index_left] + B[index]*comp[index_right]
		    + C[index]*comp[index_back] + E[index]*comp[index_front];
	  comp[index] /= (-F[index]);
	}
      }
    }
    error = compute_error(A, B, C, E, F, comp);
    if (fabs(error) < tol) {
      break;
    }
  }



    sendtomaster(taskid);
     
    free(comp);
  }
}

void mpiexchange(int taskid) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      MPI_Send(&comp[end],    MESHY, MPI_DOUBLE,  right_node, LTAG,     MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG; 
      MPI_Recv(&comp[end+1],  MESHY, MPI_DOUBLE,  source,     msgtype,  MPI_COMM_WORLD, &status);
    }
    MPI_Send(&comp[start],    MESHY, MPI_DOUBLE,  left_node,  RTAG,     MPI_COMM_WORLD);
    source  = left_node;
    msgtype = LTAG; 
    MPI_Recv(&comp[0],        MESHY, MPI_DOUBLE,  source,     msgtype,  MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG; 
       MPI_Recv(&comp[0],     MESHY, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&comp[start], MESHY, MPI_DOUBLE, left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG; 
      MPI_Recv(&comp[end+1],  MESHY, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&comp[end],    MESHY, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
    }
  }
}

void sendtomaster(int taskid) {
  dest = MASTER;
  
  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&comp[0], rows*MESHY, MPI_DOUBLE,         dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&comp[1], rows*MESHY, MPI_DOUBLE,         dest, WRITE, MPI_COMM_WORLD);
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
    MPI_Recv(&comp[offset], rows*MESHY,    MPI_DOUBLE,      source, WRITE, MPI_COMM_WORLD, &status);
  }
}

void writetofile() {
  long i,j,index;
 
  char name[1000];
  sprintf(name, "voltage.dat");
  fp = fopen(name, "w");
  for (i=0; i < (MESHX); i++) {
    for (j=0;j<MESHY;j++) {
      index = i*MESHY +j;
      
    fprintf(fp, "%le %le %le\n", i*deltax,j*deltay,comp[index]);
  }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
void apply_boundary_conditions(double *phi,double *comp) {
   copyYZ(0,             1,  phi);
   copyYZ(MESHX-1, MESHX-2,  phi);
   
   copyXZ(0,             1,  phi);
   copyXZ(MESHY-1, MESHY-2,  phi);
  
  apply_DIRICHLET_X_0(comp);
   apply_DIRICHLET_X_END(comp);
   apply_DIRICHLET_Y_0(comp);
   apply_DIRICHLET_Y_END(comp);
}

void initialize(double *phi, double *comp) {
  long i,j,index;
  for (i=0; i < MESHX; i++) {
   for(j=0;j<MESHY;j++) {
     index = i*MESHY + j;
#ifdef CIRCLE
      if ((i-centerX)*(i-centerX) + (j-centerY)*(j-centerY) <= radius*radius) {
	phi[index] = 0;
      } else {
	phi[index] = 1;
      }
#endif
     comp[index] = 0.0;
    }
  }
}

void apply_DIRICHLET_X_0(double *comp) {
  long j, index;
  for (j=0; j < MESHY; j++) {
    index    = j;
    comp[index] = CONSTANT_V_back;
  }
}
void apply_DIRICHLET_X_END(double *comp) {
  long j, index;
  for (j=0; j < MESHY; j++) {
    index    = (MESHX-1)*MESHY + j;
    comp[index] = CONSTANT_V_front;
  }
}
void apply_DIRICHLET_Y_0(double *comp) {
  long i, index;
  for (i=0; i < MESHX; i++) {
    index    = i*MESHY;
    comp[index] = CONSTANT_V_back + i*(POTENTIAL_GRADIENT);
  }
}
void apply_DIRICHLET_Y_END(double *comp) {
  long i, index;
  for (i=0; i < MESHX; i++) {
    index    = i*MESHY + MESHY-1;
    comp[index] = CONSTANT_V_back + i*(POTENTIAL_GRADIENT);
  }
}

double zeta(double phi) {
  return(zeta_alpha*phi + zeta_beta*(1.0-phi));
}

void copyYZ(long copyto, long copyfrom, double *phi) {
  long j, index_to, index_from;
  for (j=0; j < MESHY; j++) {
    index_to      = copyto*MESHY   + j;
    index_from    = copyfrom*MESHY + j;
    phi[index_to] = phi[index_from];
    
  }
}
void copyXZ(long copyto, long copyfrom, double *phi) {
  long i, index_to, index_from;
  for (i=0; i < MESHX; i++) {
    index_to      = i*MESHY   + copyto;
    index_from    = i*MESHY   + copyfrom;
    phi[index_to] = phi[index_from];
    
  }
}


double compute_error(double *A, double *B, double *C, double *E, double *F, double *comp) {
  double error=0.0;
  long index, index_right, index_front, index_left, index_back, i, j;
  for(i=1; i < MESHX-1; i++) {
    for(j=1;j < MESHY-1; j++) {
      index  = i*MESHY    + j;
      index_right         = index      + 1;
      index_front         = index      + MESHY;
      index_left          = index      - 1;
      index_back          = index      - MESHY;
      
      error += fabs(A[index]*comp[index_left] + B[index]*comp[index_right] 
	     + C[index]*comp[index_back] + E[index]*comp[index_front] + F[index]*comp[index]);
    }
  }
  return(error);
}


















