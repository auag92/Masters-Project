#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "constants.h"
#include "variables.h"

//--------------------------------------------------------------
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

long i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck, indx_ax, indx_ay;
double error;
double tol=1.0e-6;
int iter  = 0;

void main(int argc, char *argv[]){
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  numworkers = numtasks - 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
}
