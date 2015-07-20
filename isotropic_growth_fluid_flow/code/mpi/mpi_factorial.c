#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MASTER 0
#define BEGIN 111
#define WRITE 222
int my_rank, numtasks, numworkers, rank;
int source, dest, tag;
MPI_Status status;
long n;
long fact, ans;
long start, end;

long mult(long start, long end){
  if (start != end)
    return start*mult(start-1, end);
  else if(end == 0)
    return 1;
  else
    return end;
}
void get_arg(long *n){
  printf("enter the argument:\n");
  scanf("%ld", n);
//  printf("\n");
}
int main(int argc, char *argv[]){

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  numworkers  = numtasks - 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if(my_rank == 0){
    get_arg(&n);
    int q, rem;
    q   =   n/numworkers;
    rem =   n%numworkers;
    end    =    0;
    start  =    0;
    for (rank=1; rank <= (numworkers); rank++) {
     dest   =    rank;
     end    =    start+1;
     start  =    start + (long)(q);
     if(  rem >= rank ){
       start++;
     }
     MPI_Send(&start,    1,    MPI_LONG,    dest,   BEGIN,  MPI_COMM_WORLD);
     MPI_Send(&end,      1,    MPI_LONG,    dest,   BEGIN,  MPI_COMM_WORLD);
   }
   fact = 1;
   for (rank=1; rank <= numworkers; rank++) {
     source = rank;
     MPI_Recv(&ans,      1,    MPI_LONG,    source, WRITE,  MPI_COMM_WORLD, &status);
     fact = fact*ans;
   }
   printf("The factorial is:%ld\n",fact);
  }
  else{
    source = MASTER;
    MPI_Recv(&start,      1,   MPI_LONG,    source,  BEGIN,  MPI_COMM_WORLD, &status);
    MPI_Recv(&end,        1,   MPI_LONG,    source,  BEGIN,  MPI_COMM_WORLD, &status);
    ans = mult(start, end);
    MPI_Send(&ans,        1,   MPI_LONG,    MASTER,   WRITE,  MPI_COMM_WORLD);
  }
  /*Shutdown MPI*/
  MPI_Finalize();
  return 0;
}
