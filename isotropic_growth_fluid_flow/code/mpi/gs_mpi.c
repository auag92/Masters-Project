#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"


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
