//------------------------------------------------------------------------------
void sendtomaster(int taskid) {
  dest = MASTER;

  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&P[0],     rows*pmesh, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&P[pmesh], rows*pmesh, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  }
}
void receivefrmworker() {
  int rank;
  for (rank=1; rank <= numworkers; rank++) {
    source = rank;
    MPI_Recv(&offset,             1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,               1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,          1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node,         1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&P[offset*pmesh],    rows*pmesh,    MPI_DOUBLE,    source,   WRITE,  MPI_COMM_WORLD, &status);
  }
}
void mpiexchange(int taskid) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      MPI_Send(&P[end*pmesh],       pmesh, MPI_DOUBLE,  right_node, LTAG,      MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*pmesh],  pmesh, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
    }
    MPI_Send(&P[start*pmesh],      pmesh, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    source  = left_node;
    msgtype = LTAG;
    MPI_Recv(&P[0],                pmesh, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG;
       MPI_Recv(&P[0],             pmesh, MPI_DOUBLE,   source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&P[start*pmesh],   pmesh, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*pmesh],  pmesh, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&P[(end)*pmesh],    pmesh, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
    }
  }
}
