void Gauss_siedel(double *P, double *fn, double *a_x, double *a_y) {
  long i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck;
  double P_old;
  double error;
  double tol=1.0e-6;

  for(;;) {
    for(i=1; i < pmesh-1; i++) {
      for(j=1;j < pmesh-1;j++) {
        indx         = i*pmesh      + j;
        indx_rght    = indx         + 1;
        indx_frnt    = indx         + MESHY;
        indx_lft     = indx         - 1;
        indx_bck     = indx         - MESHY;
        indx_ax      = i*(pmesh-1)  +j;
        indx_ay      = i*(pmesh-2)   +j;

        if (((i+j)%2) == 0) {
          P_old     = P[indx];

          P[indx]  = a_x[indx_ax-1]*P[indx_lft] + a_x[indx_ax]*P[indx_rght]
          + a_y[indx_ay-(pmesh-2)]*P[indx_bck] + a_y[indx_ay]*P[indx_frnt] + deltax2*fn[indx];
          P[indx] /= a_x[indx_ax-1]+a_x[indx_ax]+a_y[indx_ay-(pmesh-2)]+a_y[indx_ay];
        }
      }
    }
    for(i=1; i < pmesh-1; i++) {
      for(j=1;j < pmesh-1;j++) {
        indx         = i*pmesh      + j;
        indx_rght    = indx         + 1;
        indx_frnt    = indx         + MESHY;
        indx_lft     = indx         - 1;
        indx_bck     = indx         - MESHY;
        indx_ax      = i*(pmesh-1)  +j;
        indx_ay      = i*(pmesh-2)   +j;

        if (((i+j)%2) != 0) {
          P_old     = P[indx];

          P[indx]  = a_x[indx_ax-1]*P[indx_lft] + a_x[indx_ax]*P[indx_rght]
          + a_y[indx_ay-(pmesh-2)]*P[indx_bck] + a_y[indx_ay]*P[indx_frnt] + deltax2*fn[indx];
          P[indx] /= a_x[indx_ax-1]+a_x[indx_ax]+a_y[indx_ay-(pmesh-2)]+a_y[indx_ay];
        }
      }
    }
    error = compute_error(a_x,a_y,P);
    if (fabs(error) < tol) {
      break;
    }
  }
}
double compute_error(double *A, double *B, double *C, double *E, double *F, double *V) {
  double error=0.0;
  long indx, indx_rght, indx_frnt, indx_lft, indx_bck, i, j;
  for(i=1; i < pmesh-1; i++) {
    for(j=1;j < pmesh-1; j++) {
      indx  = i*pmesh    + j;
      indx_right         = indx      + 1;
      indx_front         = indx      + MESHY;
      indx_left          = indx      - 1;
      indx_back          = indx      - MESHY;
      indx_ax      = i*(pmesh-1)  +j;
      indx_ay      = i*(pmesh-2)   +j;

      error += fabs(a_x[indx_ax-1]*P[indx_lft] + a_x[indx_ax]*P[indx_rght]
      + a_y[indx_ay-(pmesh-2)]*P[indx_bck] + a_y[indx_ay]*P[indx_frnt] +
      a_x[indx_ax-1]+a_x[indx_ax]+a_y[indx_ay-(pmesh-2)]+a_y[indx_ay] - deltax2*fn[indx]);
    }
  }
  return(error);
}
for (m=1;m<4;m++){
  d1 = c - t*(t+1);
  d2 = c;
  t = t/2;
  for (i=0;i<t+1;i++){
    for (j=0;j<t;j++){

      x = d2 + i*t+j;
      x1 = d1 + 2*t*2*i + 2*j;
      x2 = d1 + 2*i*2*t + 2*j +1;
      a_x[x] = 0.5*(a_x[x1] + a_x[x2]);

      y = d2 +  + j*(t+1)+i;
      y1 = d1 + 2*j*(2*t+1) + 2*i;
      y2 = d1 + (2*j+1)*(2*t+1) + 2*i;
      a_y[y] = 0.5*(a_y[y1] + a_y[y2]);
    }
  }
}
}
