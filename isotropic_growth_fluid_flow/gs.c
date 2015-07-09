void Gauss_siedel(double *P, double *fn, double *a_x, double *a_y) {
  long i, j, indx, indx_rght, indx_frnt, indx_lft, indx_bck;
  double P_old;
  double error;
  double tol=1.0e-6;
  // for(i=1; i < pmesh-2; i++) {
  //   for(j=1; j< pmesh-2; j++) {
  //     indx          = i*MESHX    + j;
  //     indx_rght     = indx      + 1;
  //     indx_frnt     = indx      + MESHX;
  //     indx_lft      = indx      - 1;
  //     indx_bck      = indx      - MESHX;
  //
  //     A[indx]      = 1.0 - 0.5*(phi[indx] + phi[indx_frnt];  //
  //     B[indx]      = (zeta(phi[indx_bck])  + zeta(phi[indx]));  //back
  //     E[indx]      = (zeta(phi[indx_frnt]) + zeta(phi[indx])); //front
  //     F[indx]      = -(A[indx] + B[indx] + C[indx] + E[indx]);
  //   }
  // }
  for(;;) {
    for(i=1; i < pmesh-1; i++) {
      for(j=1;j < pmesh-1;j++) {
        indx         = i*pmesh    + j;
        indx_rght    = indx      + 1;
        indx_frnt    = indx      + MESHY;
        indx_lft     = indx      - 1;
        indx_bck     = indx      - MESHY;

        if (((i+j)%2) == 0) {
          P_old     = P[indx];

          P[indx]  = P[indx_lft] + P[indx_rght]
          + P[indx_bck] + P[indx_frnt] + deltax2*fn[indx];
          P[indx] /= (-F[indx]);
        }
      }
    }
    for(i=1; i < MESHX-1; i++) {
      for(j=1;j < MESHY-1; j++) {
        indx         = i*MESHY    + j;
        indx_rght    = indx      + 1;
        indx_frnt    = indx      + MESHY;
        indx_lft     = indx      - 1;
        indx_bck     = indx      - MESHY;

        if (((i+j)%2) != 0) {
          V_old     = V[indx];

          P[indx]  = a_x[]P[indx_lft] + P[indx_rght]
          + P[indx_bck] + P[indx_frnt] - deltax2*fn[indx];
          V[indx] /= (-F[indx]);
        }
      }
    }
    error = compute_error(A, B, C, E, F, V);
    if (fabs(error) < tol) {
      break;
    }
  }
}
double compute_error(double *A, double *B, double *C, double *E, double *F, double *V) {
  double error=0.0;
  long indx, indx_rght, indx_frnt, indx_lft, indx_bck, i, j;
  for(i=1; i < MESHX-1; i++) {
    for(j=1;j < MESHY-1; j++) {
      indx  = i*MESHY    + j;
      indx_right         = indx      + 1;
      indx_front         = indx      + MESHY;
      indx_left          = indx      - 1;
      indx_back          = indx      - MESHY;

      error += fabs(A[indx]*V[indx_lft] + B[indx]*V[indx_rght]
	     + C[indx]*V[indx_bck] + E[indx]*V[indx_frnt] + F[indx]*V[indx]);
    }
  }
  return(error);
}
