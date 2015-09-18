#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define deltat .01
#define tsteps 10000

double fn1(double x, double t);
void rk4(double *x, int t, int h);
void write2file(double *x, int t, int h);
main() {
  double *x;
  double *v;
  x = (double *)(tsteps*sizeof(double));
  v = (double *)(tsteps*sizeof(double));
  x[0] = 0.0; // initial value
  v[0] = 1.0; // initial value
  rk4(x, tsteps, deltat);
  write2file(x, tsteps, deltat);
  free(x);
}
void rk4(double *x, int t, int h) {
  int i;
  double a, b, c, d;
  for (i = 0; i < (tsteps-1); i++) {
    a = fn1( x[i], i*h);
    b = fn1( x[i]+h*0.5*a, i*h+0.5*h);
    c = fn1( x[i]+0.5*h*b, i*h+0.5*h);
    d = fn1( x[i]+h*c, (i+1)*h)
    x[i+1] = x[i] + h/6*(a + 2*b + 2*c + d);
  }
}
double fn1(double x, double t) {
  double ans;
  /*Calculate the derivative*/
  return ans;
}
void write2file(double *x, int t, int h){
  FILE *fp;
  fp = fopen("evolution.dat","w");
  int i;
  for (i = 0; i<t; i++){
    fprintf(fp, "%le %le\n",h*i,x[i]);
  }
  fclose(fp);
}
