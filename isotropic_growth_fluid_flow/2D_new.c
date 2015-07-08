#include"stdio.h"
#include"math.h"
#include"stdlib.h"

#define MESHX 100
#define MESHY 50

#define deltat (0.005)
#define deltax (1.0)
#define deltay (1.0)
#define deltax2 (deltax*deltax)
#define deltay2 (deltay*deltay)
#define D (1.0)
#define M (1.0)
#define omega (1.0)

#define gamma (1.0)        /*Surface energy*/
#define epsilon (3.0)     /*Interface-width*/
#define tau (1.0)          /*relaxation constant*/
#define mu_eq (0.0)

#define ntimesteps (200001)
#define saveT (5000)

#define CONSTANT_V_left  (0.1)
#define CONSTANT_V_right (0.1)
#define CONSTANT_V_back  (0.0)
#define CONSTANT_V_front (10.0)

#define POTENTIAL_GRADIENT ((CONSTANT_V_front - CONSTANT_V_back)/MESHX)

 #define NEUMANN
//#define PERIODIC
#define DIMENSION 2

#define X 0
#define Y 1

#define K 1.2
#define mu_0 1.0

#define e 1.0
#define Z 1.0

#define zeta_alpha (1.0)
#define zeta_beta (2.0)

#define CIRCLE

#ifdef CIRCLE
#define centerX (MESHX/2)
#define centerY (MESHY/2)
#define radius  (15.0)
#endif

#define mu_initial (-gamma/radius)

double phi[MESHX*MESHY], deltaphi[MESHX*MESHY] ,deltamu[MESHX*MESHY], V[MESHX*MESHY];
double mu[MESHX*MESHY],  lap_phi[MESHX*MESHY], grad_flux[MESHX*MESHY][DIMENSION], Ddcdmu[MESHX*MESHY][DIMENSION];

double inv_deltax2 = (1.0/deltax2);
double inv_deltay2 = (1.0/deltay2);

void initialize(double *phi, double *mu, double *V);
void boundary(double *phi, double *mu, double *V);
void write2file (double *phi, double *mu, double *V, long t);
void update(double *phi, double *mu);
void laplacian(double *f, double *lap);
double dc_dmu(double phi, double mu);
double dc_dphi(double phi, double mu);
double zeta(double phi);
void calculate_flux(double *phi, double *mu, double *V, double (*grad_flux)[DIMENSION], double (*Ddcdmu)[DIMENSION]);
void calculate_Ddcdmu(double *phi, double *mu, double (*Ddcdmu)[DIMENSION]);
void Gauss_siedel(double *phi, double *mu,double *V);

void copyYZ(long copyto, long copyfrom, double *phi, double *mu);
void copyXZ(long copyto, long copyfrom, double *phi, double *mu);
void apply_DIRICHLET_X_0(double *V);
void apply_DIRICHLET_X_END(double *V);
void apply_DIRICHLET_Y_0(double *V);
void apply_DIRICHLET_Y_END(double *V);


double compute_error(double *A, double *B, double *C, double *E, double *F, double *V);

void main() {
  long i, j, index, index_left, index_back, t;
  double div_flux;

  initialize(phi, mu, V);
  boundary(  phi, mu, V);
  write2file(phi, mu, V, 0);


  for (t=0; t < ntimesteps; t++) {
    if (t%5==0) {
      printf("Iteration=%ld\n", t);
    }
    Gauss_siedel(phi, mu, V);

    //Finding the update in phi
    laplacian(phi,     lap_phi);
    for (i=0; i < (MESHX); i++) {
      for (j=0; j < (MESHY); j++) {
	index           = i*MESHY    + j;
        deltaphi[index] = (deltat/(tau*epsilon))
	                *(2.0*(gamma*epsilon)*lap_phi[index]
	                - 18.0*(gamma/epsilon)*(phi[index]*(1.0-phi[index])*(1.0-2.0*phi[index]))
			+ (mu[index]-mu_eq)*6.0*(phi[index]*(1.0-phi[index])));
      }
    }
    //Ending update of phi

    //Finding update in mu
    calculate_Ddcdmu(phi, mu, Ddcdmu);
    calculate_flux(phi, mu, V, grad_flux, Ddcdmu);


    for (i=1; i < (MESHX-1); i++) {
      for (j=1; j < (MESHY-1); j++) {
	index       = i*MESHY    + j;
	index_left     = index - 1;
	index_back     = index - MESHY;

	div_flux       = (grad_flux[index][Y] - grad_flux[index_left][Y])*inv_deltay2;
	div_flux      += (grad_flux[index][X] - grad_flux[index_back][X])*inv_deltax2;

	deltamu[index] = (deltat*div_flux  -
	               dc_dphi(phi[index], mu[index])*(deltaphi[index]))/dc_dmu(phi[index],mu[index]);
      }
    }
    //update in mu
    update(phi, mu);
    if(t%saveT==0) {
     write2file(phi, mu, V, t);
    }
  }
}

double dc_dmu(double phi, double mu) {
  return(1.0 + ( K - 1.0)*(phi*phi)*(3.0 - 2.0*phi));
}

double dc_dphi(double phi, double mu) {
  return(((K*mu + mu_0) - mu)*6.0*phi*(1.0 - phi));
}

void calculate_Ddcdmu(double *phi, double *mu, double (*Ddcdmu)[DIMENSION]) {
  long i, j,index, index_right, index_front;
  for (i=0; i < (MESHX-1); i++) {
    for (j=0; j < (MESHY-1); j++) {
      index               = i*MESHY    + j;
      index_right         = index      + 1;
      index_front         = index      + MESHY;

      Ddcdmu[index][X]    = D*(0.5/omega)*(phi[index_front]*(1.0-phi[index_front])*dc_dmu(phi[index_front], mu[index_front])
			  + phi[index]*(1.0-phi[index])*dc_dmu(phi[index_front], mu[index_front]));
      Ddcdmu[index][Y]    = D*(0.5/omega)*(phi[index_right]*(1.0-phi[index_right])*dc_dmu(phi[index_right], mu[index_right])
			  + phi[index]*(1.0-phi[index])*dc_dmu(phi[index_right], mu[index_right]));
    }
  }
}

void calculate_flux(double *phi, double *mu, double *V, double (*grad_flux)[DIMENSION], double (*Ddcdmu)[DIMENSION]) {
  long i, j,index, index_right, index_front;
  for (i=1; i < (MESHX-1); i++) {
    for (j=1; j < (MESHY-1); j++) {
      index               = i*MESHY    + j;
      index_right         = index      + 1;
      index_front         = index      + MESHY;

      grad_flux[index][X] = Ddcdmu[index][X]*(omega*(mu[index_front] - mu[index]) - e*Z*(V[index_front] - V[index]));
      grad_flux[index][Y] = Ddcdmu[index][Y]*(omega*(mu[index_right] - mu[index]) - e*Z*(V[index_right] - V[index]));
    }
  }
}

double zeta(double phi) {
  return(zeta_alpha*phi + zeta_beta*(1.0-phi));
}
void laplacian(double *f, double *lap) {
  long i, j, index, index_right, index_front, index_left, index_back;
   for(i=1; i < MESHX-1; i++) {
    for(j=1; j< MESHY-1; j++) {
      index               = i*MESHY    + j;
      index_right         = index      + 1;
      index_front         = index      + MESHY;
      index_left          = index      - 1;
      index_back          = index      - MESHY;

      lap[index]          = (f[index_right] + f[index_front] + f[index_back] + f[index_left] -4.0*f[index])*inv_deltax2;
    }
  }
}
void Gauss_siedel(double *phi, double *mu, double *V) {
   long i, j, index, index_right, index_front, index_left, index_back;
   double V_old,A[MESHX*MESHY],B[MESHX*MESHY],C[MESHX*MESHY],E[MESHX*MESHY],F[MESHX*MESHY];
   double error;
   double tol=1.0e-6;

   for(i=1; i < MESHX-1; i++) {
      for(j=1; j< MESHY-1; j++) {
	indx            = i*MESHY    + j;
	indx_rght       = indx      + 1;
	indx_frnt       = indx      + MESHY;
	indx_lft        = indx      - 1;
	indx_bck        = indx      - MESHY;

	A[index]        = (zeta(phi[index_left])  + zeta(phi[index]));  //left
	B[index]        = (zeta(phi[index_right]) + zeta(phi[index]));  //right
	C[index]             = (zeta(phi[index_back])  + zeta(phi[index]));  //back
        E[index]             = (zeta(phi[index_front]) + zeta(phi[index])); //front
        F[index]             = -(A[index] + B[index] + C[index] + E[index]);
      }
   }
   for(;;) {
    for(i=1; i < MESHX-1; i++) {
      for(j=1;j < MESHY-1;j++) {
	index               = i*MESHY    + j;
	index_right         = index      + 1;
	index_front         = index      + MESHY;
	index_left          = index      - 1;
	index_back          = index      - MESHY;

	if (((i+j)%2) == 0) {
	  V_old     = V[index];

	  V[index]  = A[index]*V[index_left] + B[index]*V[index_right]
		    + C[index]*V[index_back] + E[index]*V[index_front];
	  V[index] /= (-F[index]);
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
	  V_old     = V[index];

	  V[index]  = A[index]*V[index_left] + B[index]*V[index_right]
		    + C[index]*V[index_back] + E[index]*V[index_front];
	  V[index] /= (-F[index]);
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
  long index, index_right, index_front, index_left, index_back, i, j;
  for(i=1; i < MESHX-1; i++) {
    for(j=1;j < MESHY-1; j++) {
      index  = i*MESHY    + j;
      index_right         = index      + 1;
      index_front         = index      + MESHY;
      index_left          = index      - 1;
      index_back          = index      - MESHY;

      error += fabs(A[index]*V[index_left] + B[index]*V[index_right]
	     + C[index]*V[index_back] + E[index]*V[index_front] + F[index]*V[index]);
    }
  }
  return(error);
}

void initialize(double *phi, double *mu, double *V) {
#ifdef CIRCLE
  long i, j, index;
  for(i=0; i < MESHX; i++) {
    for(j=0; j < MESHY; j++) {
      index = i*MESHY + j;
      if ((i-centerX)*(i-centerX) + (j-centerY)*(j-centerY) <= radius*radius) {
	phi[index] = 0;
      } else {
	phi[index] = 1;
      }
      V[index]  = 0;
      mu[index] = mu_initial;
    }
  }
#endif
}
void copyYZ(long copyto, long copyfrom, double *phi, double *mu) {
  long j, index_to, index_from;
  for (j=0; j < MESHY; j++) {
    index_to      = copyto*MESHY   + j;
    index_from    = copyfrom*MESHY + j;
    phi[index_to] = phi[index_from];
    mu[index_to]  = mu[index_from];
  }
}
void copyXZ(long copyto, long copyfrom, double *phi, double *mu) {
  long i, index_to, index_from;
  for (i=0; i < MESHX; i++) {
    index_to      = i*MESHY   + copyto;
    index_from    = i*MESHY   + copyfrom;
    phi[index_to] = phi[index_from];
    mu[index_to]  = mu[index_from];
  }
}
void apply_DIRICHLET_X_0(double *V) {
  long j, index;
  for (j=0; j < MESHY; j++) {
    index    = j;
    V[index] = CONSTANT_V_back;
  }
}
void apply_DIRICHLET_X_END(double *V) {
  long j, index;
  for (j=0; j < MESHY; j++) {
    index    = (MESHX-1)*MESHY + j;
    V[index] = CONSTANT_V_front;
  }
}
void apply_DIRICHLET_Y_0(double *V) {
  long i, index;
  for (i=0; i < MESHX; i++) {
    index    = i*MESHY;
    V[index] = CONSTANT_V_back + i*(POTENTIAL_GRADIENT);
  }
}
void apply_DIRICHLET_Y_END(double *V) {
  long i, index;
  for (i=0; i < MESHX; i++) {
    index    = i*MESHY + MESHY-1;
    V[index] = CONSTANT_V_back + i*(POTENTIAL_GRADIENT);
  }
}
void boundary(double *phi, double *mu, double *V) {
   copyYZ(0,             1,  phi, mu);
   copyYZ(MESHX-1, MESHX-2,  phi, mu);

   copyXZ(0,             1,  phi, mu);
   copyXZ(MESHY-1, MESHY-2,  phi, mu);

   apply_DIRICHLET_X_0(V);
   apply_DIRICHLET_X_END(V);
   apply_DIRICHLET_Y_0(V);
   apply_DIRICHLET_Y_END(V);
}
void update(double *phi, double *mu) {
  long i, j, index;
  for(i=0; i < MESHX; i++) {
    for(j=0; j < MESHY; j++) {
      index = i*MESHY + j;
      phi[index]+= deltaphi[index];
      mu[index] += deltamu[index];
    }
  }
}
void write2file (double *phi, double *mu, double *V, long t) {
   int i,j,index;
  FILE *fp;
  char filename[1000];

  sprintf(filename,"2D_new_%ld.dat",t);
  fp = fopen(filename,"w");
  if (j=MESHY/2) {
    for(i=0;i<MESHX;i++) {
      index = i*MESHY + j;
      fprintf(fp,"%d %d %le %le %le\n",i,j,phi[index],mu[index],V[index]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
