
/*Phasefield module*/
double phi_new[MESHX2], phi_old[MESHX2];
double mu_new[MESHX2], mu_old[MESHX2];
double lap_phi[MESHX2], lap_mu[MESHX2];
double conc[MESHX2];
//----------------------------------------------------------
/*Fluid module*/
double P[2*pmesh2]; //Pressure
double u_old[MESHX2], u_now[MESHX2]; // velocity in x direction
double v_old[MESHX2], v_now[MESHX2]; // velocity in y direction
double v_str[MESHX2], u_str[MESHX2];
double a_x[2*MESHX2], a_y[2*MESHX2];
double rhs_fn[2*pmesh2]; // rhs function that is passed on to the multigrid module
double Hx[MESHX2],Hy[MESHX2]; // saving the Hx and Hy stuff from the fluid based calculations