
#define pmesh 129
#define pmesh2 (pmesh*pmesh)
#define MESHX (pmesh + 1)
#define MESHX2 (MESHX*MESHX)
#define deltax (1.0)
#define inv_deltax (1.0/deltax)
#define inv_deltax2 (inv_deltax*inv_deltax)
#define deltat (0.01)

//----------------------------
#define K (0.5)             /*Partition Coefficient*/
#define G (1.0)             /*Surface Energy*/
#define Mob (1.0)             /*Mobility*/
#define E (4.0)            /*epsilon - dimensions of length [m]*/
#define tau (1.0)
#define deltaMu (0.4) //
#define Mu (1.0)     //  chemical potential
//-----------------------------------------

#define radius2 100
//---------------------------------------------
//--------fluid States-------------------------
//---------------------------------------------
#define phi_tol 0.5
#define dirichlet_pressure
#define LDC
//#define PipeFlow


#ifdef PipeFlow
#define inv_Re (1)
#define Ue (0.0) // east wall velocity
#define Uw (0.0) // west wall velocity
#define Un (0.0) // north wall velocity
#define Us (0.0) // south wall velocity
#define Ve (0.0) // east wall velocity
#define Vw (0.0) // west wall velocity
#define Vn (0.0) // north wall velocity
#define Vs (0.0) // south wall velocity
//---------------------------------------------
#define p_up      (0)
#define p_down    (0)
#define p_left    (5)
#define p_right   (0)
#endif

#ifdef LDC
#define inv_Re (1)
#define Ue (0.0) // east wall velocity
#define Uw (0.0) // west wall velocity
#define Un (3.0) // north wall velocity
#define Us (0.0) // south wall velocity
#define Ve (0.0) // east wall velocity
#define Vw (0.0) // west wall velocity
#define Vn (0.0) // north wall velocity
#define Vs (0.0) // south wall velocity
//----------------------------------------
#define p_up      (0)
#define p_down    (0)
#define p_left    (0)
#define p_right   (0)
#endif
//------------------------------------------
#define gs_tol 10e-6
