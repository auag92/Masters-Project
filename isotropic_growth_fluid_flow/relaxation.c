double relax(double *y,double *f,double *a_x,double *a_y,int mesh_tracker)
{
	double error=0.0;

	

	//printf("At level %d\n",mesh_tracker);
	
	
	//============================================================================//
	//solving the equation for the even update only(see "even_red_black_solver.c")
		
	even_rb_solve(y,f,a_x,a_y,&error,mesh_tracker);
	//===========================================================================//	

		
	//============================================================================//
	//solving the equation for the odd update only(see "odd_red_black_solver.c")
		
	odd_rb_solve(y,f,a_x,a_y,&error,mesh_tracker);
	//===========================================================================//	
	
	//===========================================================================//
	error=sqrt(error*d_sp[mesh_tracker].x*d_sp[mesh_tracker].y);
	//===========================================================================//
	
	return error;
	
}
