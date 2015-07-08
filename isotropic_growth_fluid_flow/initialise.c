void init(double *y)
{
	long i;//this is going to iterate over the flattened array indices
	
	
	//in principle I could have initialised anything for the finest grid and kept the other grids at zero

	//initialising the array in the finest grid 

	//we already know that the level zero is the finest geid

	for(i=start[0];i<=end[0];i++) //I have initialised everything including boundary points for the finest grid; the boundary conditions will bw set correctly when I call the next function from main 
	{
		#ifdef ZERO //I am going to fill it with zeroes 
 		y[i]=0.0;
		#endif

		#ifdef RANDOM //I am going to fill it with random nos. which can be useful when I solve y''=0 
		y[i]=(2.0*ran2(&SEED)-1.0)*NOISE_AMPLITUDE;
		#endif
	}

	//initialising the other portions of the solution array; they must be initialised to zero regardless of the exact problem I am trying to solve

	//they are just solving the residual equation so the initial guess as well as the boundary condition should be zero

	for(i=start[1];i<=end[NO_OF_LEVELS-1];i++) //we are starting at the grid which is second finest and we go till the end of the array
	{
		y[i]=0.0;
	} 
	 
	  	

}