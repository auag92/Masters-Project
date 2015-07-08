void rhs(double *f,double *g)
{	
	long i,j;//they will iterate through the flattened array indices

	long array_index;

	//we are going to create the RHS for the finest grid only 

	#ifdef POINT_SOURCE

	long x_index= (long)(LOCATE_POINT_SOURCE_X*nodes[0].x);
		
	long y_index= (long)(LOCATE_POINT_SOURCE_Y*nodes[0].y);	

	printf("The point source is located at=(%ld,%ld)",x_index,y_index);

	#endif
			

	for(j=0;j<=nodes[0].y-1;j++)
	{
		for(i=0;i<=nodes[0].x-1;i++)
		{
			//I have initialised everywhere but I wouldn't solve for the Dirichlet points 
	
			array_index=(start[0]+i)+nodes[0].x*j;//this is a general form 
		
			#ifdef ZERO_RHS
			f[array_index]=0.0;
			#endif

			#ifdef POINT_SOURCE
			if(i==x_index && j==y_index)
			{	
				f[array_index]=POINT_SOURCE_MAGNITUDE;
			}
			else
			{
				f[array_index]=0.0;
			}
			#endif

			#ifdef FUNCTION_RHS
			f[array_index]=1.0;
			#endif
			
			g[array_index]=f[array_index];

		}
	}
}						
