void dirich_boundary(double *y)
{
	int i,j;

	long array_index;

	//I am going to set my boundaries at the finest grid

	//setting the dirichlet boundary condition
	//currently the dirichlet boundary conditions are at 'right' and 'up'
	
	//setting the one at 'right' 
	
	i=nodes[0].x-1;

	for(j=0;j<=nodes[0].y-1;j++)
	{
		array_index=(i+start[0])+nodes[0].x*j;
		
		y[array_index]=RIGHT;
	}

	//setting the one at 'left' 
	
	i=0;

	for(j=0;j<=nodes[0].y-1;j++)
	{
		array_index=(i+start[0])+nodes[0].x*j;
		
		y[array_index]=LEFT;
	}



	//setting the one at "up"
	j=nodes[0].y-1;							
	
	for(i=1;i<=nodes[0].x-2;i++)
	{
		array_index=(i+start[0])+nodes[0].x*j;
		
		y[array_index]=UP;
	}   

	//setting the one at "below"
	j=0;							
	
	for(i=1;i<=nodes[0].x-2;i++)
	{
		array_index=(i+start[0])+nodes[0].x*j;
		
		y[array_index]=BELOW;
	}   
}		
  
