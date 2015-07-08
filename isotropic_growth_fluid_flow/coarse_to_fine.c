void ctof(double *y,int mesh_tracker)
{
	long i,j;
	long array_index_2i_2j,array_index_2i_plus_1_2j,array_index_2i_2j_plus_1,array_index_2i_plus_1_2j_plus_1;//these are related to the finer grid

	long array_index_i_j,array_index_i_plus_1_j,array_index_i_j_plus_1,array_index_i_plus_1_j_plus_1;//these are related to thre coarser grid

	int finer_mesh_tracker=mesh_tracker-1;


	//interpolating the solutions


	//================================================================================================//
	//the first loop will map the entries over the entire coarse mesh to array-indices which are twice that of them in the fine mesh in both x and y, i.e., it is going to tackle the "2i,2j"case in the fine grid

	//it requires no interpolation (so including boundary points in both dimensions)

	for(j=0;j<nodes[mesh_tracker].y;j++)//I am including boundary points (in a domain having even no. of sub-intervals the boundaries are always at even locations) but they shouldn't get updated because the boundaries of coarser meshes are set to zero (this is only true for a Dirichlet boundary condition) 
	//and I need to update from the coarser meshes for the neumann condition 			  
	{

		for(i=0;i<nodes[mesh_tracker].x;i++) //I am including boundary points (in a domain having even no. of sub-intervals the boundaries are always at even locations) but they shouldn't get updated because the boundaries of coarser meshes are set to zero (this is only true for a Dirichlet boundary condition)  			  
		//and I need to update from the coarser meshes for the neumann condition 
		{
			//for the coarser grid
			array_index_i_j=(i+start[mesh_tracker])+nodes[mesh_tracker].x*j;

			//for the finer grid
			array_index_2i_2j=((2*i)+start[finer_mesh_tracker])+nodes[finer_mesh_tracker].x*(2*j);
		
			//interpolating		
			y[array_index_2i_2j]+=y[array_index_i_j];
		} 

	}
	//============================================================================================//
		

	//============================================================================================//
	//this second loop will just tackle the "2i+1,2j" case in the finer grid 

	//requires interpolation only along the x-direction, so including the right hand boundary in the y-direction

	for(j=0;j<nodes[mesh_tracker].y;j++)	
	{

		for(i=0;i<nodes[mesh_tracker].x-1;i++)//this is going over the coarser grid 
		{
			//computing the different array indices
	
			//for the coarser grid
			array_index_i_j=(i+start[mesh_tracker])+nodes[mesh_tracker].x*j;
			array_index_i_plus_1_j=((i+1)+start[mesh_tracker])+nodes[mesh_tracker].x*j;
				

			//for the finer grid
			array_index_2i_plus_1_2j=((2*i+1)+start[finer_mesh_tracker])+nodes[finer_mesh_tracker].x*(2*j);		

			//interpolating to the finer grid
			y[array_index_2i_plus_1_2j]+=0.5*(y[array_index_i_j]+y[array_index_i_plus_1_j]);

		}
	}

	//=============================================================================================//

	//============================================================================================//
	//this third loop will just tackle the "2i,2j+1" case in the finer grid 

	//requires interpolation only along the y-direction, so including the right boundary in the x-direction

	for(j=0;j<nodes[mesh_tracker].y-1;j++)	
	{

		for(i=0;i<nodes[mesh_tracker].x;i++)//this is going over the coarser grid 
		{
			//computing the different array indices
	
			//for the coarser grid
			array_index_i_j=(i+start[mesh_tracker])+nodes[mesh_tracker].x*j;
			array_index_i_j_plus_1=(i+start[mesh_tracker])+nodes[mesh_tracker].x*(j+1);
				

			//for the finer grid
			array_index_2i_2j_plus_1=((2*i)+start[finer_mesh_tracker])+nodes[finer_mesh_tracker].x*(2*j+1);		

			//interpolating to the finer grid
			y[array_index_2i_2j_plus_1]+=0.5*(y[array_index_i_j]+y[array_index_i_j_plus_1]);

		}
	}

	//=============================================================================================//

	//============================================================================================//
	//this third loop will just tackle the "2i+1,2j+1" case in the finer grid 

	//requires interpolation along both y-direction, so leaving the right boundary

	for(j=0;j<nodes[mesh_tracker].y-1;j++)	
	{

		for(i=0;i<nodes[mesh_tracker].x-1;i++)//this is going over the coarser grid 
		{
			//computing the different array indices
	
			//for the coarser grid
			array_index_i_j=(i+start[mesh_tracker])+nodes[mesh_tracker].x*j;
			array_index_i_j_plus_1=(i+start[mesh_tracker])+nodes[mesh_tracker].x*(j+1);
			array_index_i_plus_1_j=((i+1)+start[mesh_tracker])+nodes[mesh_tracker].x*j;
			array_index_i_plus_1_j_plus_1=((i+1)+start[mesh_tracker])+nodes[mesh_tracker].x*(j+1);	

			//for the finer grid
			array_index_2i_plus_1_2j_plus_1=((2*i+1)+start[finer_mesh_tracker])+nodes[finer_mesh_tracker].x*(2*j+1);		

			//interpolating to the finer grid
			y[array_index_2i_plus_1_2j_plus_1]+=0.25*(y[array_index_i_j]+y[array_index_i_j_plus_1]+y[array_index_i_plus_1_j]+y[array_index_i_plus_1_j_plus_1]);

		}
	}

	//=============================================================================================//

} 
			
