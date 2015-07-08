void even_rb_solve(double *y,double *f,double *a_x,double *a_y,double *error,int mesh_tracker)
{
	long i,j;
	double y_star,new_val;
	
	long array_index,array_index_up,array_index_below,array_index_left,array_index_right;

	long mesh_inter_index_x,mesh_inter_index_y,mesh_inter_index_up,mesh_inter_index_below,mesh_inter_index_left,mesh_inter_index_right;


	long even_check;

	for(j=1;j<nodes[mesh_tracker].y-1;j++) //leaving out the boundary points for the Dirichlet boundary condition 
	{
		for(i=1;i<nodes[mesh_tracker].x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 	
		{
			even_check=i+j;

			if(even_check%2==0)//the sum of the indices are even 	
			{
				//computing the different array indices
				//this has to be done considering that 'i' includes the starting indices 
				array_index=start[mesh_tracker]+i+nodes[mesh_tracker].x*j;
				array_index_left=(start[mesh_tracker]+i-1)+nodes[mesh_tracker].x*j;
				array_index_right=(start[mesh_tracker]+i+1)+nodes[mesh_tracker].x*j;
				array_index_up=start[mesh_tracker]+i+nodes[mesh_tracker].x*(j+1);
				array_index_below=start[mesh_tracker]+i+nodes[mesh_tracker].x*(j-1);


				mesh_inter_index_x=i+start_a[mesh_tracker].x+mesh_intervals[mesh_tracker].x*j;//it corresponds to the current grid point and one to the right of it 

				mesh_inter_index_y=j+start_a[mesh_tracker].y+mesh_intervals[mesh_tracker].y*i;//it corresponds to the current grid point and one to the above of it 

				mesh_inter_index_left=(i+start_a[mesh_tracker].x-1)+mesh_intervals[mesh_tracker].x*j;//it corresponds to the current grid point and one to the left of it 

				mesh_inter_index_right=(i+start_a[mesh_tracker].x+1)+mesh_intervals[mesh_tracker].x*j;//it corresponds to the point to the right of the current grid point and one to the right of it 
 
				mesh_inter_index_up=(j+start_a[mesh_tracker].y+1)+mesh_intervals[mesh_tracker].y*i;//it corresponds to the point to the above of the current grid point and one to the above of it 

				mesh_inter_index_below=(j+start_a[mesh_tracker].y-1)+mesh_intervals[mesh_tracker].y*i;//it corresponds to the current grid point and one to the below of it 


			
				//if(i==0 && j==0)
				//	y_star=	(a_x[mesh_inter_index_x]*y[array_index_right]+a_y[mesh_inter_index_y]*y[array_index_up]-f[array_index]*delta_space_sq[mesh_tracker])/(a_x[mesh_inter_index_x]+a_y[mesh_inter_index_y]);

				//else if(i==0 && j!=0)				
				//	y_star=(2.0*a_x[mesh_inter_index_x]*y[array_index_right]+a_y[mesh_inter_index_y]*y[array_index_up]+a_y[mesh_inter_index_below]*y[array_index_below]-delta_space_sq[mesh_tracker]*f[array_index])/(2.0*a_x[mesh_inter_index_x]+a_y[mesh_inter_index_y]+a_y[mesh_inter_index_below]);

				//else if(i!=0 && j==0)				
				//	y_star=(y[array_index_up]*2.0*a_y[mesh_inter_index_y]+y[array_index_left]*a_x[mesh_inter_index_left]+y[array_index_right]*a_x[mesh_inter_index_x]-delta_space_sq[mesh_tracker]*f[array_index])/(2.0*a_y[mesh_inter_index_y]+a_x[mesh_inter_index_left]+a_x[mesh_inter_index_x]);
			
				//else
					y_star=(a_x[mesh_inter_index_left]*y[array_index_left]+a_x[mesh_inter_index_x]*y[array_index_right]+a_y[mesh_inter_index_y]*y[array_index_up]+a_y[mesh_inter_index_below]*y[array_index_below]-delta_space_sq[mesh_tracker]*f[array_index])/(a_x[mesh_inter_index_left]+a_x[mesh_inter_index_x]+a_y[mesh_inter_index_y]+a_y[mesh_inter_index_below]);

				
				new_val=y[array_index]+OMEGA*(y_star-y[array_index]);	

				*error+=(new_val-y[array_index])*(new_val-y[array_index]);

				y[array_index]=new_val;
			}
		}
	}

	
}				
