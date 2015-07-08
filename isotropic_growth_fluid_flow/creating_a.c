void cr_a(double *a_x,double *a_y)
{
	int i,j,m;

	long mesh_arr_index,finer_mesh_array_index_x,finer_mesh_array_index_x_right,finer_mesh_array_index_y,finer_mesh_array_index_y_up;

	



	//'a' is known a-priori

	
	for(m=0;m<NO_OF_LEVELS;m++)
	{	

		//====================================================================================//
		//setting a_x
		for(j=0;j<nodes[m].y;j++)
		{	
			for(i=0;i<mesh_intervals[m].x;i++)
			//I am going to set it everywhere regardless of the BCs   
			{
				mesh_arr_index=i+start_a[m].x+mesh_intervals[m].x*j;

				if(m==0)
				{
					#ifdef A_ONE//mimics the simpler problem we have started with
	 				a_x[mesh_arr_index]=1.0;
					#endif

					#ifdef A_SIX_TIMES_X_SQ
					a_x[mesh_arr_index]=6.0*(SYS_LEFT_END+((i+0.5)*delta_space[m]))*(SYS_LEFT_END+((i+0.5)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	

					#ifdef A_SIX_TIMES_Y_SQ
					a_x[mesh_arr_index]=6.0*(SYS_BOTTOM_END+((j)*delta_space[m]))*(SYS_BOTTOM_END+((j)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	


					#ifdef A_1000_TIMES_X_SQ
					a_x[mesh_arr_index]=1000.0*(SYS_LEFT_END+((i+0.5)*delta_space[m]))*(SYS_LEFT_END+((i+0.5)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	

					#ifdef A_1000_TIMES_Y_SQ
					a_x[mesh_arr_index]=1000.0*(SYS_BOTTOM_END+((j)*delta_space[m]))*(SYS_BOTTOM_END+((j)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	

				}

				else
				{
					finer_mesh_array_index_x=2*i+start_a[m-1].x+mesh_intervals[m-1].x*2*j;

					finer_mesh_array_index_x_right=2*i+1+start_a[m-1].x+mesh_intervals[m-1].x*2*j;

					a_x[mesh_arr_index]=0.5*(a_x[finer_mesh_array_index_x]+a_x[finer_mesh_array_index_x_right]);

				}


			}	
					
		}
		//====================================================================================//
	

		//====================================================================================//
		//setting a_y
		for(i=0;i<nodes[m].x;i++)
		{	
			for(j=0;j<mesh_intervals[m].y;j++)
			//I am going to set it everywhere regardless of the BCs   
			{
				mesh_arr_index=j+start_a[m].y+mesh_intervals[m].y*i;
			
				if(m==0)
				{
					#ifdef A_ONE//mimics the simpler problem we have started with
	 				a_y[mesh_arr_index]=1.0;
					#endif

					#ifdef A_SIX_TIMES_X_SQ
					a_y[mesh_arr_index]=6.0*(SYS_LEFT_END+((i)*delta_space[m]))*(SYS_LEFT_END+((i)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	

					#ifdef A_SIX_TIMES_Y_SQ
					a_y[mesh_arr_index]=6.0*(SYS_BOTTOM_END+((j+0.5)*delta_space[m]))*(SYS_BOTTOM_END+((j+0.5)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	


					#ifdef A_1000_TIMES_X_SQ
					a_y[mesh_arr_index]=1000.0*(SYS_LEFT_END+((i)*delta_space[m]))*(SYS_LEFT_END+((i)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	

					#ifdef A_1000_TIMES_Y_SQ
					a_y[mesh_arr_index]=1000.0*(SYS_BOTTOM_END+((j+0.5)*delta_space[m]))*(SYS_BOTTOM_END+((j+0.5)*delta_space[m]));//as the '0' node coincides with the left end of the system
					#endif	

				}


				else
				{
					finer_mesh_array_index_y=2*j+start_a[m-1].y+mesh_intervals[m-1].y*2*i;

					finer_mesh_array_index_y_up=2*j+1+start_a[m-1].y+mesh_intervals[m-1].y*2*i;

					a_y[mesh_arr_index]=0.5*(a_y[finer_mesh_array_index_y]+a_y[finer_mesh_array_index_y_up]);
				}

			}	
		}
		//====================================================================================//
	}
}		
