double comp_res(double *y,double *f,double *a_x,double *a_y,int mesh_tracker,long temp_x,long temp_y,long i,long j)
{
	
	
	long array_index,array_index_up,array_index_below,array_index_left,array_index_right;

	long mesh_inter_index_x,mesh_inter_index_y,mesh_inter_index_up,mesh_inter_index_below,mesh_inter_index_left,mesh_inter_index_right;


	double residual;

	int tag;

	
	//this are all calculations in the fine meshes
	array_index=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*temp_y;

	array_index_up=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*(temp_y+1);
	
	array_index_below=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*(temp_y-1);

	array_index_left=(start[mesh_tracker]+(temp_x-1))+nodes[mesh_tracker].x*temp_y;
	
	array_index_right=(start[mesh_tracker]+(temp_x+1))+nodes[mesh_tracker].x*temp_y;

	mesh_inter_index_x=temp_x+start_a[mesh_tracker].x+mesh_intervals[mesh_tracker].x*temp_y;//it corresponds to the current grid point and one to the right of it 

	mesh_inter_index_y=temp_y+start_a[mesh_tracker].y+mesh_intervals[mesh_tracker].y*temp_x;//it corresponds to the current grid point and one to the above of it 

	mesh_inter_index_left=(temp_x+start_a[mesh_tracker].x-1)+mesh_intervals[mesh_tracker].x*temp_y;//it corresponds to the current grid point and one to the left of it 

	mesh_inter_index_right=(temp_x+start_a[mesh_tracker].x+1)+mesh_intervals[mesh_tracker].x*temp_y;//it corresponds to the point to the right of the current grid point and one to the right of it 
 
	mesh_inter_index_up=(temp_y+start_a[mesh_tracker].y+1)+mesh_intervals[mesh_tracker].y*temp_x;//it corresponds to the point to the above of the current grid point and one to the above of it 

	mesh_inter_index_below=(temp_y+start_a[mesh_tracker].y-1)+mesh_intervals[mesh_tracker].y*temp_x;//it corresponds to the current grid point and one to the below of it 


	//if(temp_x==0 && temp_y==0)
	//{
	//	tag=1;	
	//	residual=f[array_index]-(((a_x[mesh_inter_index_x]*(y[array_index_right]-y[array_index]))/delta_space_sq[mesh_tracker]) + ((a_y[mesh_inter_index_y]*(y[array_index_up]-y[array_index]))/delta_space_sq[mesh_tracker]));
	//}

	//else if(temp_x==0 && temp_y!=0)
	//{
	//	tag=2;
	//	residual=f[array_index]-(((-2.0*a_x[mesh_inter_index_x]*y[array_index]+2.0*a_x[mesh_inter_index_x]*y[array_index_right])/delta_space_sq[mesh_tracker]) + ((a_y[mesh_inter_index_y]*y[array_index_up]-(a_y[mesh_inter_index_y]+a_y[mesh_inter_index_below])*y[array_index]+a_y[mesh_inter_index_below]*y[array_index_below])/delta_space_sq[mesh_tracker]));	
	//}

	//else if(temp_x!=0 && temp_y==0)
	//{
	//	tag=3;
	//	residual=f[array_index]-(((a_x[mesh_inter_index_left]*y[array_index_left]-(a_x[mesh_inter_index_left]+a_x[mesh_inter_index_x])*y[array_index]+a_x[mesh_inter_index_x]*y[array_index_right])/delta_space_sq[mesh_tracker]) + ((2.0*a_y[mesh_inter_index_y]*y[array_index_up]-2.0*a_y[mesh_inter_index_y]*y[array_index])/delta_space_sq[mesh_tracker]));
	//}
		
	//else
	//{
	//	tag=4;			 
		residual=f[array_index]-(((a_x[mesh_inter_index_left]*y[array_index_left]-(a_x[mesh_inter_index_left]+a_x[mesh_inter_index_x])*y[array_index]+a_x[mesh_inter_index_x]*y[array_index_right])/delta_space_sq[mesh_tracker]) + ((a_y[mesh_inter_index_y]*y[array_index_up]-(a_y[mesh_inter_index_y]+a_y[mesh_inter_index_below])*y[array_index]+a_y[mesh_inter_index_below]*y[array_index_below])/delta_space_sq[mesh_tracker]));
	//}


	//printf("res=%lf\n",residual);

	//getchar();

	//if(fabs(residual)>1e-4)
	//{
	//	printf("res=%lf at i=%ld and j=%ld tag=%d\n",residual,i,j,tag);
	//	printf("temp_x=%ld temp_y=%ld\n",temp_x,temp_y);
	//	getchar();
	//}	

	return residual;
}
