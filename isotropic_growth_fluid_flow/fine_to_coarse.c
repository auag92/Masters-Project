void ftoc(double *y,double *f,double *a_x,double *a_y,int mesh_tracker)
{
	//the value of the residuals that has to be found 
	double res_2i_2j,res_2i_min_1_2j,res_2i_plus_1_2j;
	double res_2i_2j_min_1,res_2i_min_1_2j_min_1,res_2i_plus_1_2j_min_1;	
	double res_2i_2j_plus_1,res_2i_min_1_2j_plus_1,res_2i_plus_1_2j_plus_1;	
	
	
		
	long i,j;	

	long array_index_i_j;//the coarse grid array index 

	long temp_x,temp_y;//local coordinates for the finer grid

	int coarser_mesh_tracker=mesh_tracker+1;
	
	//I have to first initialise the coarse grid array to zero (because it might have a solution from a previous iteration
	for(j=0;j<=nodes[coarser_mesh_tracker].y-1;j++)
	{	
		for(i=0;i<=nodes[coarser_mesh_tracker].x-1;i++)//I am iterating over the entire coarse grid array including the boundary locations
		{
			array_index_i_j=(i+start[coarser_mesh_tracker])+nodes[coarser_mesh_tracker].x*j;		
		 
			y[array_index_i_j]=0.0;
		}
	}		
				
	//I am going to employ full weighting for the interpolation
	//the range of iterations is selected by the algorithm 
	//I am iterating over the points in the coarser grid
	//I am simply leaving out the boundary points as in a coarser grid for a Dirichlet problem they are always zero
	//If I don't leave out the boundary points I'll end up accessing array indices which are out of bounds  	

	for(j=1;j<=nodes[coarser_mesh_tracker].y-2;j++)
	{ 
		for(i=1;i<=nodes[coarser_mesh_tracker].x-2;i++)	
		{
			//I have left out those nodes which correspond to the Dirichelt situation

			//===========================================================================//
			//computing the residuals in the finer grid (see "compute_residual.c")
			//we have to compute 9 in total at every location in the grid

			//at the locations where no-flux is implemented
			//if(i==0 && j==0)	
			//{

				
				//computing res_2i_2j

			//	temp_x=2*i;				
			//	temp_y=2*j;
		
			//	res_2i_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);

				//computing res_2i_plus_1_2j
			
			//	temp_x=2*i+1;				
			//	temp_y=2*j;
		
			//	res_2i_plus_1_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_2j_plus_1
			
			//	temp_x=2*i;				
			//	temp_y=2*j+1;
		
			//	res_2i_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_plus_1_2j_plus_1
			
			//	temp_x=2*i+1;				
			//	temp_y=2*j+1;
		
			//	res_2i_plus_1_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);


				
				//========================================================================//
				//moving the errors to the coarser grid using the full weighing operator  		
				//array_index_i_j=(i+start[coarser_mesh_tracker])+nodes[coarser_mesh_tracker].x*j;

				//f[array_index_i_j]=0.0625*(res_2i_plus_1_2j_plus_1);

				//f[array_index_i_j]+=0.125*(res_2i_2j_plus_1+res_2i_plus_1_2j);

				//f[array_index_i_j]+=0.25*res_2i_2j;	

				//printf("f=%lf\n",f[array_index_i_j]);	
				//========================================================================//

			//}


			//else if (i==0 && j!=0)
			//{
				//computing res_2i_2j

				//temp_x=2*i;				
				//temp_y=2*j;
		
				//res_2i_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_plus_1_2j
			
				//temp_x=2*i+1;				
				//temp_y=2*j;
		
				//res_2i_plus_1_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_2j_min_1
				//temp_x=2*i;				
				//temp_y=2*j-1;
		
				//res_2i_2j_min_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	
				

				//computing res_2i_plus_1_2j_min_1
			
				//temp_x=2*i+1;				
				//temp_y=2*j-1;
		
				//res_2i_plus_1_2j_min_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);

				//computing res_2i_2j_plus_1
			
				//temp_x=2*i;				
				//temp_y=2*j+1;
		
				//res_2i_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_plus_1_2j_plus_1
			
				//temp_x=2*i+1;				
				//temp_y=2*j+1;
		
				//res_2i_plus_1_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);

				
				//========================================================================//
				//moving the errors to the coarser grid using the full weighing operator  		
				//array_index_i_j=(i+start[coarser_mesh_tracker])+nodes[coarser_mesh_tracker].x*j;

				//f[array_index_i_j]=0.0625*(res_2i_plus_1_2j_min_1+res_2i_plus_1_2j_plus_1);

				//f[array_index_i_j]+=0.125*(res_2i_2j_min_1+res_2i_2j_plus_1+res_2i_plus_1_2j);

		//		f[array_index_i_j]+=0.25*res_2i_2j;	

				//printf("f=%lf\n",f[array_index_i_j]);	
				//========================================================================//
	

//			}

//			else if (i!=0 && j==0)
//			{
				//computing res_2i_2j

//				temp_x=2*i;				
//				temp_y=2*j;
		
//				res_2i_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_min_1_2j
//				temp_x=2*i-1;				
//				temp_y=2*j;
		
//				res_2i_min_1_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_plus_1_2j
			
//				temp_x=2*i+1;				
//				temp_y=2*j;
		
//				res_2i_plus_1_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);

				//computing res_2i_2j_plus_1
			
//				temp_x=2*i;				
//				temp_y=2*j+1;
		
//				res_2i_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);


				//computing res_2i_min_1_2j_plus_1
			
//				temp_x=2*i-1;				
//				temp_y=2*j+1;
		
//				res_2i_min_1_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_plus_1_2j_plus_1
			
//				temp_x=2*i+1;				
//				temp_y=2*j+1;
		
//				res_2i_plus_1_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

	
				
				//=======================================================================//
				//moving the errors to the coarser grid using the full weighing operator  		
//				array_index_i_j=(i+start[coarser_mesh_tracker])+nodes[coarser_mesh_tracker].x*j;

//				f[array_index_i_j]=0.0625*(res_2i_min_1_2j_plus_1+res_2i_plus_1_2j_plus_1);

//				f[array_index_i_j]+=0.125*(res_2i_2j_plus_1+res_2i_min_1_2j+res_2i_plus_1_2j);

//				f[array_index_i_j]+=0.25*res_2i_2j;	

				//printf("f=%lf\n",f[array_index_i_j]);	
				//========================================================================//
					

//			}
	
			//else
			//{	
				//computing res_2i_2j

				temp_x=2*i;				
				temp_y=2*j;
		
				res_2i_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	


				//computing res_2i_min_1_2j
						
				temp_x=2*i-1;				
				temp_y=2*j;
		
				res_2i_min_1_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	
			

				//computing res_2i_plus_1_2j
			
				temp_x=2*i+1;				
				temp_y=2*j;
		
				res_2i_plus_1_2j=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_2j_min_1
				temp_x=2*i;				
				temp_y=2*j-1;
		
				res_2i_2j_min_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	
			


				//computing res_2i_min_1_2j_min_1
				temp_x=2*i-1;				
				temp_y=2*j-1;
		
				res_2i_min_1_2j_min_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);
				

				//computing res_2i_plus_1_2j_min_1
			
				temp_x=2*i+1;				
				temp_y=2*j-1;
		
				res_2i_plus_1_2j_min_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);

				//computing res_2i_2j_plus_1
			
				temp_x=2*i;				
				temp_y=2*j+1;
		
				res_2i_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	


				//computing res_2i_min_1_2j_plus_1
			
				temp_x=2*i-1;				
				temp_y=2*j+1;
		
				res_2i_min_1_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);	

				//computing res_2i_plus_1_2j_plus_1
			
				temp_x=2*i+1;				
				temp_y=2*j+1;
		
				res_2i_plus_1_2j_plus_1=comp_res(y,f,a_x,a_y,mesh_tracker,temp_x,temp_y,i,j);

				//========================================================================//
				//moving the errors to the coarser grid using the full weighing operator  		
				array_index_i_j=(i+start[coarser_mesh_tracker])+nodes[coarser_mesh_tracker].x*j;

				f[array_index_i_j]=0.0625*(res_2i_min_1_2j_min_1+res_2i_min_1_2j_plus_1+res_2i_plus_1_2j_min_1+res_2i_plus_1_2j_plus_1);

				f[array_index_i_j]+=0.125*(res_2i_2j_min_1+res_2i_2j_plus_1+res_2i_min_1_2j+res_2i_plus_1_2j);

				f[array_index_i_j]+=0.25*res_2i_2j;	

				//printf("f=%lf\n",f[array_index_i_j]);	
				//=======================================================================//

			//}	
		}	

	}
}				
		
