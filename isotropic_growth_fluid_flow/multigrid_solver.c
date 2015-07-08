void mult_solve(double *y,double *f,double *a_x,double *a_y)
{
	int mesh_tracker=0;//we will first start our computations from the finest grid
		
	int iter=0;//this will count the no. of iterations at a particular mesh level 	

	int tot_iter=0;//this will count the total no. of iterations at all the mesh levels combined

	double error=0.0,old_error;
	
	int count=0;//this variable keeps track of the no. of times we change grids 

	printf("At the level=%d\n",mesh_tracker);

	//relaxing on the grid using red-black gauss-seidel (see "relaxation.c")
	for( ; ; ) //I can break out of this loop only under a particular situation 
	{
		old_error=error;//transferring the old error to the current error		
				
		//solving the equations on the current grid
		//see "relaxation.c"
		error=relax(y,f,a_x,a_y,mesh_tracker);//this is a single iteration of the red-black gauss-seidel
		
		iter++;//updating the no. of iterations at the current level

		tot_iter++;//updating the total no. of iterations

		printf("tot_iter=%d\titer=%d\terror=%lf\n",tot_iter,iter,error);

		printf("ftoc cond=%lf\told_err=%lf\n",old_error*TOLER_RATE,old_error);	

		//there is a minimum no. of iterations to be performed on every grid
		if(iter<MIN_NUM_ITER)
		{
			continue;
		}
							
		//this is the exit condition
		//if after performing enough iterations on the finest the error is below the tolerance we exit the loop
		else if(error<TOLER && mesh_tracker==0)
		{
			//recording the iteration details (see "write_iterations.c")
			write_iter(count,iter,mesh_tracker);		 
		
	
			break;
		}

		//condition for moving to fine grid
		//convergence fine, enough iterations done (see the subtle difference from the previous condition)
		else if(error<TOLER)
		{

			
			//recording the iteration details (see "write_iterations.c")
			write_iter(count,iter,mesh_tracker);	

			//keeping track of the grid exchanges
			count++;	  

			//see "coarse_to_fine.c"
			//this routine will serve an initial guess for the finer grid relaxation
			ctof(y,mesh_tracker);
			
			//updating a new value of mesh tracker
			mesh_tracker--;//a finer grid is the exisiting value -1 

			//re-setting the no. of iterations
			iter=0;

			//re-setting the error
			error=0.0;			
	
			printf("At the level=%d\n",mesh_tracker);
		}

		//condition for moving onto a coarser grid
		//convergence slow, have coarser grids left

		else if(old_error*TOLER_RATE<=error && mesh_tracker<NO_OF_LEVELS-1)			
		{
			//recording the iteration details (see "write_iterations.c")
			write_iter(count,iter,mesh_tracker);	

			//keeping track of the grid exchanges
			count++;	  


			//see "fine_to_coarse.c"						
			ftoc(y,f,a_x,a_y,mesh_tracker);	

			//updating a new value of mesh tracker
			mesh_tracker++;//a coarser grid is the exisiting value +1 

			//re-setting the no. of iterations
			iter=0;

			//re-setting the error
			error=0.0;

			printf("At the level=%d\n",mesh_tracker);
		}

		//getchar();
	}//closing the 'for' loop

}

void uni_solve(double *y,double *f,double *a_x,double *a_y)
{
  int mesh_tracker=0;//we will first start our computations from the finest grid		
  int iter=0;//this will count the no. of iterations at a particular mesh level 	
  int tot_iter=0;//this will count the total no. of iterations at all the mesh levels combined
  double error=0.0,old_error;
  int count=0;//this variable keeps track of the no. of times we change grids 

  printf("At the level=%d\n",mesh_tracker);

  //relaxing on the grid using red-black gauss-seidel (see "relaxation.c")
  for( ; ; ) //I can break out of this loop only under a particular situation 
  {
    old_error=error;//transferring the old error to the current error		
				
    //solving the equations on the current grid
    //see "relaxation.c"
    error=relax(y,f,a_x,a_y,mesh_tracker);//this is a single iteration of the red-black gauss-seidel		

    iter++;//updating the no. of iterations at the current level
    tot_iter++;//updating the total no. of iterations
    printf("tot_iter=%d\titer=%d\terror=%lf\n",tot_iter,iter,error);
    printf("ftoc cond=%lf\told_err=%lf\n",old_error*TOLER_RATE,old_error);	

    //there is a minimum no. of iterations to be performed on every grid
    if(iter<MIN_NUM_ITER)
    {
      continue;      
    }							
    //this is the exit condition
    //if after performing enough iterations on the finest the error is below the tolerance we exit the loop
    else if(error<TOLER)
    {
      //recording the iteration details (see "write_iterations.c")
      write_iter(count,iter,mesh_tracker);		 		
      break;
    }
  }//closing the 'for' loop
}