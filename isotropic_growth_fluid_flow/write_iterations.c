void write_iter(int count,int iter,int mesh_tracker)	
{
	//==============================================================================//
	char fname[300];
	FILE *fp;
	//==============================================================================//


	sprintf(fname,"./DATA/iter.dat");
	if((fp=fopen(fname,"a"))==NULL)
	{
		printf("iteration info can't be written\n");
		exit(1);
	}

	fprintf(fp,"%d\t%d\t%d\n",count,-mesh_tracker,iter);//I am writing down the negatives of mesh_tracker to get the correct appearance of V and W cycles
	fclose(fp);
}
	

