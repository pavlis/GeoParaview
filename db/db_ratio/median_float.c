/* Sort compare function for qsort function */
int compare_float(x1,x2)
float *x1,*x2;
{
        if(*x1<*x2)return (-1);
        else if(*x1==*x2) return (0);
        else return(1);
}


/*  Median for float array function :
   Inputs:  
	y - input vector of length ny for which median is to be calculated
  Returns median of y array values 
*/
float median_float(float *y,int ny)
{

	
	/* This sorts y */
	qsort((char *)y,ny,sizeof(float),compare_float);


	/* These are done exactly using formulas appropriate for small samples*/
	if(ny%2)  /* this is case for odd number */
		return(*(y + ny/2));
	else
		return((y[ny/2 - 1]+y[ny/2])/2.0);

}
