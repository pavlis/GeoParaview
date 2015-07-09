#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

/* Sort compare function for qsort function */
int compare_float(x1,x2)
float *x1,*x2;
{
        if(*x1<*x2)return (-1);
        else if(*x1==*x2) return (0);
        else return(1);
}

/* Main subroutine like function to return statistical parameters 
calculated from a vector of floats */

void calc_statistics(y,ny,mean,median,low_quartile,high_quartile,ymin,ymax)
float *y;
int ny;
float *mean,*median,*low_quartile,*high_quartile,*ymin,*ymax;

/* Inputs:  
	y - input vector of length ny for which statistics are to 
		be calculated
   Outputs:
	mean,median - as name implies
	low_quartile, high_quartile - 1/4 and 3/4 quartile points
	ymin,ymax - low and high value of y respectively.

Algorithm calculates median, quartiles, min, and max by sorting input
y array.  Note this is destructive and the y array is rearranged on
return.  mean is calculate straight up by the normal formula.
*/
{
	int i;  /* standard counter */

	/* first let's go ahead and calculate the mean */
	*mean = 0.0;
	for(i=0;i<ny;++i) *mean += *(y+i);
	*mean /= (float)ny;
	
	/* Now go ahead and sort y using the qsort function */
	qsort((char *)y,ny,sizeof(float),compare_float);


	/* These are done exactly using formulas appropriate for small samples*/
	if(ny%2)  /* this is case for odd number */
		*median = *(y + ny/2);
	else
		*median = (y[ny/2 - 1]+y[ny/2])/2.0;
	switch((ny-1)%4)
	{
	case(0):
		*low_quartile = y[(ny-1)/4];
		*high_quartile = y[3*(ny-1)/4];
		break;
	case(1):
		*low_quartile = (0.75*y[(ny-1)/4]+0.25*y[((ny-1)/4)+1]);
		*high_quartile = (0.75*y[3*(ny-1)/4]+0.25*y[(3*(ny-1)/4)+1]);
		break;
	case(2):
		*low_quartile = (0.5*y[(ny-1)/4]+0.5*y[((ny-1)/4)+1]);
		*high_quartile = (0.5*y[3*(ny-1)/4]+0.5*y[(3*(ny-1)/4)+1]);
		break;
	case(3):
		*low_quartile = (0.25*y[(ny-1)/4]+0.75*y[((ny-1)/4)+1]);
		*high_quartile = (0.25*y[3*(ny-1)/4]+0.75*y[(3*(ny-1)/4)+1]);
		break;

	}
	*ymin = *y;
	*ymax = *(y+ny-1);
}
