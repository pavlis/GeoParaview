#include <stdio.h>
#include <string.h>
#include "stock.h"
#include "arrays.h"
#include "pf.h"

/* compare function for pure integer keys for associative array*/
int cmpint(int *i1,int *i2)
{
	if(i1<i2)
		return(-1);
	else if(i1==i2)
		return(0);
	else
		return(1);
}
/* This function creates an associative array keyed by an integer
event id from a parameter file based descriptor &Tbl like this example:

pmel_calibration_events &Tbl{
	10 xyz
	11 xyzt
}
where 10 and 11 are event ids and the string defines coordinates to fix.

Author:  Gary Pavlis
*/
Arr *load_calibration_events(Pf *pf)
{
	Tbl *t;
	int evid;
	char fix[5],*line;
	Arr *a;
	int i;

	a = newarr(cmpint);

	t = pfget_tbl(pf,"pmel_calibration_events");
	if(t==NULL) 
	{
		elog_complain(0,
		  "pmel_calibration_events Tbl not in parameter file\nAssuming no calibration events exist\n");
	}
	else
	{
		for(i=0;i<maxtbl(t);++i)
		{
			line = (char *)gettbl(t,i);
			sscanf(line,"%d%s",&evid,fix);
			setarr(a,(char *)evid,fix);
		}
			
	}
	return(a);
}
char *get_fixlist(Arr *a,int evid)
{
	char *o;
	o=(char *)getarr(a,(char *)evid);
	return(o);
}
