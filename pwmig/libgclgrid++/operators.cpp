#include "gclgrid.h"
GCLscalarfield3d GCLscalarfield3d::operator+=(const GCLscalarfield3d& g)
{
	int i,j,k;
	double valnew;

	reset_index(); 

	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			for(k=0;k<n3;++j)
			{
				try{
					lookup(g.x1[i][j][k],g.x2[i][j][k],g.x3[i][j][k]);
				}
				catch(int err) { GCLlookup_error_handler(err); };
				valnew=interpolate(g.x1[i][j][k],g.x2[i][j][k],g.x3[i][j][k]);

				val[i][j][k]+=valnew;
			}
		}
	}
	return(*this);
}

GCLvectorfield3d GCLvectorfield3d::operator += (const GCLvectorfield3d& g)
{
	int i,j,k,l;
	double *valnew;

	reset_index(); 

	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			for(k=0;k<n3;++j)
			{
				try{
					lookup(g.x1[i][j][k],g.x2[i][j][k],g.x3[i][j][k]);
				}
				catch(int err) { GCLlookup_error_handler(err); };
				valnew=interpolate(g.x1[i][j][k],g.x2[i][j][k],g.x3[i][j][k]);
				for(l=0;l<nv;++nv) val[i][j][k][l]=valnew[l];
				delete valnew;
			}
		}
	}
	return(*this);
}
GCLscalarfield GCLscalarfield::operator += (const GCLscalarfield& g)
{
	int i,j,k;
	double valnew;

	reset_index(); 

	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			try{
				lookup(g.lat[i][j],g.lon[i][j]);
			}
			catch(int err) { GCLlookup_error_handler(err); };
			valnew=interpolate(g.x1[i][j],g.x2[i][j],g.x3[i][j]);
			val[i][j]+=valnew;
		}
	}
	return(*this);
}

GCLvectorfield GCLvectorfield::operator += (const GCLvectorfield& g)
{
	int i,j,k,l;
	double *valnew;

	reset_index(); 

	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			try{
				lookup(g.lat[i][j],g.lon[i][j]);
			}
			catch(int err) { GCLlookup_error_handler(err); };
			valnew=interpolate(g.x1[i][j],g.x2[i][j],g.x3[i][j]);
			for(l=0;l<nv;++nv) val[i][j][l]=valnew[l];
			delete valnew;
		}
	}
	return(*this);
}
GCLscalarfield3d GCLscalarfield3d::operator *= (double c1)
{
	int i,j,k;
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k)
				val[i][j][k]*=c1;
	return(*this);
}
GCLvectorfield3d GCLvectorfield3d::operator *= (double c1)
{
	int i,j,k,l;
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k)
				for(l=0;l<nv;++l)
					val[i][j][k][l]*=c1;
	return(*this);
}
GCLscalarfield GCLscalarfield::operator *= (double c1)
{
	int i,j,k;
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			val[i][j]*=c1;
	return(*this);
}
GCLvectorfield GCLvectorfield::operator *= (double c1)
{
	int i,j,k,l;
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(l=0;l<nv;++l)
				val[i][j][l]*=c1;
	return(*this);
}
