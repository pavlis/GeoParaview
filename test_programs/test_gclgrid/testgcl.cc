// Test program for basic functionality of gclgrid library
#include <iostream>
#include "stock.h"
#include "pf.h"
#include "dmatrix.h"
#include "gclgrid.h"

void load_constant(GCLscalarfield3d& g,double c)
{
	int i,j,k;
	for(i=0;i<g.n1;++i)
		for(j=0;j<g.n2;++j)
			for(k=0;k<g.n3;++k)
				g.val[i][j][k]=c;
}


int main(int arc, char **argv)
{
	cout << "testing constructors with no arguments" << endl;
	GCLgrid *g2d;
	GCLgrid3d *g3d;
	GCLscalarfield *s2d;
	GCLvectorfield *v2d;
	GCLscalarfield3d *s3d;
	GCLvectorfield3d *v3d;
	
	g2d=new GCLgrid();
	g3d=new GCLgrid3d();
	s2d=new GCLscalarfield();
	v2d=new GCLvectorfield();
	s3d=new GCLscalarfield3d();
	v3d=new GCLvectorfield3d();

	cout << "construction successful.  Destroying" << endl;

	delete g2d;
	delete g3d;
	delete s2d;
	delete v2d;
	delete s3d;
	delete v3d;

	cout << "Trying integer size constructors " << endl;

	
	g2d=new GCLgrid(100,50);
	g3d=new GCLgrid3d(100,50,20);
	s2d=new GCLscalarfield(100,50);
	v2d=new GCLvectorfield(100,50,5);
	s3d=new GCLscalarfield3d(100,50,20);
	v3d=new GCLvectorfield3d(100,50,20,3);

	cout << "construction successful.  Destroying" << endl;

	delete g2d;
	delete g3d;
	delete s2d;
	delete v2d;
	delete s3d;
	delete v3d;

	// dbroutines are tested in makegclgrid so we don't mess with them here
	// some of the lower level items like lat,lon, gtoc, etc. are not tested
	// directly as they are used implicitly in some higher level member functions
	// we test here

	cout << "Creating test grids for summation test" << endl;
	cout << " test grid 1.  simple exact match test";


	string name("testgrid1");
	g3d=new GCLgrid3d(20,20,10,const_cast<char *>(name.c_str()),
		9.0,-63.0,6371.0,
		100.0,50.0,50.0,30.0,
		10,10);
	name="testgrid2";
	GCLgrid3d *grid2=new GCLgrid3d(*g3d);
	strcpy(grid2->name,name.c_str());
	cout << "grids created.  Creating 3d scalar field" << endl;
	GCLscalarfield3d testgrid1(*g3d);
	GCLscalarfield3d testgrid2(*grid2);
	delete g3d;
	delete grid2;
	cout << "loading test pattern 1 into test grids"<<endl;
	int i,j,k,l;

	load_constant(testgrid1,1.0);
	load_constant(testgrid2,1.0);

	cout << "testing += operator" << endl;

	testgrid1 += testgrid2;

	cout << testgrid1;

	cout << "Creating test grids for summation test" << endl;
	cout << " test grid 2.  lat lon offset test";


	name="testgrid1";
	g3d=new GCLgrid3d(20,20,10,const_cast<char *>(name.c_str()),
		9.0,-63.0,6371.0,
		100.0,50.0,50.0,30.0,
		10,10);
	name="testgrid3";
	grid2=new GCLgrid3d(20,20,10,const_cast<char *>(name.c_str()),
		10.0,-64.0,6371.0,
		100.0,50.0,50.0,30.0,
		10,10);
	cout << "grids created.  Creating 3d scalar field" << endl;
	GCLscalarfield3d testgrid3(*g3d);
	GCLscalarfield3d testgrid4(*grid2);
	delete g3d;
	delete grid2;
	cout << "loading test pattern 1 into test grids"<<endl;

	load_constant(testgrid3,1.0);
	load_constant(testgrid4,1.0);

	cout << "testing += operator" << endl;

	testgrid3 += testgrid4;

	cout << testgrid3;


}
	

