// Test program for basic functionality of gclgrid library
#include <iostream>
#include <fstream>
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
void fill_pattern(GCLgrid *g)
{
	for(int i=0;i<g->n1;++i)
		for(int j=0;j<g->n2;++j)
		{
			g->x1[i][j]=i;
			g->x2[i][j]=j;
			g->x3[i][j]=i+j;
		}
}
void fill_pattern(GCLgrid3d *g)
{
	for(int i=0;i<g->n1;++i)
		for(int j=0;j<g->n2;++j)
		for(int k=0;k<g->n3;++k)
		{
			g->x1[i][j][k]=i;
			g->x2[i][j][k]=j;
			g->x3[i][j][k]=k;
		}
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

	
	g2d=new GCLgrid(3,3);
	g3d=new GCLgrid3d(3,3,2);
	s2d=new GCLscalarfield(3,3);
	v2d=new GCLvectorfield(3,3,2);
	s3d=new GCLscalarfield3d(3,3,2);
	v3d=new GCLvectorfield3d(3,3,2,3);
	fill_pattern(dynamic_cast<GCLgrid*>(s2d));
	fill_pattern(dynamic_cast<GCLgrid*>(v2d));
	fill_pattern(dynamic_cast<GCLgrid3d*>(s3d));
	fill_pattern(dynamic_cast<GCLgrid3d*>(v3d));
// temp
//cout << *s2d;
//cout << *v2d;
//cout << *s3d;
//cout << *v3d;

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
	g3d=new GCLgrid3d(20,20,10,name,
		rad(9.0),rad(-63.0),6371.0,
		rad(100.0),50.0,50.0,30.0,
		10,10);
	name="testgrid2";
	GCLgrid3d *grid2=new GCLgrid3d(*g3d);
	//strcpy(grid2->name,name.c_str());
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
	cout << "Saving results to testgrid1.out" << endl;
	ofstream out1;
	out1.open("testgrid1.out",ios::out);

	out1 << testgrid1;
	out1.close();

	cout << "Creating test grids for summation test" << endl;
	cout << " test grid 2.  lat lon offset test";


	name="testgrid1";
	g3d=new GCLgrid3d(20,20,10,name,
		rad(9.0),rad(-63.0),6371.0,
		rad(100.0),50.0,50.0,30.0,
		10,10);
	name="testgrid3";
	grid2=new GCLgrid3d(20,20,10,name,
		rad(10.0),rad(-64.0),6371.0,
		rad(110.0),50.0,50.0,30.0,
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
	cout << "testing assignment operator for scalar field" << endl;
	GCLscalarfield3d testgrid5=testgrid4;
	cout << testgrid5;
	cout << "testing == operator" << endl;
	if(testgrid4==testgrid5)
		cout << "Test passed.  equal grids yield equal"<<endl;
	else
		cout << "Test failed.  equal grids resolve no equal"<<endl;
	cout << "testing remap functions" << endl;
	GCLgrid3d remaptest=remap_grid(dynamic_cast<GCLgrid3d&>(testgrid3),
		dynamic_cast<GCLgrid3d&>(testgrid4));
	GCLscalarfield3d rmapfield(remaptest);
	cout << rmapfield << endl;
}
	

