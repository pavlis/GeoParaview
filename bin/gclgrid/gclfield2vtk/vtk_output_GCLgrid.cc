#include <iostream>
#include <fstream>
#include <string>
#include "gclgrid.h"
#include "GCLMasked.h"
using namespace std;
void write_header_block(string nametag,ofstream& out)
{
    // Header block
    out << "# vtk DataFile Version 2.0"<<endl
    	<< nametag <<endl
    	<< "ASCII" <<endl
    	<< "DATASET POLYDATA"<<endl;
}
void write_points(GCLgrid& g,ofstream& out)
{
    int i,j;
    out << "POINTS " << (g.n1)*(g.n2) <<" float"<<endl;
    for(j=0;j<g.n2;++j)
    	for(i=0;i<g.n1;++i)
    	{
       	    out << g.x1[i][j]
        		<< " "
        		<< g.x2[i][j]
        		<< " "
        		<< g.x3[i][j] 
        		<< endl;
        }
}
int vtk_output_GCLgrid(GCLgrid& g,ofstream& out, string data_tag)
{
	int icell,i,j;
        write_header_block(g.name,out);
        write_points(g,out);
	// Documenation says number of polygons and then
	// "sizeof the cell list".  Since polygons here are 
	// all squares for a full grid this seems to reduce to 5*number polygons
	out << "POLYGONS " << (g.n2 - 1)*(g.n1 - 1) 
		<< " " << (g.n2 - 1)*(g.n1 - 1)*5 <<endl;

	for(j=0,icell=0;j<g.n2-1;++j)
		for(i=0;i<g.n1-1;++i,++icell)
		{
			out << "4 "
				<< j*g.n1 + i << " "
				<< j*g.n1 + i + 1<< " "
				<< (j+1)*g.n1 + i + 1<< " "
				<< (j+1)*g.n1 + i << endl;
		}
        /* Finally the data section.  For the applications this was created
           this surface is always either a DEM or a surface picked by
           interpretation of a horizon in the subsurface.  Hence, we 
           dogmatically force the data field to be elevation. To convert
         to a depth use filters in paraview to flip the polarity*/
		out << "POINT_DATA "<< (g.n1)*(g.n2) <<endl
		<< "SCALARS "<<data_tag<< " double"<<endl
		<< "LOOKUP_TABLE default"<<endl;

	for(j=0;j<g.n2;++j)
		for(i=0;i<g.n1;++i)
		{
		    out <<	-g.depth(i,j) <<endl;
		}
	
	return(icell);
}
/* Polymorphic version to save something other than elevation/depth in
   the data section of the ouput file */
int vtk_output_GCLgrid(GCLscalarfield& g,ofstream& out, string data_tag)
{
	int icell,i,j;
        write_header_block(g.name,out);
        write_points(g,out);
	// Documenation says number of polygons and then
	// "sizeof the cell list".  Since polygons here are 
	// all squares for a full grid this seems to reduce to 5*number polygons
	out << "POLYGONS " << (g.n2 - 1)*(g.n1 - 1) 
		<< " " << (g.n2 - 1)*(g.n1 - 1)*5 <<endl;

	for(j=0,icell=0;j<g.n2-1;++j)
		for(i=0;i<g.n1-1;++i,++icell)
		{
			out << "4 "
				<< j*g.n1 + i << " "
				<< j*g.n1 + i + 1<< " "
				<< (j+1)*g.n1 + i + 1<< " "
				<< (j+1)*g.n1 + i << endl;
		}
        /* Finally the data section.  Here we just use the val array data 
           without any editing*/
		out << "POINT_DATA "<< (g.n1)*(g.n2) <<endl
		<< "SCALARS "<<data_tag<< " double"<<endl
		<< "LOOKUP_TABLE default"<<endl;

	for(j=0;j<g.n2;++j)
		for(i=0;i<g.n1;++i)
		{
		    out <<g.val[i][j] <<endl;
		}
	return(icell);
}
int vtk_output_GCLgrid(GCLMaskedScalarField& g,ofstream& out,string data_tag)
{
    /* This should be like above but does not write polygon data when a point
       is masked.   */
	int icell,i,j;
        write_header_block(g.name,out);
        write_points(g,out);
        /* For a masked grid the number of actual polygons has to be computed.
           To do that we basically build the POLYGONS section first as a list of
           strings.  We then just ask the size of the list to get the size field
           and then dump this list of strings. Algorithm is further complicated 
           by the fact that a valid polygon requires testing all 4 corners.  
           If the masking algorithm were expensive this could be a problem, but
           since this implementation is simply testing a boolean this is not a 
           concern. Point this out as an extension to a generalization could
         make this an issue*/
        list<string> polygons;
        char linebuf[512];
        for(j=0,icell=0;j<g.n2-1;++j)
            for(i=0;i<g.n1-1;++i,++icell)
            {
                if( g.point_is_valid(i,j) && g.point_is_valid(i+1,j)
                        &&  g.point_is_valid(i,j+1) &&  g.point_is_valid(i+1,j+1))
                {
                    /* Because we define all the points above and not just those
                       with valid data we can use the same point indexing formulas
                       as the full grid writer */
                    sprintf(linebuf,"4 %d %d %d %d\n",j*g.n1 + i,j*g.n1 + i + 1,
                            (j+1)*g.n1 + i + 1, (j+1)*g.n1 + i);
                    polygons.push_back(string(linebuf));
                }
            }
        out << "POLYGONS " << polygons.size()<< " " << 5*polygons.size()<<endl;
        list<string>::iterator lptr;
        for(lptr=polygons.begin();lptr!=polygons.end();++lptr)
            out << *lptr<<endl;
        /* Finally write out the data section.   I believe the concept is 
           that values are linked to points.   Any masked points could be 
           deleted, but here we include them and set them to the value 
           defined as a null value below */
        const double masked_val(-9999999.);
		out << "POINT_DATA "<< (g.n1)*(g.n2) <<endl
		<< "SCALARS "<<data_tag<< " double"<<endl
		<< "LOOKUP_TABLE default"<<endl;
	for(j=0;j<g.n2;++j)
		for(i=0;i<g.n1;++i)
		{
		    out <<	g.val[i][j]<<endl;
		}
	
	return(icell);
}
