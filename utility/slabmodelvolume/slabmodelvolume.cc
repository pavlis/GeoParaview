/* This is a first order program for the specialized purpose of taking the output
   of the slabmodel program and properly projecting flow lines computed by the 
   program to approximately project the curve a specified depth.  The algorithm is
   an approximation because it does not project out of the plane formed by the 
   vertical and the local tangent to each flow line curve.   In areas with high dip 
   perpendicular to the flow lines this is a poor approximation. */
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include "coords.h"
#include "seispp.h"
#include "GeoTriMeshSurface.h"
#include "GeoSplineSurface.h"
#include "PLGeoPath.h"
#include "PlateBoundaryPath.h"
using namespace std;
using namespace SEISPP;
/* Convenient typedef for primary procedure of this code. */
typedef vector<PLGeoPath> PathArray;
/* Small dangerous helper function to normalize a 3 vector to unit length. 
   Returns normalization constant which should be tested for machine zero */
double unit_vector_normalization(double *x)
{
    /* Note frozen 3 vector assumption - very dangerous but appropriate here */
    double length;
    int i;
    for(length=0.0,i=0;i<3;++i) length+=x[i]*x[i];
    double scale=1.0/sqrt(length);
    for(i=0;i<3;++i) x[i] *= scale;
    return(length);
}
/* This program takes the output of the slabmodel program and 
   builds a surface projected downward to a specified depth 
   (radially) below the input surface.   The normal use for this
   is to build a parallel surface for the crust or lithosphere 
   base consistent with the give top surface.   Default thickness
   argument is 100 km.`
   */
PathArray BuildLithosphereSurface(PathArray& topsurface,double dz)
{
    const string base_error("BuildLithosphereSurface procedure of slabmodelvolume:  ");
    PathArray base;
    int npaths=topsurface.size();
    base.reserve(npaths);
    int i,j,k,jj;
    for(i=0;i<npaths;++i)
    {
        int npts;
        npts=topsurface[i].number_points();
        dmatrix normals(3,npts);
        dmatrix tangents(3,npts);
        Cartesian_point thisnode,lastnode;
        /* First comput tangents using a forward difference.
           Tangent vectors are indexed starting at 0.  Last point is 
           made a copy of the second to last */
        for(j=0,jj=-1;j<(npts-1);++j,++jj)
        {
            thisnode=topsurface[i].node_position_xyz(j);
            if(j==0) 
            {
                lastnode=thisnode;
            }
            else
            {
                tangents(0,jj)=thisnode.x1-lastnode.x1;
                tangents(1,jj)=thisnode.x2-lastnode.x2;
                tangents(2,jj)=thisnode.x3-lastnode.x3;
            }
        }
        for(j=0;j<3;++j) tangents(j,npts-1)=tangents(j,npts-2);
        for(j=0;j<npts;++j)
        {
            /* Need the direction up to compute the normal vector.  
            Note because we are using forward difference this is referenced
            with jj instead of j */
            Geographic_point gpthis=topsurface[i].node_position(j);
            Geographic_point gpdz=gpthis;
            gpdz.r += 10.0;  // Arbitrary radial upward step size.
            Cartesian_point cpr=topsurface[i].coordxyz.cartesian(gpthis);
            Cartesian_point cpdz=topsurface[i].coordxyz.cartesian(gpdz);
            double upvector[3];
            upvector[0]=cpdz.x1-cpr.x1;
            upvector[0]=cpdz.x2-cpr.x2;
            upvector[0]=cpdz.x3-cpr.x3;
            /* now get the normal vector as double cross product */
            double horizontal[3],tangent_unit[3],rawnormal[3];
            dcopy(3,tangents.get_address(0,j),1,tangent_unit,1);
            dr3cros(tangent_unit,upvector,horizontal);
            double norm;
            // normalize this vector so we can test cleanly for parallel vectors
            norm=unit_vector_normalization(horizontal);
            if(norm==0.0) throw GeoCoordError(base_error
                       + "cross product of curve tangent and local vertical is zero\n"
                       + "Cannot handle vertical 90 degree surface dip");
            /* Note this normal points upward because down is a bit
            confusing to work through for me */
            dr3cros(horizontal,tangent_unit,raw_normal);
            norm=unit_vector_normalization(raw_normal);
            /* This number is a conservative measure of equality */
            if(fabs(norm)<0.00000001)
                throw GeoCoordError(base_error
                        + "Computed transverse direction vector has zero length\n"
                        + "This should not happen - use a debugger to find problem");
            /* Replace both tangent and normal with unit vectors */
            dcopy(3,raw_normal,1,normals.get_address(0,j),1);
            dcopy(3,tangent_unit,1,tangents.get_address(0,j),1);
        }
        /* Now project blindly along normals.  Simultaneously
           compute second derivative of tangent.  We need that below to
           avoid overlaps when curvature is convex up */
        dmatrix rawbase(3,npts);
        /* Note this may be misnamed for differential geometry definitions.
           In differential geometry curvature has a cross relationship to the
           path.  What this really here is vector change in tangent vectors */
        dmatrix curvature(3,npts);
        dmatrix tangentbase(3,npts);
        for(j=0;j<npts;++j)
        {
            thisnode=topsurface[i].node_position_xyz(j);
            double nodexyz[3];
            nodexyz[0]=thisnode.x1;
            nodexyz[1]=thisnode.x2;
            nodexyz[2]=thisnode.x3;
            for(k=0;k<3;+k)
                rawbase(k,j)=nodexyz[k]-dz*normals(k,j);
            if(j==npts-1)
            {
                for(k=0;k<3;++k)
                {
                    curvature(k,j)=0.0;
                    tangentbase(k,j)=0.0;
                }
            }
            else
            {
                for(k=0;k<3;++k) 
                {
                    curvature(k,j)=tangents(k,j+1)
                                        - tangents(k,j);
                    /* this is tangent vector of uncorrected basal curve */
                    tangentbase(k,j)=rawbase(k,j+1)-rawbase(k,j);
                }
            }
        }
        /* We only care about the direction of the curvature and tangentbase
           vectors so we normalize them here.  We do this, however, in a 
           special way to handle the common cases when the computed curvature
           is zero.   Less common case to handle is when the tangent base is
           an exact zero vector */
        const double zero_size_test(0.00000001);  // conservative test for double
        for(j=0;j<npts;++j)
        {
            double length;
            length=unit_vector_normalization(curvature.get_address(0,j));
            /* When parallel we set the curvature vector to point upward */
            if(fabs(length)<zero_size_test)
                for(k=0;k<3;++k) curvature(k,j)=-normals(k,j);
            endif
            length=unit_vector_normalization(tangentbase.get_address(0,j));
            /* We force this negative to make sure this is treated like the case
               where the projected points are retrograde */
            if(fabs(length)<zero_size_test)
                for(k=0;k<3;++k) tangentbase(k,j)=-tangents(k,j);
            endif
        }
        dmatrix fixedbase(3,npts);
        /* Now the algorithm we use is this.  when tangent \dot tangentbase is positive
           we assume the curve is not backing up.   When that dot produce it negative we
           scan forward until we find a point whose dot product with the tangent is positive. 
           We then divide the computed vector into that many points and continue */
        j=0;
        int jstart,jend;  // range of indices when handling overlaps
        double basetest[3];
        double frac;   //computed fraction as delta of basetest vector
        while(j<npts)
        {
            if(j==0)
                for(k=0;k<3;++k) fixedbase(k,0)=rawbase(k,0);
            /* This is a test for concave up.   Works because we set normals to point
               downward along flow lines and curvature is turning of tangent vector.
               That vector points up if curving upward and down if curving downward*/
            else if(ddot(3,normals.get_address(0,j),1,curvature.get_address(0,j),1)<=0.0) 
            {
                // When the curvature is concave up assume we can just keep the projected point 
                for(k=0;k<3;++k) fixedbase(k,j)=rawbase(k,j);
            }
            else 
            {
                if(ddot(3,tangents.get_address(0,j),1,tangentbase.get_address(0,j),1)>0.0)
                {
                    for(k=0;k<3;++k) fixedbase(k,j)=rawbase(k,j);
                }
                else
                {
                    /* Here is the hard part - scan forward to find point that does   
                    not overlap this one */
                    jstart=j;
                    do {
                        ++j;
                        if(j==(npts-1))break;
                        for(k=0;k<3;++k) basetest[k]=rawbase(k,j)-rawbase(k,jstart);
                    }while(ddot(3,tangents.get_address(0,jstart),1,basetest,1)<0.0);
                    jend=j;
                    frac=1.0/static_cast<double>(jend-jstart-1);
                    /* WORKING ISSUE:   do not think this handles large
                       slowly changing curves well */ 
                    for(int jj=jstart;jj<=jend;++jj)
                    {
                        for(k=0;k<3;++k)
                            fixedbase(k,jj)=fixedbase(k,jstart)
                                    +basetest[k]*static_cast<double>(jj-jstart)*frac;
                    }
                }
            }
            ++j;
        }
        /* The last thing we need to do is convert the fixedbase points to vector of Geographic_points
           that are then pushed to the base object we are constructing. */
        vector<Geographic_point> fbgeo;
        fbgeo.reserve(npts);
        for(j=0;j<npts;++j)
        {
            Geographic_point ptgeo; 
            Cartesian_point cpt;
            cpt.x1=fixedbase(0,j);
            cpt.x2=fixedbase(1,j);
            cpt.x3=fixedbase(2,j);
            ptgeo=topsurface[i].coordxyz.geographic(cpt);
            fbgeo.push_back(ptgeo);
        }
        PLGeoPath thisbasecurve(fbgeo,0,0.0);
        base.push_back(thisbasecurve);
    }
    return(base);
}
void usage()
{
    cerr << "slabmodelvolume [-t thickness] < in > out "<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
try {
    int i;
    double surface_offset(100.0);
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-t")
        {
            ++i;
            std::istringstream sstr(argv[i]);
            sstr>>surface_offset;
            cerr << "Projecting input surface downward by "
               << surface_offset << " km"<<endl;
        }
        else
            usage();
    }
    typedef vector<Geographic_point> Gvec;
    Gvec thisflowline;
    PathArray topsurface;
    char line[256];
    /* Initialize these variables at top of the main read loop */
    int segsize(0);
    i=0; 
    int nseg(0);
    while(!cin.eof())
    {
        double latin,lonin,zin;
        Geographic_point gpoint;
        cin.getline(line,256);
        if(line[0]=='#') continue;
        if(cin.eof()) break;
        if(line[0]=='>')
        {
            segsize=thisflowline.size();
            if(segsize==0)
                cerr << "Warning:  zero length line segment number "<<nseg
                    <<endl << "Ignored."<<endl;
            else
            {
                PLGeoPath PLthisfl(thisflowline,0,0.0);
                topsurface.push_back(PLthisfl);
                thisflowline.clear();
                ++nseg;
            }
        }
        else
        {
            sscanf(line,"%lf%lf%lf",&lonin,&latin,&zin);
	    gpoint.lat=rad(latin);
	    gpoint.lon=rad(lonin);
	    gpoint.r=r0_ellipse(gpoint.lat)-zin;
  	    thisflowline.push_back(gpoint);
        }
    }
    // Data are now loaded.  Now the for the projection
    PathArray lowersurface=BuildLithosphereSurface(topsurface,surface_offset);
    /* Now write the new surface to stdout */
    PathArray::iterator pathptr;
    int j;
    for(pathptr=lowersurface.begin(),j=0;pathptr!=lowersurface.end();++pathptr,++j)
    {
    	int npts=pathptr->number_points();
    	for(i=0;i<npts;++i)
    	{
    		Geographic_point nodepoint=pathptr->node_position(i);
    		cout << deg(nodepoint.lon) << " "
    		     << deg(nodepoint.lat) << " "
    		     << r0_ellipse(nodepoint.lat)-nodepoint.r<<endl;
    	}
    	cout << "<"<<endl;
    }
} catch(exception& err)
{
    cerr << "Fatal error:  "<<endl
        << err.what()<<endl;
    exit(-1);
}
}

