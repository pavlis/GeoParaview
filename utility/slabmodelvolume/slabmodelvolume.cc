/* This is a first order program for the specialized purpose of taking the output
   of the slabmodel program and properly projecting flow lines computed by the 
   program to approximately project the curve a specified depth.  The algorithm is
   an approximation because it does not project out of the plane formed by the 
   vertical and the local tangent to each flow line curve.   In areas with high dip 
   perpendicular to the flow lines this is a poor approximation. */
#include <float.h>
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
typedef vector<dmatrix> MatrixVector;
/* Small dangerous helper function to normalize a 3 vector to unit length. 
   Returns normalization constant which should be tested for machine zero */
double unit_vector_normalization(double *x)
{
    /* Note frozen 3 vector assumption - very dangerous but appropriate here */
    double length;
    int i;
    for(length=0.0,i=0;i<3;++i) length+=x[i]*x[i];
    // Do nothing if length is 0
    if(length<DBL_EPSILON) return(0.0);
    length=sqrt(length);
    double scale=1.0/length;
    for(i=0;i<3;++i) x[i] *= scale;
    return(length);
}
MatrixVector smooth_normals(MatrixVector& n)
{
    /* We use a 11 point smoother.  5 along flow lines and only 3 points 
         on each cross line.  Equal weight with that geometry is 1/11.   For 
         edges it is 1/8 */
    const double wt(1.0/9.0);
    const double ewt(1.0/6.0);
    MatrixVector result;
    int i,j,k;
    MatrixVector::iterator nptr,nptrp1;
    dmatrix n0;
    for(nptr=n.begin(),i=0;nptr!=n.end();++nptr,++i)
    {
        dmatrix nm1,np1;
        if(i!=0)
            nm1=n0;
        n0=(*nptr);
        int npts=n0.columns();
        if(i<(npts-1))
        {
            nptrp1=nptr+1;
            np1=(*nptrp1);
        }

        // Intentionally make output copy.  This allow us to 
        // simply do nothing to some boundary points
        dmatrix snormals(n0); 
        /* along flow line section does not do an end correction
           This means first 2 and last 2 points are not smoothed. */
        int nend;
        int jj;
        nend=n0.columns();
        nend=min(nend,nm1.columns());
        nend=min(nend,np1.columns());
        nend-=2;
        for(j=2;j<nend;++j)
        {
            double sum;
            if( (i==0) || (i==(npts-1)))
            {
                for(k=0;k<3;++k)
                {
                    sum=0.0;
                    for(jj=j-1;jj<=j+1;++jj)
                    {
                        if(i!=0)
                            sum+=wt*nm1(k,jj);
                        else
                            sum+=wt*np1(k,jj);
                    }
                    for(jj=j-2;jj<=j+2;++jj)
                        sum+=wt*n0(k,jj);
                    snormals(k,j)=sum;
                }
            }
            else
            {
                for(k=0;k<3;++k)
                {
                    sum=0.0;
                    for(jj=j-1;jj<=j+1;++jj)
                    {
                        sum+=wt*np1(k,jj);
                        sum+=wt*nm1(k,jj);
                    }
                    for(jj=j-2;jj<=j+2;++jj)
                        sum+=wt*n0(k,jj);
                    snormals(k,j)=sum;
                }
            }
            double length=unit_vector_normalization(snormals.get_address(0,j));
        }
        result.push_back(snormals);
    }
    return(result);
}
/* This program takes the output of the slabmodel program and 
   builds a surface projected downward to a specified depth 
   (radially) below the input surface.   The normal use for this
   is to build a parallel surface for the crust or lithosphere 
   base consistent with the give top surface.   Default thickness
   argument is 100 km.`

   If smooth is set true a (frozen size) smoother is applied to the
   normal vectors before projection.
   */
PathArray BuildLithosphereSurface(PathArray& topsurface,double dz,
        bool smooth)
{
    const string base_error("BuildLithosphereSurface procedure of slabmodelvolume:  ");
    PathArray base;
    int npaths,npts;
    npaths=topsurface.size();
    base.reserve(npaths);
    int i,j,k;
    MatrixVector tangents,crosslines,normals,ds;
    Cartesian_point thisnode,lastnode,nextnode;
    double length;
    for(i=0;i<npaths;++i)
    {
        int npts;
        npts=topsurface[i].number_points();
        dmatrix pathtangents(3,npts);
        dmatrix pathds(1,npts);
        /* First comput tangents using a forward difference.
           Tangent vectors are indexed starting at 0.  Last point is 
           made a copy of the second to last */
//DEBUG
  //      cerr << "Start tangents for flow line "<<i<<endl;
        for(j=0;j<npts;++j)
        {
            thisnode=topsurface[i].node_position_xyz(j);
            //DEBUG
            //cerr << "x,y,z="<<thisnode.x1<<" "<<thisnode.x2<<" "<<thisnode.x3<<endl;
            if(j==0) 
            {
                lastnode=thisnode;
                pathds(0,j)=0.0;
            }
            else
            {
                pathtangents(0,j-1)=thisnode.x1-lastnode.x1;
                pathtangents(1,j-1)=thisnode.x2-lastnode.x2;
                pathtangents(2,j-1)=thisnode.x3-lastnode.x3;
                length=unit_vector_normalization(pathtangents.get_address(0,j-1));
                pathds(0,j)=length;
                //DEBUG
                //cerr << "i j distance:  "<<i<<" "<<j<<" "<<length<<endl;
                /* The size of this test is dependent upon assumed cartesian
                   coordinates in units of km.  i.e. this means within 1 m */
                if(fabs(length)<0.001) throw GeoCoordError(base_error
                        + "Duplicate points in input.  Fix data file");
                lastnode=thisnode;
            }
            //DEBUG
            //cerr << "tangent " << j-1 << " length="<<length<<endl;
        }
        /* Above loop leaves last point unset.  Make it a copy of the second to last */
        for(j=0;j<3;++j) pathtangents(j,npts-1)=pathtangents(j,npts-2);
        tangents.push_back(pathtangents);
        ds.push_back(pathds);
    }
    /* Now get the cross-flowline dip direction by computing tangent vectors
       in the crossline direction of the flow lines.  This is complicated by 
       the fact that the flow lines will (for most recent versions) have irregular
       lengths. Because this is a locked sibling to slabmodel we can assume that 
       all points start at the trench (0) and the differences are only how long
       each curve is.   The algorithm is to compute crossline dips by finite differences
       in the + mesh direction and if path 2 (more +) is shorter to use the last 
       valid crossdip vector to the end of the current path.  The npts-1 path is
       a copy of the npts-2 path crossdip.
    
    A Huge complication was created by using PLGeoPath objects.  Each path in
    these objects define it's own local coordinate system with the origin at 
    the first point in the path.  What we do here handle this in forward 
    difference schemes to compute the crossline vectors.  They will be correct
    in each line segments coordinates then.  */
    Geographic_point gpwork;
    /* We need this outside the loop to simplify copying of last path crossline vectors */
    dmatrix pathcross;
    for(i=0;i<npaths-1;++i)
    {
        int np1=topsurface[i].number_points();
        int np2=topsurface[i+1].number_points();
        int np=min(np1,np2);
        pathcross=dmatrix(3,np1);
//DEBUG
  //      cerr << "Start crossline computation for flow line "<<i<<endl;
        for(j=0;j<np;++j)
        {
            gpwork=topsurface[i+1].node_position(j);
            nextnode=topsurface[i].coordxyz.cartesian(gpwork);
            thisnode=topsurface[i].node_position_xyz(j);

            pathcross(0,j)=nextnode.x1-thisnode.x1;
            pathcross(1,j)=nextnode.x2-thisnode.x2;
            pathcross(2,j)=nextnode.x3-thisnode.x3;
            length=unit_vector_normalization(pathcross.get_address(0,j));
            if(fabs(length)<0.001) throw GeoCoordError(base_error
                    + "Intersecting flow line curves found. Increase crossline spacing");
            //DEBUG
            //cerr << "node "<<j<<" length "<<length<<endl;
        }
        for(j=np;j<np1;++j)
        {
            for(k=0;k<3;++k) pathcross(k,j)=pathcross(k,np-1);
        }
        crosslines.push_back(pathcross);
    }
    crosslines.push_back(pathcross);
    /* Now we compute normals for each node point from cross product of 
       tangent and crossline.   There is a nontrivial sign issue here.  
       The generalized coordinates in the slabmodel grid can be left or
       right handed.  We do this by computing a vector at each node point
       that points up and using the sign of the dot product to assume our
       normal here point downward (the direction we always project here). */
    MatrixVector::iterator tptr,cptr;  //iterators for tangent and crossline resp
    //DEBUG remove 
    //dmatrix testmat;
    for(i=0,tptr=tangents.begin(),cptr=crosslines.begin();i<npaths;
            ++i,++tptr,++cptr)
    {
        //DEBUG
        /*cerr <<"Tangents for path "<<i<<endl
            << *tptr <<endl;
        cerr<< "Crossline vectors for same path"<<endl
            << *cptr <<endl;
            */
        npts=topsurface[i].number_points();
        dmatrix pathnormals(3,npts);
        for(j=0;j<npts;++j)
        {
            Geographic_point gpthis=topsurface[i].node_position(j);
            Geographic_point gpdz=gpthis;
            gpdz.r += 10.0;  // Arbitrary radial upward step size.
            Cartesian_point cpr=topsurface[i].coordxyz.cartesian(gpthis);
            Cartesian_point cpdz=topsurface[i].coordxyz.cartesian(gpdz);
            double upvector[3];
            upvector[0]=cpdz.x1-cpr.x1;
            upvector[1]=cpdz.x2-cpr.x2;
            upvector[2]=cpdz.x3-cpr.x3;
            /* now get the normal vector as cross product of tangent and crossline 
               vectors */
            double rawnormal[3];
            dr3cros(tptr->get_address(0,j),cptr->get_address(0,j),rawnormal);
            length=unit_vector_normalization(rawnormal);
            /*This size is a bit arbitrary (about float epslison).   
              Likely conservative.  Intentionally made to not be a
              fatal error but the warning seems appropriate */
            if(length<0.000001) 
            {
                cerr << base_error 
                    << "WARNING:  cross product of inline and crossline tangent "
                    << "vectors used to compute projection direction is very small"
                    <<endl
                    << "Length computed is "<<length<<endl
                    << "Working on flow line number "<<i <<" at node "<<j<<endl
                    << "Check for geometry anomalies at this point"<<endl;;
            }
            /* Now force the normal to be downward in the earth */
            double testdot=ddot(3,rawnormal,1,upvector,1);
            //DEBUG
            //cerr << j << " "<<testdot<<endl;
            /* Again the test for zero is a bit arbitrary and a bit conservative */
            if(fabs(testdot)<0.0000001) 
            {
                cerr << base_error
                    << "WARNING:  near vertical dip detected."<<endl
                    << "Computed normal vector is very close to horizonal. "
                    << "Dot product of up direction and tangent unit vectors is only "
                    << length<<endl
                    << "Working on flow line number "<<i <<" at node "<<j<<endl
                    << "This may create geometry anomalies in the output"
                    <<endl;
            }
            /* This assures normal used points down */
            if(testdot>0.0)
                for(k=0;k<3;++k) pathnormals(k,j)=-rawnormal[k];
            else
                for(k=0;k<3;++k) pathnormals(k,j)=rawnormal[k];
        }
        normals.push_back(pathnormals);
        //DEBUG
        //testmat=tr(pathnormals);
        //cerr << "normals for path "<<i<<endl;
        //cerr << testmat<<endl; 
    }
    if(smooth) normals=smooth_normals(normals);
    /* Now project along normals.   In the same loop we compute sin of the angle
       between tangents which is used to test for overlaps when the curvature is strongly
       concave down. */
    MatrixVector::iterator nptr,dsptr;
    for(i=0,tptr=tangents.begin(),nptr=normals.begin(),dsptr=ds.begin();
            i<npaths;++i,++tptr,++nptr,++dsptr)
    {
        npts=topsurface[i].number_points();
        dmatrix rawbase(3,npts);
        /* This vector is used as a test to deal with retrograde line segments
           in projecting lines with strong concave down curvature. */
        dvector xtest(npts); 
        /* It is a bit inefficient to make these copies of path vectors, but
           I do not expect a performance problem for these simple calculations
           and this makes the syntax SOOO much cleaner. (alternative is
           ugly thins like uptr->(k,j)  */
        dmatrix pathtangents(*tptr);
        dmatrix pathnormals(*nptr);
        dmatrix  pathds(*dsptr);
        /* Now project this path */
        for(j=0;j<npts;++j)
        {
            double tangent_dot_product,phi;
            double dtds[3],nudotdtds;
            thisnode=topsurface[i].node_position_xyz(j);
            rawbase(0,j)=thisnode.x1+dz*pathnormals(0,j);
            rawbase(1,j)=thisnode.x2+dz*pathnormals(1,j);
            rawbase(2,j)=thisnode.x3+dz*pathnormals(2,j);
            if(j==npts-1)
            {
                xtest(j)=0.0;
            }
            else
            {
                tangent_dot_product=ddot(3,pathtangents.get_address(0,j),1,
                        pathtangents.get_address(0,j+1),1);
                /* This test seems necessary.  Found occasional nans otherwise
                   from roundoff errors.  Using conservative 32 bit epsilon */
                if(fabs(tangent_dot_product-1.0)<FLT_EPSILON)
                    phi=0.0;
                else
                    phi=acos(tangent_dot_product);
                xtest(j)=dz*sin(phi);
                /* Set xtest negative (which forces it to be ignored) for 
                   the case of concave up.  Determined by dot product of 
                   tangent vector difference with normal.  Sign is 
                   because normals point down */
                for(k=0;k<3;++k) 
                    dtds[k]=pathtangents(k,j+1)-pathtangents(k,j);
                nudotdtds=ddot(3,dtds,1,pathnormals.get_address(0,j),1);
                if(nudotdtds<0) xtest(j) = (-xtest(j));
                //DEBUG
                //cerr << "xtest for "<<j<<" is "<<xtest(j)<<endl;
            }
        }
        //DEBUG
        /*
        if(i==8) {
        testmat=tr(rawbase);
        cerr << "rawbase"<<endl<<testmat<<endl;
        testmat=tr(pathnormals);
        cerr << "pathnormals"<<endl<<testmat<<endl;
        testmat=tr(pathtangents);
        cerr << "pathtangents"<<endl<<testmat<<endl;
        testmat=tr(pathds);
        cerr << "pathds"<<endl<<testmat<<endl;
        }
        */
        dmatrix fixedbase(3,npts);
        /* We use the geometry defined by the xtest vector to decide if we need to deal
           with the complications of a strongly concave down curve segment.   The 
           theory behind this follows from simple trigonometry with a right triangle 
           with angles phi and 90-phi with phi computed from local turning rate of the 
           curve.   Note xtest can and will be negative if the curve is turning upward 
           so the test is simply for ds values larger than xtest. */
        j=0;
        int jstart,jend;  // range of indices when handling overlaps
        double basetest[3];
        double frac;   //computed fraction as delta of basetest vector
        for(j=0;j<npts;++j)
        {
            double dssum,sumxtest;
            if(j==0)
                for(k=0;k<3;++k) fixedbase(k,0)=rawbase(k,0);
            else if(pathds(0,j)>=xtest(j))
            {
                /* Case with no overlap - just copy rawbase to fixedbase*/
                for(k=0;k<3;++k) fixedbase(k,j)=rawbase(k,j);
            }
            else 
            {
                jstart=j;
                dssum=0.0;
                do{
                    double phi;
                    ++j;
                    if(j==npts)break;
                    dssum+=pathds(0,j);
                    double tangent_dot_product=ddot(3,pathtangents.get_address(0,jstart),1,
                            pathtangents.get_address(0,j),1);
                    if(fabs(tangent_dot_product-1.0)<FLT_EPSILON)
                        phi=0.0;
                    else
                        phi=acos(tangent_dot_product);
                    sumxtest=dz*sin(acos(phi));
                }while(dssum<sumxtest);
                if(j==npts)
                {
                        /* A crude recovery for this condition.  This is entered if
                           the end of the curve is retrograde.   What we do then is 
                           project along the last valid tangent a tiny distance 
                           where tiny is a scale dependent for this application */
                        frac=1.0; /* Since tangents are unit vectors this means 100 m */
                        double increment;
                        for(int jfix=jstart,increment=frac;jfix<npts;++jfix,increment+=frac)
                        {
                            for(k=0;k<3;++k) 
                            {
                                fixedbase(k,jfix)=fixedbase(k,jstart-1)
                                                + increment*(pathtangents(k,jstart-1));
                            }
                        }
                }
                else
                {
                        jend=j;
                        // set basetest to vector between jstart and jend
                        for(k=0;k<3;++k) basetest[k]=rawbase(k,jend)-rawbase(k,jstart);
                        frac=1.0/static_cast<double>(jend-jstart-1);
                        for(int jj=jstart;jj<=jend;++jj)
                        {
                            for(k=0;k<3;++k)
                                fixedbase(k,jj)=fixedbase(k,jstart-1)
                                        +basetest[k]*static_cast<double>(jj-jstart)*frac;
                        }
                }
            }
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
    cerr << "slabmodelvolume [-t thickness -s] < in > out "<<endl
        << "  where -s enables smoothing option (off by default)"<<endl;;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
try {
    int i;
    double surface_offset(100.0);
    bool smooth(false);
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
        else if(sarg=="-s")
            smooth=true;
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
    PathArray lowersurface=BuildLithosphereSurface(topsurface,
            surface_offset,smooth);
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
    	cout << ">"<<endl;
    }
} catch(exception& err)
{
    cerr << "Fatal error:  "<<endl
        << err.what()<<endl;
    exit(-1);
}
}

