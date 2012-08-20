#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include "stock.h"
#include "pf.h"
#include "seispp.h"
#include "Metadata.h"
#include "gclgrid.h"
using namespace SEISPP;
void usage()
{
	cerr << "rpitomo2db tomofile db modelname [-S -vpvs x -savegrid]\n"
                << "tomofile is input"<<endl
                << "db is Antelope database name" <<endl
                << "modelname is base name (_P and/or _S are appended to this name"<<endl
		<< "Options: "<<endl
		<< "-S save the S model (default is P)"<<endl
		<< "-vpvs save a version of P model converted to S with vpvs"
                << "-savegrid saves grid file (default does not - normally run once per tomography grid definition)"
		<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	//int nx,ny,nz;
	FILE *fp;
	string fname;
	int ntotal,ntest;
	double lat,lon;
	//double dx,dy,dz;
	//double dxutm,dyutm;
	char *tomofile,*dbname;
	double vpvs=0.707;
	bool save_S=false;
        bool save_grid(false);

	ios::sync_with_stdio();
	if(argc<3) usage();
	tomofile=argv[1];
	dbname=argv[2];
        string modelbasename(argv[3]);
	bool usevpvs(false);
	for(int iarg=4;iarg<argc;++iarg)
	{
		string sarg(argv[iarg]);
		if(sarg=="-vpvs")
		{
			++iarg;
			vpvs=atof(argv[iarg]);
			usevpvs=true;
		}
		if(sarg=="-S")
			save_S=true;
                if(sarg=="-savegrid")
                    save_grid=true;
		else
			usage();
	}
	if(usevpvs && save_S) 
	{
		cerr << "Inconsistent parameters."
			<<endl
			<< "-vpvs and -S cannot be used together"<<endl;
		usage();
	}
        cout << "Program "<<argv[0]<<" starting"<<endl
            << "This is version for newer format grid files"<<endl;
	
	Pf *pf;
        string pfname("rpitomo2db");
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
        {
            cerr << "pfread failed for pf file="<<pfname<<endl;
            exit(-1);
        }
	Dbptr db;
	if(dbopen(dbname,"r+",&db)==dbINVALID)
	{
		cerr << "Cannot open output database "<<dbname<<endl;
		exit(-1);
	}
	
	try {
		int i,j,k,kk;
		Metadata control(pf);
		double lat0=control.get_double("grid3d_lat0");
		double lon0=control.get_double("grid3d_lon0");
		double r0=control.get_double("grid3d_r0");
                double gridaz=control.get_double("grid3d_azimuth");
		string gridname=control.get_string("output_gridname_base");
		string fielddir=control.get_string("field_directory");
		string gclgdir=control.get_string("grid_directory");
                // Hard code this.  Thought about putting it in pf
                // but since this is locked to rpi code make it const
                const int headersize(232);
		// convert these to radians
		lat0=rad(lat0);
		lon0=rad(lon0);
		gridaz=rad(gridaz);
		FILE *fp=fopen(tomofile,"r");
		if(fp==NULL)
		{
			cerr << "Open of file "<<tomofile
				<<"failed.  Exiting\n";
			exit(-1);
		}
                /* Read and crack the header */
                char *head=new char[headersize];
                fread(head,1,headersize,fp);
                // This should be removed once the SPHR and grid type test validated
                cout <<"Text section of heder read from file"<<endl;
                const int hdumpsize(60);
                for(i=0;i<hdumpsize;++i)cout << head[i];
                cout << endl;
                /* Important to test that this is a spherical grid.
                   Cartesian grid in Steve's code is a different beast.*/
                char gtest[5];
                for(i=0;i<4;++i)gtest[i]=head[i+8];
                gtest[4]='\0';
                if(string(gtest)!="SPHR")
                {
                    cerr << "File does not appear to be an RPI Spherical Grid"
                        <<endl
                        <<"Cannot continue"<<endl;
                    exit(-1);
                }
                for(i=0;i<4;++i)gtest[i]=head[i+12];
                cout << "Grid content (quantity) tag from header="<<gtest<<endl;
                for(i=0;i<4;++i)gtest[i]=head[i+4];
                gtest[4]='\0';
                string gridtype(gtest);
                delete [] head;
                if(gridtype=="CORS")
                    cout << "This appears to be a coarse grid"<<endl;
                else if(gridtype=="FINE")
                    cout << "This appears to be a fine grid"<<endl;
                else
                {
                    cout << "Unknown grid type ="<<gridtype
                        << "  Check file.  Exiting"<<endl;
                    exit(-1);
                }
                /* New format has the following attributes 
                   stored in the header section.  This was determined
                   from looking and c2f_sph and the spherical version
                   of the eikonal solver.*/
                const size_t base_offset(120);
                double fxs,fys,fzs;  // only relevant for fine grid - source 
                /* Using pointer arithmetic for this is perhaps a bit
                   obtuse, but should work */
                unsigned char *vhptr=reinterpret_cast<unsigned char *>(head);
                double *ptr;
                vhptr=vhptr+base_offset;
                ptr=reinterpret_cast<double*>(vhptr);
                fxs=*ptr;
                ++ptr;
                fys=*ptr;
                ++ptr;
                fzs=*ptr;
                ++ptr;
                double clat,clon,cz;
                clat=*ptr;
                ++ptr;
                clon=*ptr;
                ++ptr;
                cz=*ptr;
                ++ptr;
                double x0,y0,z0,dx,dy,dz;
                x0=*ptr;
                ++ptr;
                y0=*ptr;
                ++ptr;
                z0=*ptr;
                ++ptr;
                dx=*ptr;
                ++ptr;
                dy=*ptr;
                ++ptr;
                dz=*ptr;
                ++ptr;
                float *az;
                az=reinterpret_cast<float *>(ptr);
                // az will always be zero so we largely ignore it.
                ++az;
                int *iptr;
                iptr=reinterpret_cast<int *>(az);
                int nx,ny,nz;
                nx=*iptr;
                ++iptr;
                ny=*iptr;
                ++iptr;
                nz=*iptr;
                ++iptr;
                cout << "Header attributes extracted from file:"<<endl
                    << "fxs, fys, fzs="<<fxs<<", "<<fys<<", "<<fzs<<endl
                    << "clat, clon, cz="<<clat<<", "<<clon<<", "<<cz<<endl
                    << "x0, y0, z0="<<x0<<", "<<y0<<", "<<z0<<endl
                    << "dx(df), dy(dq), dz(h)="<<dx<<", "<<dy<<", "<<dz<<endl
                    << "nx, ny, nz="<<nx<<", "<<ny<<", "<<nz<<endl;
                /* I am pretty sure Steve stores lat,lon in degrees so
                   need to convert them to radians here */
                x0=rad(x0);
                y0=rad(y0);
                /* dx and dy are degrees.  Convert to radians */
                dx=rad(dx);
                dy=rad(dy);
                /* depths are km so I think that requires no change */
                
                /* New format uses these spacing arrays BUT only in
                the coarse grid mode  */
                int *igridx,*igridy,*igridz;
                if(gridtype=="CORS")
                {
                    igridx=new int[nx-1];
                    if(igridx==NULL) 
                        cerr << "Alloc failed for igridx of size="<<nx-1<<endl;
                    igridy=new int[ny-1];
                    if(igridy==NULL) 
                        cerr << "Alloc failed for igridy of size="<<ny-1<<endl;
                    igridz=new int[nz-1];
                    if(igridz==NULL) 
                        cerr << "Alloc failed for igridz of size="<<nz-1<<endl;
                    fread(igridx,sizeof(int),nx-1,fp);
                    fread(igridy,sizeof(int),ny-1,fp);
                    fread(igridz,sizeof(int),nz-1,fp);
                }
		int ntotal=nx*ny*nz;
		float *buffer=new float[ntotal];
		ntest = fread(buffer,sizeof(float),ntotal,fp);
		if(save_S)  fread(buffer,sizeof(float),ntotal,fp);
		if(ntest!=ntotal)
		{
			cerr << "Read error:  read "<<ntest
			   <<" numbers but expected" <<ntotal << endl;
			exit(-1);
		}
                /* These hold lon, lat, and depth of grid lines 
                   for both fine or coarse grids. Size different
                   in the two cases as can be seen how they 
                   are defined below */
                vector<double>xoffset,yoffset,zoffset;
		double x,y,z;
                /* The following all assume the origin (x0,y0,z0)
                   is in the upper left corner (top of NW corner) 
                   of the grid */
                if(gridtype=="CORS")
                {
                    cout << "xoffset 0="<<x0<<" deg->"<<deg(x0)<<endl;
                    xoffset.push_back(x0);
                    for(x=x0,i=0;i<(nx-1);++i)
                    {
                        x+=dx*static_cast<double>(igridx[i]);
                    cout << "xoffset "<<i<<"="<<x<<" deg->"<<deg(x)<<endl;
                        xoffset.push_back(x);
                    }
                    cout << "yoffset 0="<<y0<<" deg->"<<deg(y0)<<endl;
                    yoffset.push_back(y0);
                    for(y=y0,i=0;i<(ny-1);++i)
                    {
                        y-=dy*static_cast<double>(igridy[i]);
                    cout << "yoffset "<<i<<"="<<y<<" deg->"<<deg(y)<<endl;
                        yoffset.push_back(y);
                    }
                    cout << "zoffset 0="<<z0<<" deg->"<<deg(z0)<<endl;
                    zoffset.push_back(z0);
                    for(z=z0,i=0;i<(nz-1);++i)
                    {
                        z+=dz*static_cast<double>(igridz[i]);
                    cout << "zoffset "<<i<<"="<<z<<endl;
                        zoffset.push_back(z);
                    }
                }
                else
                {
                    // this is fine grid
                    for(i=0;i<nx;++i)
                    {
                        x=x0+dx*static_cast<double>(i);
                        xoffset.push_back(x);
                    }
                    for(i=0;i<ny;++i)
                    {
                        y=y0-dy*static_cast<double>(i);
                        yoffset.push_back(y);
                    }
                    for(i=0;i<nz;++i)
                    {
                        z=z0+dz*static_cast<double>(i);
                        zoffset.push_back(z);
                    }
                }
                // These will be wrong, but should be harmless
		int nx0,ny0;  
		nx0=nx/2;
		ny0=ny/2;
                /* Build the grid template. We will rebuild it below.
                   We use lat0, lon0, r0,and gridaz from Metadata
                   that sets the coordinate system.   Normally will
                   want that for a standard scene so best to build it
                   in to avoid a remap. */
		GCLgrid3d *g=new GCLgrid3d(nx,ny,nz,gridname,
			lat0,lon0,r0,gridaz,dx,dy,dz,nx0,ny0);
		// use this grid to create an empty scalar field that will hold
		// the velocity model we are going to input.
		GCLscalarfield3d vel(*g);
		// don't need g any longer
		delete g;
		
                /* Steve Roecker's grid format scans from the 
                   NW corner at the top.  First scan line is
                   then from NW to NE at the top of the grid.  
                   Next is south by one increment, etc.  Final
                   scan is top down.  This is almost completely
                   inverted from GCLgrid convention that wants the
                   origin in the lower left corner of the 
                   (distorted) box.  Handled here by counting 
                   z and y indices backward. */
		int offset=0;
                int jj;
		for(k=0,kk=vel.n3-1;k<vel.n3;++k,--kk)
			  for(j=0,jj=vel.n2-1;j<vel.n2;++j,--jj) 
				for(i=0;i<vel.n1;++i,++offset)

		{
			Geographic_point geo;
			Cartesian_point cart;
                        //WARNING:  unit skew likely fix
                        geo.lon=xoffset[i];
                        geo.lat=yoffset[j];
                        geo.r=r0_ellipse(geo.lat)-zoffset[k];
			vel.val[i][jj][kk]=buffer[offset];
			cart=vel.gtoc(geo);
			vel.x1[i][jj][kk]=cart.x1;
			vel.x2[i][jj][kk]=cart.x2;
			vel.x3[i][jj][kk]=cart.x3;
		}
		// save field to database
                string modelname;
		if(save_S)
			modelname=modelbasename+"_S";
		else
			modelname=modelbasename+"_P";
		if(usevpvs)
		{
			modelname=modelbasename+"_S";
			vel *= vpvs;
		}
		// Use the model name for the output data file name
                if(save_grid)
                {
                    cout << "Saving grid geometry (-savegrid) option."
                        <<endl
                        << "WARNING:  additional runs with this geometry should use default (no save)"
                        <<endl;
		    vel.dbsave(db,gclgdir,fielddir,modelname,modelname);
                }
                else
                {
                    cout<< "Warning:  grid geometry assumed already defined (default)"<<endl;
                    vel.dbsave(db,string(""),fielddir,modelname,modelname);
                }
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
		exit(-1);
	}
	catch (...)
	{
		cerr << "Something threw an exception that wasn't handled." << endl;
	}
}
