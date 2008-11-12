#include <stdio.h>
#include <float.h>
#include <string>
#include <list>
#include <fstream>
#include "gclgrid.h"
#include "dbpp.h"
#include "perf.h"
#include "seispp.h"
#include "GridScratchFileHandle.h"
#include "GridStackPenaltyFunction.h"
void VectorField3DWeightedStack(GCLvectorfield3d& result,GridScratchFileHandle& handle,
		list<MemberGrid>& mgl)
{
	const string base_error("VectorField3DWeightedStack(Error):  ");
	list<MemberGrid>::iterator mptr;
	handle.rewind();
	result.zero();
	int nsize=result.n1*result.n2*result.n3*result.nv;
	double *buffer=new double[nsize];
	double totalweight,sumwt;
	int nread;
	for(sumwt=0.0,mptr=mgl.begin();mptr!=mgl.end();++mptr)
	{
		if(handle.at_eof())
		{
			delete [] buffer;
			throw SeisppError(base_error
			  + string("EOF encountered prematurely on scratch file"));
		}
		nread=handle.load_next(buffer);
		if(nread!=nsize)
		{
			delete [] buffer;
			throw SeisppError(base_error
				+ string("short read from scratch file"));
		}
		totalweight=(mptr->baseweight)*(mptr->reswt);
		daxpy(nsize,totalweight,buffer,1,result.val[0][0][0],1);
		sumwt+=totalweight;
	}
	delete [] buffer;
	if(!handle.at_eof())
		throw SeisppError(base_error
		  + string("Scratch file mismatch.")
		  + string("  Not at EOF after working through member list"));
	/* Always normalize by sumwt*/
	dscal(nsize,1.0/sumwt,result.val[0][0][0],1);
}
/*! \brief Decimate a GCLgrid3d.

We sometimes want to decimate a grid. This procedure does this for
a 3D grid with variable decimation for each generalized coordinate 
axis.

\param g parent grid that is to be decimated
\param dec1 decimation factor for x1
\param dec2 decimation factor for x2
\param dec3 decimamtion factor for x3

\return pointer to decimated grid object 
*/
GCLgrid3d *decimate(GCLgrid3d& g,int dec1, int dec2, int dec3)
{
	int n1,n2,n3;
	n1=(g.n1)/dec1;
	n2=(g.n2)/dec2;
	n3=(g.n3)/dec3;
	if( (n1<2) || (n2<2) || (n3<2) )
	  throw SeisppError(string("grid decimator:  decimation factor larger than grid dimension"));
	/* Call an appropriate constructor here.  May need to reset
	some attributes of the grid object */
	GCLgrid3d *result=new GCLgrid3d(n1,n2,n3);;

	int i,j,k,ii,jj,kk;
	for(i=0,ii=0;i<g.n1 && ii<n1;i+=dec1,++ii)
	{
		for(j=0,jj=0;j<g.n2 && jj<n2;j+=dec2,++jj)
		{
			for(k=0,kk=0;k<g.n3 && kk<n3; k+=dec3,++kk)
			{
				result->x1[ii][jj][kk]=g.x1[i][j][k];
				result->x2[ii][jj][kk]=g.x2[i][j][k];
				result->x3[ii][jj][kk]=g.x3[i][j][k];
			}
		}
	}
	/* We clone other attribures of the parent grid.*/
	result->lat0=g.lat0;
	result->lon0=g.lon0;
	result->r0=g.r0;
	result->azimuth_y=g.azimuth_y;
	result->x1low=g.x1low;
	result->x1high=g.x1high;
	result->x2low=g.x2low;
	result->x2high=g.x2high;
	result->x3low=g.x3low;
	result->x3high=g.x3high;
	result->dx1_nom=g.dx1_nom*(static_cast<double>(dec1));
	result->dx2_nom=g.dx2_nom*(static_cast<double>(dec2));
	result->dx3_nom=g.dx3_nom*(static_cast<double>(dec3));
	result->i0=g.i0/dec1;
	result->j0=g.j0/dec2;
	result->k0=g.k0/dec3;
	for(i=0;i<3;++i)
		for(j=0;j<3;++j) 
		    result->gtoc_rmatrix[i][j]=g.gtoc_rmatrix[i][j];
	for(i=0;i<3;++i)result->translation_vector[i]
			 = g.translation_vector[i];

	return(result);
}
/* This uses a very crude and slow algorithm to find the 
grid with the highest weight.  It is about the only way we
can do this, however, because the scratch file as currently implemented
is only accessed sequentially.  Could add a directory of foff values
if this proves too slow.
 f will contain the result on return
 gsfh  is the scratch file handle
 mgl is the list that contains attribute info on the members in the scratch
	file.

\return size of largest weight of member chosen
*/
double LoadLargestWeightMember(GCLvectorfield3d& f,
	GridScratchFileHandle& gsfh,
		list<MemberGrid>& mgl)
{
	int nval=(f.n1)*(f.n2)*(f.n3)*(f.nv);
	double *buffer=new double[nval];
	// better safe than sorry
	gsfh.rewind();
	list<MemberGrid>::iterator mglptr;
	double maxwgt=0.0;
	for(mglptr=mgl.begin();mglptr!=mgl.end();++mglptr)
	{
		if(gsfh.load_next(buffer) != nval)
		{
			delete [] buffer;
			throw SeisppError(string("LoadLargestWeight:  ")
			  + string("error reading scratch file"));
		}
		if(mglptr->baseweight > maxwgt)
		{
			maxwgt=mglptr->baseweight;
			// This depends on GCLgrid field values being
			// stored on a contiguous block.  WARNING this
			// is not a requirement of that interface.
			memcpy(static_cast<void *>(f.val[0][0][0]),
				static_cast<const void *>(buffer),
				nval*sizeof(double));
		}
	}
	return(maxwgt);
}
/* Computes coherence grid smoothed by decimation variables.  That is,
coherence at each point in the output grid is produced from field values
within + and - of decimatio range.  At edges the boxes are still 
dec1 X dec2 X dec3, but the averaging boxes are not allowed to pass beyond
the boundaries of the grid.  The points on the outer edge of the result
are thus somewhat distorted, but for normal uses this should not be a
large issue.   Note that boxes with no data will have a value of zero 
on output

Arguments:
	f - stack
	gsfh - scratch file holding migrated grids to be stacked.
	coh - output coherence grid
	dec1, dec2, and dec3 are decimation variables for x1,x2, and x3 
		coordinates (coh is assumed decimated by these amounts from
		f.  No checking is done.)
	component - vector component to use for computing coherence.  
		Normally should be 0,1, or 2.  If negative or greater than
		2 this routine SILENTLY switches to a three-component mode
		in which residuals from all three components are used to 
		compute coherence. 
*/
void ComputeGridCoherence(GCLvectorfield3d& f, 
	GridScratchFileHandle& gsfh,
		list<MemberGrid>& mgl,
			GCLscalarfield3d& coh,
				int dec1, 
					int dec2,
						int dec3,
							int component)
{
	/* decide if this is to be a 3c mode run */
	bool mode3c;
	if( (component<0) || (component>2) ) 
		mode3c=true;
	else
		mode3c=false;
	/* This is a cutoff constant on weight component of each grid point.
	Using this is safer than a test against zero to skip points with no data */
	const double minweight(0.00001);
	/* These may not be necessary, but cost for maintenance advantage is
	small */
	gsfh.rewind();
	coh.zero();
	/* indices for f */
	int i,j,k,l;
	int nsize=f.n1*f.n2*f.n3*f.nv;
	int npoints=f.n1*f.n2*f.n3;
	double ***sumsqr, ***sumsqd;
	sumsqr=create_3dgrid_contiguous(f.n1,f.n2,f.n3);
	sumsqd=create_3dgrid_contiguous(f.n1,f.n2,f.n3);
	/* This is a huge memory model to do this calculation, but pwmig has
	more problems so this shouldn't be an issue for any pwmig output */
	double ****buffer;
	buffer=create_4dgrid_contiguous(f.n1,f.n2,f.n3,f.nv);
	/* This loop is similar to stacking loop, but we do not check dimension
	here assuming we can't get here if such an error is present.  Beware if
	you ever extract this function for another program.   On the other hand
	the purpose of this loop is quite different.  Here the total purpose
	is to accumulate sum of squared residuals an sum of squared data values.*/
	int nread;  /* Ignored return */
	/* Loop over all grids in the scratch file accumulating weighted sum of squares of
	residuals and weighted sum of squared data values.  These are used later to 
	compute coherence in larger boxes.  This is an efficient algorithm as it 
	allows passing through the scratch file only once. */
	double weight,sumwt,sumr3sq,sumd3sq,val;
//DEBUG
vector<double> ddbg;
	for(i=0;i<f.n1;++i)
	 for(j=0;j<f.n2;++j)
	  for(k=0;k<f.n3;++k)
	  {
	    sumsqr[i][j][k]=0.0;
	    sumsqd[i][j][k]=0.0;
	  }
	list<MemberGrid>::iterator mptr; 
	for(mptr=mgl.begin();mptr!=mgl.end();++mptr)
	{
		/*Tricky use of a pointer here to get address of start of buffer */
		nread=gsfh.load_next(buffer[0][0][0]);
		double totalweight=(mptr->baseweight)*(mptr->reswt);
		sumwt+=totalweight;
		if(mode3c)
		{
		    for(i=0;i<f.n1;++i)
		    {
			for(j=0;j<f.n2;++j)
			{
				for(k=0;k<f.n3;++k)
				{
				    if(f.val[i][j][k][4]>minweight)
				    {
					for(l=0,sumd3sq=0.0,sumr3sq=0.0;l<3;++l)
					{
						val=buffer[i][j][k][l];
						sumd3sq=val*val;
						/* subtract stack value */
						val-=f.val[i][j][k][l];
						sumr3sq=val*val;
					}
					sumsqr[i][j][k]+=sumr3sq*totalweight*totalweight;
					sumsqd[i][j][k]+=sumd3sq*totalweight*totalweight;
				    }
				}
			}
		    }
		}
		else
		{
		    for(i=0;i<f.n1;++i)
		    {
			for(j=0;j<f.n2;++j)
			{
				for(k=0;k<f.n3;++k)
				{
				    if(f.val[i][j][k][4]>minweight)
				    {
					val=buffer[i][j][k][component];
					sumd3sq=val*val;
//DEBUG
ddbg.push_back(val);
					/* subtract stack value */
					val-=f.val[i][j][k][component];
					sumr3sq=val*val;
					sumsqr[i][j][k]+=sumr3sq*totalweight*totalweight;
					sumsqd[i][j][k]+=sumd3sq*totalweight*totalweight;
				    }
				}
			}
		    }
		}
	}
	/* We can release this now as what we need is now in the sumsqr an sumsqd grids */
	free_4dgrid_contiguous(buffer,f.n1,f.n2,f.n3);
	/* Now we loop over coh grid computing coherence of values averaged over boxes*/
	/* ranges for computing each box (computed for each output grid point*/
	int imin,imax,jmin,jmax,kmin,kmax;
	/* indices or coh */
	int ic,jc,kc;
	double x,y,z;
	coh.zero();
	int index[3];
	int dec1o2=dec1/2;
	int dec2o2=dec2/2;
	int dec3o2=dec3/2;
	for(ic=0;ic<coh.n1;++ic)
	{
		for(jc=0;jc<coh.n2;++jc)
		{
			for(kc=0;kc<coh.n3;++kc)
			{
				/* because the coh grid is obtained by decimation
				of the f grid we can assume the cartesian coordinates
				are congruent.  WARNING:  potential problem if this
				assumption gets violated as code evolves. */;
				x=coh.x1[ic][jc][kc];
				y=coh.x2[ic][jc][kc];
				z=coh.x3[ic][jc][kc];
				/* This method returns 0 on success and an error
				code if failure of any kind happens. We do nothing
				for failure case as that will leave the associated
				cell zero */
				if(!f.lookup(x,y,z))
				{
					f.get_index(index);
					/* safely get top left corner of averaging box */
					if(index[0]-dec1o2<0)
						imin=0;
					else
						imin=index[0]-dec1o2;
					if(index[1]-dec2o2<0)
						jmin=0;
					else
						jmin=index[1]-dec2o2;
					if(index[2]-dec3o2<0)
						kmin=0;
					else
						kmin=index[2]-dec3o2;
					/* now get opposite corner, resettng min values
					when hitting the upper wall.  When resetting we
					assume decimation factors are smaller than the 
					grid size in all directions. */
					imax=imin+dec1;
					if(imax>=f.n1)
					{
						imax = f.n1 - 1;
						imin = imax - dec1 + 1;
					}
					jmax=jmin+dec2;
					if(jmax>=f.n2)
					{
						jmax = f.n2 - 1;
						jmin = jmax - dec2 + 1;
					}
					kmax=kmin+dec3;
					if(kmax>=f.n3)
					{
						kmax = f.n3 - 1;
						kmin = kmax - dec3 + 1;
					}
				}
				for(i=imin;i<=imax;++i)
				  for(j=jmin;j<=jmax;++j)
				    for(k=kmin,sumr3sq=0.0,sumd3sq=0.0;k<=kmax;++k)
				    {
					sumr3sq+=sumsqr[i][j][k];
					sumd3sq+=sumsqd[i][j][k];
				    }
				if(fabs(sumd3sq)<DBL_EPSILON) 
					val=0.0;
				else
					val=sqrt(sumr3sq/sumd3sq);
				if(val>1.0)
					coh.val[ic][jc][kc]=0.0;
				else
					coh.val[ic][jc][kc]=1.0-val;
			}
		}
	}
	free_3dgrid_contiguous(sumsqr,f.n1,f.n2);
	free_3dgrid_contiguous(sumsqd,f.n1,f.n2);
}
void compute_reswt(GCLvectorfield3d& stack,GridScratchFileHandle& handle,
                list<MemberGrid>& mgl,GridStackPenaltyFunction& penalty_function)
{
	const string base_error("compute_reswt function:  ");
	list<MemberGrid>::iterator mptr;
	handle.rewind();
	int nsize=stack.n1*stack.n2*stack.n3*stack.nv;
	double *buffer=new double[nsize];
	int nread;
	if(SEISPP_verbose) cout << "fieldname baseweight residual_weight"<<endl;
	for(mptr=mgl.begin();mptr!=mgl.end();++mptr)
	{
		if(handle.at_eof())
		{
			delete [] buffer;
			throw SeisppError(base_error
			  + string("EOF encountered prematurely on scratch file"));
		}
		nread=handle.load_next(buffer);
		if(nread!=nsize)
		{
			delete [] buffer;
			throw SeisppError(base_error
				+ string("short read from scratch file"));
		}
		mptr->reswt=penalty_function.weight(nsize,buffer,stack.val[0][0][0]);
		if(SEISPP_verbose) cout << mptr->fieldname
				<< " " << mptr->baseweight
				<< " " << mptr->reswt <<endl;
	}
}
	
void usage()
{
	cerr << "gridstacker db [-i infile -pf pffile -norobust -V] "<<endl;
	exit(-1);
}
void inreadabort()
{
	cerr << "Error reading input file of grids to stack.  Check format"
		<<endl;
	exit(-1);
}

bool reswt_change_test(list<MemberGrid>& mgltest, list<MemberGrid>& mgl0)
{
	const double nrmcut(0.05);
	list<MemberGrid>::iterator mgptr0,mgptr;
	if(mgltest.size()!=mgl0.size())
	{
		cerr << "Size mismatch in MemberGrid lists. "
			<< "Something bad has happened.  exiting."<<endl;
		exit(-1);
	}
	double sumsqr0(0.0),sumsqr(0.0);
	for(mgptr=mgltest.begin(),mgptr0=mgl0.begin();mgptr!=mgltest.end();
		++mgptr,++mgptr0)
	{
		sumsqr0 += (mgptr0->reswt)*(mgptr0->reswt);
		sumsqr += (mgptr->reswt)*(mgptr->reswt);
	}
	sumsqr0=sqrt(sumsqr0);
	sumsqr=sqrt(sumsqr);
	double L2nrmchange=fabs(sumsqr-sumsqr0)/sumsqr0;
	if(L2nrmchange>nrmcut) return true;
	return false;
	
}
/* Used to test if a grid file name key is already in the db. 
This is a very very slow algorithm and should not be used often.
Appropriate here as it is called once and only once in a give run. */
bool GridNameDefined(DatascopeHandle dbh, string gridname, int dimensions)
{
	char line[256];
	stringstream ss(line);
	dbh.lookup("gclgdisk");
	ss << "gridname =~/" << gridname << "/ && dimensions == "<<dimensions;
	dbh.subset(ss.str());
	if(dbh.number_tuples() > 0)
		return true;
	else
		return false;
}
bool SEISPP::SEISPP_verbose(false);

int main(int argc, char **argv)
{
	//istream in=cin;
	ifstream in;
	string pffile("gridstacker");
	if(argc<2) usage();
	string dbname(argv[1]);
	bool robust(true);

	int i;
	for(i=2;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-i")
		{
			++i;
			string infile(argv[i]);
			in.open(infile.c_str(), ios::in);
			if(in.fail())
			{
				cerr << "Open failure on -i file="
					<<argv[i]<<endl;
				usage();
			}
		}
		else if(sarg=="-pf")
		{
			++i;
			pffile=string(argv[i]);
		}
		else if(sarg=="-norobust")
			robust=false;
		else if(sarg=="-V")
			SEISPP_verbose=true;
		else
			usage();
	}
	try {
		DatascopeHandle dbh(dbname,false);
		/* Build the list of grid definitions to form stack */
		list<MemberGrid> mgl;
		/* Old loop top, retain in case simpler version fails 
		while(in.good())
		{
			string g,f;
			double wgt;
			in >> g;
			if(!in.good()) break;
			in >> f;
			if(!in.good()) 
				inreadabort();
			in >> wgt;
			if(!in.good()) 
				inreadabort();
		*/
		char *inputline=new char[512];
		while(in.getline(inputline,512))
		{
			stringstream ss(inputline);
			string g,f;
			const double wgt(1.0);
			ss >> g;
			ss >> f;
			// Initialize residual weight portion to 1.0 at start
			mgl.push_back(MemberGrid(g,f,wgt,1.0));
		}
		delete [] inputline;
		/* Extract control parameters */
		Pf *pf;
		if(pfread(const_cast<char *>(pffile.c_str()),&pf))
		{
			cerr << "pfread failure for pffile="<<pffile<<endl;
			exit(-1);
		}
		Metadata control(pf);
		string mastergridname=control.get_string("mastergridname");
		string masterfieldname=control.get_string("masterfieldname");
		string baseofn=control.get_string("base_output_field_name");
		string outdir=control.get_string("output_data_directory");
		int dec1,dec2,dec3;
		dec1=control.get_int("coherence_x1_decimation_factor");
		dec2=control.get_int("coherence_x3_decimation_factor");
		dec3=control.get_int("coherence_x3_decimation_factor");
		int coherence_component;
		coherence_component=control.get_int("coherence_component");
		string cohgname,cohbasefname;
		cohgname=control.get_string("coherence_grid_name");
		string cohgdir;
		if(GridNameDefined(dbh,cohgname,3) )
		{
			cerr << "gridstacker (WARNING):  grid for stack coherence "
				<< "named " <<cohgname<<" already exists in the database"<<endl
				<< "This program assumes this grid geometry is consistent with this run"
				<<endl;
			cohgdir=string("");
		}
		else
		{
			cohgdir=outdir;
		}
		cohbasefname=control.get_string("coherence_field_base_name");
		// Default is to use the average as the initial estimate
		// If this is set true wil use member of list with largest weight as initial
		bool uselargest=control.get_bool("use_largest_weight_as_initial");
		/* This object defines the type of robust stack to use.  We construct
		it even if robust stack is not requested as at present the construction
		is trivial and doesn't waste any memory. */
		GridStackPenaltyFunction penalty_function(control);

		/* End section extracting control parameters */

		GCLvectorfield3d mastergrid(dbh.db,mastergridname,masterfieldname);
		GridScratchFileHandle gsfh(mastergrid,mgl,dbh.db);
		gsfh.rewind();
		/* First compute and save a straight stack.  This works because
		mgl.reswt values are all set to 1.0 above. */
		GCLvectorfield3d result(mastergrid);  // clone grid
		result.zero();
		VectorField3DWeightedStack(result,gsfh,mgl);
		string avgfname=baseofn + "_avg";
		result.dbsave(dbh.db,string(""),outdir,avgfname,avgfname);
		gsfh.rewind();
		/* Now compute and save a coherence attribute grid for straight stack.
		Grid is formed by decimation of master grid controlled by dec1,dec2,
		and dec3. */
		GCLgrid3d *cohgrd;
		cohgrd=decimate(dynamic_cast<GCLgrid3d&>(mastergrid),dec1,dec2,dec3);
		GCLscalarfield3d coh(*cohgrd);
		coh.zero(); 
		delete cohgrd;  // no longer needed.  coh clones its geometry
		/* Need to set the grid name */
		coh.name=cohgname;
		ComputeGridCoherence(result,gsfh,mgl,coh,dec1,dec2,dec3,coherence_component);
		// Assume this will throw an error if the grid name is already defined.
		// Probably should test first to make sure it does not exist already
		// as this will likely be a common error.  Ignored for now
		string sscohf=cohbasefname+"_avg";
		coh.dbsave(dbh.db,cohgdir,outdir,sscohf,sscohf);
		/* At this point no matter what we set cohgdir as a null string.  This is
		a weird feature of the gclgrid library that turns off saving grid file */
		cohgdir=string("");
		// Stop if robust is turned off
		if(!robust) exit(0);
		/* this replaces straight stack with largest weight member as initial
		estimate for robust stacking method */
		if(uselargest) LoadLargestWeightMember(result,gsfh,mgl);
		int robust_loop_count=0;
		const int MAXIT(20);
		
		list<MemberGrid> mgl_previous;
		mgl_previous=mgl;
		bool keep_looping;
		do {
			gsfh.rewind();
			compute_reswt(result,gsfh,mgl,penalty_function);
			gsfh.rewind();
			VectorField3DWeightedStack(result,gsfh,mgl);
			keep_looping=reswt_change_test(mgl,mgl_previous);
			mgl_previous=mgl;
			++robust_loop_count;
		} while ( (robust_loop_count<MAXIT) && keep_looping);
		string robustfname=baseofn+"_ravg";
		result.dbsave(dbh.db,string(""),outdir,robustfname,robustfname);
		gsfh.rewind();
		ComputeGridCoherence(result,gsfh,mgl,coh,dec1,dec2,dec3,coherence_component);
		string robustcohf=cohbasefname+"_ravg";
		coh.dbsave(dbh.db,cohgdir,outdir,robustcohf,robustcohf);

	} catch (int ierr)
	{
		cerr << "GCLgrid library error number "<<ierr<< " was thrown"<<endl;
		exit(-1);
	}
	catch (SeisppError serr)
	{
		serr.log_error();
		exit(-1);
	}
}
		
