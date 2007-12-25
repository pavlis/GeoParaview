#include <stdio.h>
#include <string>
#include <list>
#include <istream>
#include "GCLgrid.h"
#include "dbpp.h"
#include "perf.h"
#include "GridScratchFileHandle.h"
void VectorField3DWeightedStack(GCLvectorfield3d& result,GridScratchHandle& handle,
		list<MemberGrid> mgYl)
{
	const string base_error("VectorField3DWeightedStack(Error):  ");
	list<MemberGrid>::iterator mptr;
	handle.rewind();
	result.zero();
	int nsize=result.n1*result.n2*result.n3*result.nv;
	auto_ptr<double> buffer(new double[nsize]);
	double totalweight,sumwt;
	int nread;
	for(sumwt=0.0,mptr=mgl.begin();mptr!=mptr.end();++mptr)
	{
		if(handle.at_eof())
			throw SeisppError(base_error
			  + string("EOF encountered prematurely on scratch file"));
		nread=handle.load_next(buffer);
		if(nread!=nsize)
			throw SeisppError(base_error
				+ string("short read from scratch file"));
		totalweight=(mptr->baseweight)*(mptr->reswt);
		daxpy(nsize,totalweight,buffer,1,result.val[0][0][0],1);
		sumwt+=totalweight;
	}
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
GCLgrid3d *decimate(GCLgrid3d& g,int dec1, int dec2, int dec3);
{
	int n1,n2,n3;
	n1=(g.n1)/dec1;
	n2=(g.n2)/dec2;
	n3=(g.n3)/dec3;
	/* Call an appropriate constructor here.  May need to reset
	some attributes of the grid object */
	GCLgrid3d *result;

	int i,j,k,ii,jj,kk;
	for(i=0,ii=0;i<g.n1 && ii<n1;i+=dec1,++ii)
	{
		for(j=0,jj=0;j<g.n2 && jj<n2;j+=dec2,++jj)
		{
			for(k=0,kk=0;k<g.n3 && kk<n3; j+=dec3,++kk)
			{
				result->x1[ii][jj][kk]=g.x1[i][j][k];
				result->x2[ii][jj][kk]=g.x2[i][j][k];
				result->x3[ii][jj][kk]=g.x3[i][j][k];
			}
		}
	}
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
	GridScatchFileHandle& gsfh,
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
within + and - of decimatio range.  At edges range is reduced by 
boundaries up to 1/2 of the full box size. 

Arguments:
	f - stack
	gsfh - scratch file holding migrated grids to be stacked.
	coh - output coherence grid
	dec1, dec2, and dec3 are decimation variables for x1,x2, and x3 
		coordinates (coh is assumed decimated by these amounts from
		f.  No checking is done.)
*/
void ComputeGridCoherence(GCLvectorfield3d& f, 
	GridScatchFileHandle& gsfh,
		list<MemberGrid> mgl,
			GCLscalarfield3d& coh,
				int dec1, 
					int dec2,
						int dec3)
{
	/* These may not be necessary, but cost for maintenance advantage is
	small */
	gsfh.rewind();
	coh.zero();
	/* indices for f */
	int i,j,k,l;
	/* ranges for computing each box (computed for each output grid point*/
	int imin,imax,jmin,jmax,kmin,kmax;
	/* indices or coh */
	int ic,jc,kc;
	int nsize=f.n1*f.n2*f.n3*f.nv;
	auto_ptr<double> buffer(new double[nsize]);
	/* This loop is similar to stacking loop, but we do not check dimension
	here assuming we can't get here if such an error is present.  Beware if
	you ever extract this function for another program.*/
	int nread;  /* Ignored return */
	double weight,sumwt;
	for(mptr=mgl.begin();mptr!=mptr.end();++mptr)
	{
		nread=handle.load_next(buffer);
		totalweight=(mptr->baseweight)*(mptr->reswt);
		daxpy(nsize,totalweight,buffer,1,result.val[0][0][0],1);
		sumwt+=totalweight;
	for(ic=0,i=0;ic<coh.n1;++ic)
	{
		if(i<dec1)
			imin=0;
		else
			imin=i-dec1;
	
}
	
void usage()
{
	cerr << "gridstacker db [-i infile -pf pffile -norobust] "<<endl;
	exit(-1);
}

int main(int argc, char **argv)
{
	istream in=cin;
	string pffile("gridstacker");
	if(argc<2) usage();
	string dbname(argc[1]);
	bool robust(true);

	int i;
	for(i=2;i<argc;++argv)
	{
		string sarg(argv[i]);
		if(sarg=="-i")
		{
			++i;
			in.open(argv[i], ios::in);
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
		else
			usage();
	}
	try {
		DatascopeHandle dbh(dbname,false);
		/* Build the list of grid definitions to form stack */
		list<MemberGrid> mgl;
		while(in.good())
		{
			string g,f;
			double wgt;
			in >> g;
			in >> f;
			in >> wgt;
			// Initialize residual weight portion to 1.0 at start
			mgl.push_back(MemberGrid(g,f,wgt,1.0));
		}
		Metadata control(pffile);
		string mastergridname=control.get_string("mastergridname");
		string masterfieldname=control.get_string("masterfieldname");
		string baseofn=control.get_string("base_output_field_name");
		GCLvectorfield3d mastergrid(dbh.db,mastergridname,masterfieldname);
		GridScatchFileHandle gsfh(mastergrid,mgl,dbh.db);
		gsfh.rewind();
		/* First compute and save a straight stack.  This works because
		mgl.reswt values are all set to 1.0 above. */
		GCLvectorfield3d result(mastergrid);  // clone grid
		result.zero();
		VectorField3DWeightedStack(result,gsfh,mgl);
		string avgfname=baseofn + "_avg";
		result.dbsave(dbh.db,mastergridname,avgfname);
		gsfh.rewind();
		/* Now compute and save a coherence attribute grid for straight stack.
		Grid is formed by decimation of master grid controlled by dec1,dec2,
		and dec3. */
		int dec1,dec2,dec3;
		dec1=control.get_int("coherence_x1_decimation_factor");
		dec2=control.get_int("coherence_x3_decimation_factor");
		dec3=control.get_int("coherence_x3_decimation_factor");
		string cohgname,cohbasefname;
		cohgname=control.get_string("coherence_grid_name");
		cohbasefname=control.get_string("coherence_field_base_name");
		GCLgrid3d *cohgrd;
		cohgrd=decimate(dynamic_cast<GCLgrid3d&>(mastergrid),dec1,dec2,dec3);
		GCLscalarfield3d coh(*cohgrd);
		coh.zero(); 
		delete cohgrd;  // no longer needed.  coh clones its geometry
		ComputeGridCoherence(result,gsfh,coh);
		// Assume this will throw an error if the grid name is already defined.
		// Probably should test first to make sure it does not exist already
		// as this will likely be a common error.  Ignored for now
		string sscohf=cohbasefname+"_avg";
		coh.dbsave(dbh.db,cohgname,sscohf);
		// Stop if robust is turned off
		if(!robust) exit(0);

		// Default is to use the average as the initial estimate
		// If desired use member of list with largest weight
		bool uselargest=control.get_bool("use_largest_weight_as_initial");
		if(uselargest) LoadLargestWeightMember(result,gsfh);
		
		

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
		
			
		
	
