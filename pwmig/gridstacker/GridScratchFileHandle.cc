#include <stdio.h>
#include <list>
#include "seispp.h"
#include "GridScratchFileHandle.h"

MemberGrid::MemberGrid(const MemberGrid& parent);
{
	gridname=parent.gridname;
	fieldname=parent.fieldname;
	baseweight=parent.baseweight;
	reswt=parent.reswt;
}

MemberGrid& MemberGrid::operator= (const MemberGrid& parent);
{
	if(this!=&parent)
	{
		gridname=parent.gridname;
		fieldname=parent.fieldname;
		baseweight=parent.baseweight;
		reswt=parent.reswt;
	}
	return(*this);
}

/* More rigorous test for equality than operator== required for this code */
bool compare_equal(GCLvectorfield3d& g1, GCLvectorfield3d& g2)
{
	if(g1!=g1) return(false);
	if(g1.n1!=g2.n2) return(false);
	if(g1.n2!=g2.n2) return(false);
	if(g1.n3!=g2.n3) return(false);
	if(g1.nv!=g2.nv) return(false);
}

GridScratchFileHandle::GridScratchFileHandle(GCLvectorfield3d& mastergrid,
		list<MemberGrid> mg, Dbptr db)
{
	const string base_error("GridScratchFileHandle constructor:  ");
	fp=tmpfile();
	if(fp==NULL) throw SeisppError(base_error
			+ string("tmpfile() failed"));
	n1=mastergrid.n1;
	n2=mastergrid.n2;
	n3=mastergrid.n3;
	nv=mastergrid.nv;
	nval=n1*n2*n3*nv;
	current_member=0;
	nmembers=0;  // Could use size on list, but better to count grids loaded
	list<MemberGrid>::iterator mptr;
	for(mptr=mgl.begin();mptr!=mptr.end();++mptr)
	{
	    try{
		GCLvectorfield3d current(db,mptr->gridname,mptr->fieldname,nv);
		if(compare_equal(mastergrid,current))
		{
			if(current.nv==mastergrid.nv)
				fwrite(current.val[0][0][0],sizeof(double),nval,fp);
			else
			{
				cerr << base_error
					<< "vector field size mismatch for fieldname="
					<< mptr->fieldname
					<< endl
					<< "Required nv of master="
					<<mastergrid.nv
					<< "this member's nv="
					<< current.nv
					<< endl
					<< "Skipping data for this grid"
					<< endl;
			}
		}
		else
		{
			GCLvectorfield3d scratch(mastergrid);
			scratch.zero();
			scratch += current;
			fwrite(scratch.val[0][0][0],sizeof(double),nval,fp);
		}
		++nmembers;
		++current_member;
	    } catch (int ierr) {
		cerr << base_error
			<< "GCLvectorfield3d constructor failed for grid="
			<<mptr->gridname
			<< " fieldname="
			<< mptr->fieldname
			<<endl
			<<"These data will be dropped from the stack"
			<<endl
	    }    
	}
		
}
GridScatchFileHandle::~GridScratchFileHandle();
{
	if(fp!=NULL) fclose(fp);
}

void GridScratchFileHandle::rewind()
{
	rewind(fp);
	current_member=0;
}
// The format of this is unquestionably wrong.  NO references to consult
void GridScratchFileHandle::operator ++()
{
	double *buffer=new double[nval];
	int nread=fread(buffer,sizeof(double),nval,fp);
	delete [] buffer;
	if(nread!=nval) throw SeisppError(
		string("GridScratchFileHandle::operator++:  fread error on scratch file"));
	++current_member;

}
int GridScratchFileHandle::load_next(double *buffer)
{
	int nread=fread(buffer,sizeof(double),nval,fp);
	if(nread!=nval) throw SeisppError(
		string("GridScratchFileHandle::load_next:  fread error"));
	++current_member;
}
bool at_eof()
{
	if(current_member==nmember)
		return true;
	else
		return false;
}

