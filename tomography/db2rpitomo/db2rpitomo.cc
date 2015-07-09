#include <stdio.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <string>
#include "seispp.h"
#include "dbpp.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "db2rpitomo db [-o outfile -full]"<<endl
	  << "Program prompts for inputs"<<endl
	  << "default outfile = rpitomo.dat"<<endl
	  << "default format is restest input.  -full flag writes full model"
	  <<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	if(argc<2) usage();
	string dbname(argv[1]);
	string outfile("rpitomo.dat");
	bool writefull(false);
	int i,j,k,kk;
	for(i=0;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-o")
		{
			++i;
			if(i>=argc) usage();
			outfile=string(argv[i]);
		}
		else if (sarg=="-full")
		{
			writefull=true;
		}
		else
			usage();
	}
	FILE *fp;
	fp=fopen(outfile.c_str(),"w");
	if(fp==NULL)
	{
		cerr << "Open failed on file="<<outfile<<endl;
		usage();
	}
	try {
		string gridname,fieldname;
		DatascopeHandle dbh(dbname,false);
		cout << "db2rpitomo:  will write output to file="<<outfile<<endl
			<< "Enter GCLgrid name for model to convert->";
		cin >> gridname;
		cout << endl<<"Enter GCLgrid field name or model to be converted->";
		cin >> fieldname;
		cout << "Attempting to load "<<gridname<<":"<<fieldname<<endl;
		GCLscalarfield3d f(dbh.db,gridname,fieldname);
		int ntotal;
		ntotal=(f.n1)*(f.n2)*(f.n3);
		float *buffer;
		if(writefull) buffer=new float[ntotal];
		int offset(0);
		double z,vel;
		for(k=0,kk=f.n3-1;k<f.n3;++k,--kk)
                  for(j=0;j<f.n2;++j)
                    for(i=0;i<f.n1;++i,++offset)
		    {
			if(fabs(vel)>FLT_EPSILON)
			{
				vel=f.val[i][j][kk];
				/* Intentionally do not use flattening
				in this mode */
				if(!writefull)
				  fprintf(fp,"%d %lf\n",offset,vel);
				/* apply it now so it gets set in output*/
			    	z=f.depth(i,j,kk);
				vel=flatvel(vel,z);
				f.val[i][j][kk]=vel;
			}
			if(writefull)
		    		buffer[offset]=f.val[i][j][kk];
		    }
		if(writefull) fwrite(buffer,sizeof(float),ntotal,fp);
		delete [] buffer;
		fclose(fp);

	} 
	catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (int ierr)
	{
		cerr << "GCLgrid error:  Something threw int excpetion = "<<ierr<<endl;
	}
}


