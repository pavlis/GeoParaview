#include <cstdio>
#include <cstring>
#include "gclgrid.h"
// Constructor creating a GCLgrid object by reading data from
// an Antelope database.  Uses an extension table used to index
// GCLgrid objects.  This is the 3d version.  The 2d Version 
// is below.
//
//Arguments:   
//	db - database to access to look up this grid file under
//	gridnmae - name of requested grid object  (key of the table)
//
//Throws a nonzero integer exception when creation of the GCLgrid 
//object fails.  Every thrown error here means the process failed
//completely so this condition must be caught and handled.  The
//details are assumed to be told to the user through messages
//written through the antelope elog mechanism.
GCLgrid3d::GCLgrid3d(Dbptr db, char *gridname)
{
	Dbptr dbgrd;
	char sstring[40];
	char datatype[4];
	char dir[65], dfile[36];
	char filename[512];  
	int foff;
	char *base_message="Cannot create GCL3Dgrid object: ";
	char *read_error="fread failed";
	int nrec;
	char cdef[2],gdef[2];
	int retcode;
	int gridsize;  
	FILE *fp;

	dbgrd = dblookup(db,0,"gclgdisk",0,0); 
	if(dbgrd.record == dbINVALID) 
	{
		elog_notify(0,"%s gclgdisk table not defined in schema definition\n",base_message);
		throw 1;
	}
	sprintf(sstring,"gridname =~ /%s/ && dimensions == 3",gridname);
	dbgrd = dbsubset(dbgrd,sstring,0);
	dbquery(dbgrd,dbRECORD_COUNT,&nrec);
	if(nrec <= 0) 
	{
		elog_notify(0,"%s grid with name %s not found in database\n",
			base_message,gridname);
		throw 1;
	}
	dbgrd.record = 0;
	strcpy(name,gridname);
	/* We intentionally ignore the dimension field because the subset
	should have assured we match the right record */
	if(dbgetv(dbgrd,0,
		"lat",&(lat0),
		"lon",&(lon0),
		"radius",&(r0),
		"azimuth_y",&(azimuth_y),
		"dx1nom",&(dx1_nom),
		"dx2nom",&(dx2_nom),
		"dx3nom",&(dx3_nom),
		"n1",&(n1),
		"n2",&(n2),
		"n3",&(n3),
		"xlow",&(xlow),
		"xhigh",&(xhigh),
		"ylow",&(ylow),
		"yhigh",&(yhigh),
		"zlow",&(zlow),
		"zhigh",&(zhigh),
		"i0",&(i0),
		"j0",&(j0),
		"k0",&(k0),
		"cdefined",cdef,
		"geodefined",gdef,
		"datatype",datatype,
		"dir",dir,
		"dfile",dfile,
		"foff",&foff,
		0) == dbINVALID)
	{
		elog_notify(0,"%s dbgetv error reading gclgdisk table\n",
			base_message);
		throw 1;
	}
	/* These parameters are stored in the database in degrees but
	are converted to radians for internal use */
	lat0 = rad(lat0);
	lon0 = rad(lon0);
	azimuth_y = rad(azimuth_y);
	if(!strcmp(cdef,"n") || !strcmp(gdef,"n") )
	{
		elog_notify(0,"%s Cartesian and Geographical mapping (cdefined and geodefined attributes) must both be defined for input\n",
			base_message);
		throw 1;
	}
	else
	{
		cartesian_defined=1;
		geographic_defined=1;
	}
	if(strcmp(datatype,"t8"))
	{
		elog_notify(0,"%s data type %s not allowed.  Currently only support t8\n",
			base_message, datatype);
		throw 1;
	}
	/* Get the file name to read the gclgrid data from.*/
	if(dbextfile(dbgrd,"gclgdisk",filename) <=0)
	{
		elog_notify(0,"%s Cannot find grid file named %s in directory %s\n",
			base_message,dfile,dir);
		throw 1;
	}
	fp = fopen(filename,"r");
	if(fp == NULL)
	{
		elog_notify(0,"%s file %s cannot be openned for read\n",
			base_message,filename);
		throw 1;
	}
	fseek(fp,foff,SEEK_SET);
	/* We alloc all memory first before reading so we can call a 
	free routine in case any reads fail */
	x1 = create_3dgrid_contiguous(n1,n2,n3);
	x2 = create_3dgrid_contiguous(n1,n2,n3); 
	x3 = create_3dgrid_contiguous(n1,n2,n3); 
	lat = create_3dgrid_contiguous(n1,n2,n3);
	lon = create_3dgrid_contiguous(n1,n2,n3); 
	r = create_3dgrid_contiguous(n1,n2,n3); 

	//read data from file.  Assume the destructor will be called
	//when throw an exception
	gridsize = (n1)*(n2)*(n3);
	if(fread(x1[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading x1\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(x2[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading x2\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(x3[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading x3\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(lat[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading latitude\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(lon[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading longitude\n",
			base_message,read_error);
		throw 1;
	}

	if(fread(r[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading radius\n",
			base_message,read_error);
		throw 1;
	}
}
/* close companion to the above for 2d grid */
GCLgrid::GCLgrid(Dbptr db,char *gridname)
{
	Dbptr dbgrd;
	char sstring[40];
	char datatype[4];
	char dir[65], dfile[36];
	char filename[512];  
	int foff;
	char *base_message="Cannot create GCL2Dgrid object: ";
	char *read_error="fread failed";
	int nrec;
	char cdef[2],gdef[2];
	int retcode;
	int gridsize;  
	FILE *fp;

	dbgrd = dblookup(db,0,"gclgdisk",0,0); 
	if(dbgrd.record == dbINVALID) 
	{
		elog_notify(0,"%s gclgdisk table not defined in schema definition\n",base_message);
		throw 1;
	}
	sprintf(sstring,"gridname =~ /%s/ && dimensions == 2",gridname);
	dbgrd = dbsubset(dbgrd,sstring,0);
	dbquery(dbgrd,dbRECORD_COUNT,&nrec);
	if(nrec <= 0) 
	{
		elog_notify(0,"%s grid with name %s not found in database\n",
			base_message,gridname);
		throw 1;
	}
	dbgrd.record = 0;
	strcpy(name,gridname);
	/* We intentionally ignore the dimension field because the subset
	should have assured we match the right record */
	if(dbgetv(dbgrd,0,
		"lat",&(lat0),
		"lon",&(lon0),
		"radius",&(r0),
		"azimuth_y",&(azimuth_y),
		"dx1nom",&(dx1_nom),
		"dx2nom",&(dx2_nom),
		"n1",&(n1),
		"n2",&(n2),
		"xlow",&(xlow),
		"xhigh",&(xhigh),
		"ylow",&(ylow),
		"yhigh",&(yhigh),
		"zlow",&(zlow),
		"zhigh",&(zhigh),
		"i0",&(i0),
		"j0",&(j0),
		"cdefined",cdef,
		"geodefined",gdef,
		"datatype",datatype,
		"dir",dir,
		"dfile",dfile,
		"foff",&foff,
		0) == dbINVALID)
	{
		elog_notify(0,"%s dbgetv error reading gclgdisk table\n",
			base_message);
		throw 1;
	}
	/* These parameters are stored in the database in degrees but
	are converted to radians for internal use */
	lat0 = rad(lat0);
	lon0 = rad(lon0);
	azimuth_y = rad(azimuth_y);

	if(!strcmp(cdef,"n") || !strcmp(gdef,"n") )
	{
		elog_notify(0,"%s Cartesian and Geographical mapping (cdefined and geodefined attributes) must both be defined for input\n",
			base_message);
		throw 1;
	}
	else
	{
		cartesian_defined=1;
		geographic_defined=1;
	}
	if(strcmp(datatype,"t8"))
	{
		elog_notify(0,"%s data type %s not allowed.  Currently only support t8\n",
			base_message, datatype);
		throw 1;
	}
	/* Get the file name to read the gclgrid data from.*/
	if(dbextfile(dbgrd,"gclgdisk",filename) <=0)
	{
		elog_notify(0,"%s Cannot find grid file named %s in directory %s\n",
			base_message,dfile,dir);
		throw 1;
	}
	fp = fopen(filename,"r");
	if(fp == NULL)
	{
		elog_notify(0,"%s file %s cannot be openned for read\n",
			base_message,filename);
		throw 1;
	}
	fseek(fp,foff,SEEK_SET);
	/* We alloc all memory first before reading so we can call a 
	free routine in case any reads fail */
	x1 = create_2dgrid_contiguous(n1,n2);
	x2 = create_2dgrid_contiguous(n1,n2); 
	x3 = create_2dgrid_contiguous(n1,n2); 
	lat = create_2dgrid_contiguous(n1,n2);
	lon = create_2dgrid_contiguous(n1,n2); 
	r = create_2dgrid_contiguous(n1,n2); 

	/* read data trapping read errors */
	gridsize = (n1)*(n2);
	if(fread(x1[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading x1\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(x2[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading x2\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(x3[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading x3\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(lat[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading latitude\n",
			base_message,read_error);
		throw 1;
	}
	if(fread(lon[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading longitude\n",
			base_message,read_error);
		throw 1;
	}

	if(fread(r[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,"%s %s reading radius\n",
			base_message,read_error);
		throw 1;
	}
}

/* The following two parallel functions are the inverse of the load
routines above.  That is, they take a db pointer and a pointer to 
a grid object and save the contents to external files and a database.
dir is the directory where the results are stored.  The file name
is always created from the gridname.  

Errors of any kind leave a warning message on the error log and 
invoke a throw.  This should be caught with a simple int handler.
-1 means total failure, positive count is number of i/o errors encounterd.

Author:  G Pavlis
Date:  August 2000
Modified:  December 2002
Converted to C++ with subsequent name change.  Little of the code
changed.
*/
void GCLgrid3d::dbsave(Dbptr db, char *dir)
{
	FILE *fp;
	char *filename;
	int fssize;
	int errcount=0;
	long int fpos;
	int foff;
	int gridsize;
	int dimensions=3;
	char *fwerr="fwrite error on file %s";
	char cdef[2]="d",gdef[2]="d";

	db = dblookup(db,0,"gclgdisk",0,0);
	if(db.record == dbINVALID)
	{
		elog_notify(0,"lookup failed for gclgdisk table.  Extension table probably not defined\n");
		throw -1;
	}
	/*Save the data first so that in the event of a failure we don't 		have to patch the database afterwards. */
	fssize = strlen(dir) + strlen(name) + 1;
	allot(char *,filename,fssize);
	
	strcpy(filename,dir);
	strcat(filename,"/");
	strcat(filename,name);
	fp = fopen(filename,"a+");
	if(fp==NULL)
	{
		elog_notify(0,"Cannot open output file %s\nNothing save\n",
			filename);
		throw -1;
	}
	fseek(fp,0,SEEK_END);
	/* The use of the int cast is unnecessary on some machines,
	but may be problematic on Solaris where they are switching from
	32 to 64.   This is safe if overkill*/
	fpos = ftell(fp);
	foff = (int)fpos;
	gridsize = (n1)*(n2)*(n3);
	if(fwrite(x1[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(x2[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(x3[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(lat[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(lon[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(r[0][0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	fclose(fp);
	/* Now we write a row in the database for this grid.  Note
	some quantities have to be converted from radians to degrees.*/
	if(dbaddv(db,0,
		"gridname",name,
		"dimensions",dimensions,
		"lat",deg(lat0),
		"lon",deg(lon0),
		"radius",r0,
		"azimuth_y",deg(azimuth_y),
		"dx1nom",dx1_nom,
		"dx2nom",dx2_nom,
		"dx3nom",dx3_nom,
		"n1",n1,
		"n2",n2,
		"n3",n3,
		"xlow",xlow,
		"xhigh",xhigh,
		"ylow",ylow,
		"yhigh",yhigh,
		"zlow",zlow,
		"zhigh",zhigh,
		"i0",i0,
		"j0",j0,
		"k0",k0,
		"cdefined",cdef,
		"geodefined",gdef,
		"datatype","t8",
		"dir",dir,
		"dfile",name,
		"foff",foff,
		0) < 0)
	{
		elog_notify(0,"dbaddv error for 3d grid into gclgdisk table\n");
		throw -1;
	}
	free(filename);
	if(errcount)throw errcount;
}
/* Parallel routine for 2d*/
void GCLgrid::dbsave(Dbptr db, char *dir)
{
	FILE *fp;
	char *filename;
	int fssize;
	int errcount=0;
	long int fpos;
	int foff;
	int gridsize;
	int dimensions=2;
	char *fwerr="fwrite error on file %s";
	char cdef[2]="d",gdef[2]="d";

	db = dblookup(db,0,"gclgdisk",0,0);
	if(db.record == dbINVALID)
	{
		elog_notify(0,"lookup failed for gclgdisk table.  Extension table probably not defined\n");
		throw -1;
	}
	/*Save the data first so that in the event of a failure we don't 		have to patch the database afterwards. */
	fssize = strlen(dir) + strlen(name) + 1;
	allot(char *,filename,fssize);
	
	strcpy(filename,dir);
	strcat(filename,"/");
	strcat(filename,name);
	fp = fopen(filename,"a+");
	if(fp==NULL)
	{
		elog_notify(0,"Cannot open output file %s\nNothing save\n",
			filename);
		throw -1;
	}
	fseek(fp,0,SEEK_END);
	/* The use of the int cast is unnecessary on some machines,
	but may be problematic on current versions of solaris so I'll
	be safe */
	fpos = ftell(fp);
	foff = (int)fpos;
	gridsize = (n1)*(n2);
	if(fwrite(x1[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(x2[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(x3[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(lat[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(lon[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	if(fwrite(r[0],sizeof(double),gridsize,fp) != gridsize)
	{
		elog_notify(0,fwerr,filename);
		++errcount;
	}
	fclose(fp);
	/* Now we write a row in the database for this grid*/
	if(dbaddv(db,0,
		"gridname",name,
		"dimensions",dimensions,
		"lat",deg(lat0),
		"lon",deg(lon0),
		"radius",r0,
		"azimuth_y",deg(azimuth_y),
		"dx1nom",dx1_nom,
		"dx2nom",dx2_nom,
		"n1",n1,
		"n2",n2,
		"xlow",xlow,
		"xhigh",xhigh,
		"ylow",ylow,
		"yhigh",yhigh,
		"zlow",zlow,
		"zhigh",zhigh,
		"i0",i0,
		"j0",j0,
		"cdefined",cdef,
		"geodefined",gdef,
		"datatype","t8",
		"dir",dir,
		"dfile",name,
		"foff",foff,
		0) < 0)
	{
		elog_notify(0,"dbaddv error for 2d grid into gclgdisk table\n");
		throw -1;
	}
	free(filename);
	if(errcount) throw errcount;
}
