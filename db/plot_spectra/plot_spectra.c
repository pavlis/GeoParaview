#include <stdio.h>
#include <string.h>
#include "coords.h"
#include "stock.h"
#include "arrays.h"
#include "db.h"
#include "pf.h"
#include "spectrum.h"
/* plottype defines */
#define LINLIN 1
#define LOGLIN 2
#define LOGLOG 3
#define UNKNOWN 0
void usage()
{
        die(0,"Usage:  \nplot_spectrum db [-pf pffile]\n");
}
Spectrum read_spectrum(Dbptr);
main(int argc, char **argv)
{
	char *dbin;
	Dbptr db,dbsub;
	Spectrum *s;  /* array of spectra structures holding all data to plot */
	float *f;  /* frequencies */
	char ss[256];  /* Holds subset condition string */
	char test[4];  /* y or n test character for interactive questions */
	int nspec;  /* number of records */
	int i,j;
	Tbl *process_list;
	Pf *pf;
	char *pffile=NULL;

	/* plot variables */
        int plottype;
        float xmin=1.0,xmax=100.0;
        float ymin=0.01,ymax=1000000.0;
        float xdim=4.0,ydim=6.0;
        float xmarg=1.0, ymarg=1.0, xlow=1.5, ylow=1.5;
        char labelx[16]="Frequency (Hz)",labely[32],
                        title[32];
        int iclear=0;
        int roff=0,iclip=0;
	char psfile[256];
 
        float hue, light, sat;
        char *fmtx,*fmty;
        float dxsmal=5.0, dxnumb=20.0;
        float dysmal=5.0, dynumb=20.0;

	/* range variables */
	float fmin, fmax, pmin, pmax;
	char *ptstr;

	if(argc != 2)
	{
		usage();
		exit(1);
	}
	dbin = argv[1];
        if (dbopen (dbin, "r+", &db) == dbINVALID) {
                die(0,"%s:  Unable to open database '%s'.\n", argv[0], dbin);
        }
	for(i=2;i<argc;++i)
	{
		if(!strcmp(argv[i],"-pf"))
		{
			++i;
			pffile=argv[i];
		}
	}
	if(pffile==NULL) pffile=strdup("plot_spectra");
	if(pfread(pffile,&pf)) die(0,"pfread error\n");
	ptstr=pfget_string(pf,"plot_type");
	if(!strcmp(ptstr,"linlin"))
		plottype = LINLIN;
	else if(!strcmp(ptstr,"loglin"))
		plottype = LOGLIN;
	else if (!strcmp(ptstr,"loglog"))
		plottype = LOGLOG;
	else
		die(0,"Illegal plot_type parameter=%s\n  must be linlin, loglin, or loglog\n",
			ptstr);
	xmin = pfget_double(pf,"fmin");
	xmax = pfget_double(pf,"fmax");
	ymin = pfget_double(pf,"ymin");
	ymax = pfget_double(pf,"ymax");
	dysmal=pfget_double(pf,"yaxis_tic");
	dynumb=pfget_double(pf,"yaxis_annotate");
	fmty=pfget_string(pf,"yaxis_annotate_format");
	dxsmal=pfget_double(pf,"xaxis_tic");
	dxnumb=pfget_double(pf,"xaxis_annotate");
	fmtx=pfget_string(pf,"xaxis_annotate_format");

	process_list=pfget_tbl(pf,"dbprocess_list");
	db = dbprocess(db,process_list,0);
        dbquery(db,dbRECORD_COUNT,&nspec);
        if(nspec <= 0)
                die(0,"%s:  Cannot build working database view from db=%s\nCheck dbprocess_list Tbl in pf file\n",argv[0], dbin); 
	fprintf(stdout,"plot_spectra program\ndatabase %s open\npsdisc table has %d entries\n",
			dbin, nspec);
	do
	{
		fprintf(stdout,"Enter subset condition for database (enter none to omit subsetting\n");
		fgets(ss,256,stdin);
		if(strcmp(ss,"none"))
		{
			dbsub = dbsubset(db,ss,0);
		}
		else
		{
			dbsub =  dblookup(db,0,"psdisc",0,0);
		}
		dbquery(dbsub,dbRECORD_COUNT,&nspec);
		fprintf(stdout,"View has %d record\nType y to continue, or type any other character to change subset condition\n",
			nspec);
		fgets(test,4,stdin);
	} 
	while(test[0] != 'y');
	s = (Spectrum *) calloc(nspec,sizeof(Spectrum));
	if(s==NULL) die(1,"Cannot alloc %d element Spectrum array\n",nspec);
	for(dbsub.record=0,i=0;dbsub.record<nspec;++dbsub.record,++i)
	{
		s[i] = read_spectrum(dbsub);
		if(s[i].spec == NULL)
		{
			fprintf(stderr,"Warning:  error reading spectrum for %s %s %s %s %lf\n", 
				s[i].sta, s[i].chan, s[i].sname, s[i].time);
			fprintf(stderr,"Data from this spectrum skipped\n");
			--i;
		}
	}
	nspec = i;
	if(nspec <= 0) die(1,"No data to process after read loop\n");
	/* We need to scan all the data we have for frequency and spectral ranges */
	fmin = s[0].freq0;
	fmax = fmin;
	for(i=0;i<nspec;++i)
	{
		double f_test;
		fmin = MIN(fmin,s[i].freq0);
		f_test = (double) (s[i].nfreq - 1);
		f_test *= s[i].df;
		f_test += s[i].freq0;
		fmax = MAX(fmax,(float)f_test);
	}
	pmin = s[0].spec[0];
	pmax = pmin;
	for(i=0;i<nspec;++i)
        {
		for(j=0;j<s[i].nfreq;++j)
		{
			pmin = MIN(pmin,s[i].spec[j]);
			pmax = MAX(pmax,s[i].spec[j]);
		}
	}

	fprintf(stdout,"Enter plot title:");
	fscanf(stdin,"%s",title);
	fprintf(stdout,"Postscript output file (none for to omit):",psfile);
	fscanf(stdin,"%s",psfile);
	init_plot(psfile);


	switch(plottype)
	{
	case (LINLIN):
		axis_(&xdim,&ydim,&xmarg,&ymarg,
                                &xlow,&ylow,&xmax,&xmin,&ymax,&ymin,
				&dxsmal, &dxnumb, &dysmal, &dynumb,
				fmtx, fmty, labelx, labely, title, &iclear,
				strlen(fmtx), strlen(fmty),
                                strlen(labelx),strlen(labely),strlen(title));
		break;
	case (LOGLIN):
                lyaxis_(&xdim,&ydim,&xmarg,&ymarg,
                                &xlow,&ylow,&xmax,&xmin,&ymax,&ymin,
                                &dxsmal, &dxnumb,
                                fmtx, labelx, labely, title, &iclear,
                                strlen(fmtx),strlen(labelx),
                                strlen(labely),strlen(title));
		break;
	case (LOGLOG):
        	llaxis_(&xdim,&ydim,&xmarg,&ymarg,
			&xlow,&ylow,&xmax,&xmin,&ymax,&ymin,
                        labelx, labely, title,&iclear,
                        strlen(labelx),strlen(labely),strlen(title));
	}
	for(i=0;i<nspec;++i)
	{
		f = (float *) calloc(s[i].nfreq,sizeof(float));
		for(j=0;j<s[i].nfreq;++j)
                        f[j]=s[i].freq0+(float)j*s[i].df;
		plot1_(&(s[i].nfreq),f, s[i].spec,&roff, &iclip);
		free(f);
	}

	/* a bug in niceplot makes it necessary to replot the 
	last item to flush the buffer */
	i=nspec-1;
	f = (float *) calloc(s[i].nfreq,sizeof(float));
	for(j=0;j<s[i].nfreq;++j)
                        f[j]=s[i].freq0+(float)j*s[i].df;
        plot1_(&(s[i].nfreq),f, s[i].spec,&roff, &iclip);

	finitt_();
	exit(1);
}
