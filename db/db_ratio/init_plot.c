#include <stdio.h>

int init_plot (char *fname)
{
	int itran=0;
	float ssize=0.9;
	float xwin=0.0;
	float ywin=0.0;
	char plotfile[256];
	char display[16];
	char program[1024];
	static int nplot=1;
        float xplt, yplt;
        float angle=0.0;
        int iclip=0;
        int iref=5;
        float height=0.08;
        float ratio=1.0;
        float slant=0.0;
        int jfont=114; 
	float xdim, ydim, xlow, ylow;
	float xmin, xmax, ymin, ymax;
	float thick=0.0;
	int ithick=0;
	float x1, y1, x2, y2;
	float hue, light, sat;
	float fac;
	long itime;
	int i;

 	if (fname) {
 		if (fname[0]) {
 			strcpy (plotfile, fname);
		} else {
			sprintf (plotfile, "%s.%d%d.ps", "db_ratio", getpid(), nplot++);
		}
	} else {
		sprintf (plotfile, "%s.%d%d.ps", "db_ratio", getpid(), nplot++);
	}
	strcpy (display, " ");
	strcpy (program, "db_ratio");
	initt_ (&itran, plotfile, display, program, &ssize,
		&xwin, &ywin, strlen(plotfile), strlen(display),
		strlen(program));
	ydim = 10.0;
	xdim = 10.0;
	ylow = 0.0;
	xlow = 0.0;
	xmin = 0.0;
	xmax = 1.0;
	ymin = 0.0;
	ymax = 1.0;
	setdim_ (&xdim, &ydim, &xlow, &ylow);
	setscl_ (&xmin, &xmax, &ymin, &ymax);
	iclip = 1;
	box_ (&xmin, &xmax, &ymin, &ymax, &thick, &ithick, &iclip);
	ydim = 9.8;
	xdim = 7.3;
	ylow = 0.2;
	xlow = 0.2;
	ymin = 0.0;
	ymax = 9.8;
	xmin = 0.0;
	xmax = 7.3;
	setdim_ (&xdim, &ydim, &xlow, &ylow);
	setscl_ (&xmin, &xmax, &ymin, &ymax);
        xplt = 0.0;
        yplt = 0.0;
	jfont = 113;
	height = 0.10;
	iref = 0;
        cfont_ (&jfont);
        chrsiz_ (&height, &ratio, &slant);
        text_ (&xplt, &yplt, &angle, &iref, "JSPC", &iclip, strlen("JSPC"));
	jfont = 115;
        cfont_ (&jfont);
        itime = time(NULL);
	sprintf (program, "%s %s %s %s", "db_ratio:", plotfile, cuserid(NULL), ctime(&itime));
	program[strlen(program)-1] = '\0';
	xplt = 0.5;
        text_ (&xplt, &yplt, &angle, &iref, program, &iclip, strlen(program));
}
