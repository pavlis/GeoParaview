#include <string>
#include <math.h>
#include <algorithm>

#include "elog.h"
#include "xplot.h"
#include "gclgrid.h"
#include "perf.h"
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include "Metadata.h"
#include "seispp.h"
#include "SeismicPlot.h"

#include <vector>

#define MAXPLOTSAMPLES 100000
namespace SEISPP
{
using namespace std;
using namespace SEISPP;

/*********************************************************/
static XImage *RotImage90(Display *dpy, XImage *oldImage)
{
	int	widthpad,x1,y1,x2,y2,i,nbpr;
	unsigned char	*bits;
	XImage	*image;
	int	width1 =		oldImage->width;
	int	width2 =		oldImage->height;
	int	height2 =	oldImage->width;
	int	scr =			DefaultScreen(dpy);

	widthpad = (1 + (width2-1)/(BitmapPad(dpy)/8))*BitmapPad(dpy)/8;
	nbpr = 1 + (widthpad-1)/8;
	bits = static_cast<unsigned char *>(calloc(nbpr*height2,sizeof(unsigned char)));
	if(bits==NULL)
		throw SeisppError("SeismicPlot::RotImage90:  Cannot alloc bitmap");
	image = XCreateImage(	(Display *) dpy,
				(Visual *) DefaultVisual(dpy,scr),
				(unsigned int) 1,
				(int) XYBitmap,
				(int) 0,
				(char *) bits,
				(unsigned int) widthpad,
				(unsigned int) height2,
				(int) BitmapPad(dpy),
				(int) nbpr);

	for (i = 0; i < nbpr*height2; i++)	bits[i]=0;
	for (x2 = 0; x2 < width2; x2++) {
		y1 = x2;
		for (y2 = 0; y2 < height2; y2++) {
			x1 = width1 - 1 - y2;
			XPutPixel(image,x2,y2,XGetPixel(oldImage,x1,y1));
		}
	}
	return image;
}
//
// These are functions extracted from SU xwigb 
//

			
/* update parameters associated with zoom box */
static void zoomBox (int x, int y, int w, int h, 
	int xb, int yb, int wb, int hb,
	float x1, float x2,
	float y1, float y2,
	float *x1b, float *x2b,
	float *y1b, float *y2b,
        int style)
{
	/* if width and/or height of box are zero, just copy values */
	if (wb==0 || hb==0) {
		*x1b = x1; *x2b = x2;
		*y1b = y1; *y2b = y2;
		return;		
	} 
	
	/* clip box */
	if (xb<x) {
		wb -= x-xb;
		xb = x;
	}
	if (yb<y) {
		hb -= y-yb;
		yb = y;
	}
	if (xb+wb>x+w) wb = x-xb+w;
	if (yb+hb>y+h) hb = y-yb+h;	
	
	/* determine box limits */
        if (style == SEISMIC) {
                *x1b = x1+(xb-x)*(x2-x1)/w;
                *x2b = x1+(xb+wb-x)*(x2-x1)/w;
                *y1b = y1+(yb-y)*(y2-y1)/h;
                *y2b = y1+(yb+hb-y)*(y2-y1)/h;
        } else {
                *x2b = x2+(yb-y)*(x1-x2)/h;
                *x1b = x2+(yb+hb-y)*(x1-x2)/h;
                *y1b = y1+(xb-x)*(y2-y1)/w;
                *y2b = y1+(xb+wb-x)*(y2-y1)/w;
        }
}

void xMouseLoc(Display *dpy, Window win, XEvent event, int style, Bool show,
	int x, int y, int width, int height,
	float x1begb, float x1endb, float x2begb, float x2endb,
	float p2beg, float p2end)
{
	static XFontStruct *fs=NULL;
	static XCharStruct overall;
	static GC gc;
	int dummy,xoffset=5,yoffset=5;
	float x1,x2;
	char string[256];

	/* if first time, get font attributes and make gc */
	if (fs==NULL) {
		fs = XLoadQueryFont(dpy,"fixed");
		gc = XCreateGC(dpy,win,0,NULL);

		/* make sure foreground/background are black/white */
		XSetForeground(dpy,gc,BlackPixel(dpy,DefaultScreen(dpy)));
		XSetBackground(dpy,gc,WhitePixel(dpy,DefaultScreen(dpy)));

		XSetFont(dpy,gc,fs->fid);
		overall.width = 1;
		overall.ascent = 1;
		overall.descent = 1;
	}

	/* erase previous string */
	XClearArea(dpy,win,xoffset,yoffset,
		overall.width,overall.ascent+overall.descent,False);

	/* if not showing, then return */
	if (!show) return;

	/* convert mouse location to (x1,x2) coordinates */
	if (style==NORMAL) {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
		x2 = p2end+x2endb+(p2beg+x2begb-x2endb-p2end)*
			(event.xmotion.y-y)/height;
	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = p2beg+x2begb+(p2end+x2endb-x2begb-p2beg)*
			(event.xmotion.x-x)/width;
	}

	/* draw string indicating mouse location */
	sprintf(string,"(%0.6g,%0.6g)",x1,x2);
	XTextExtents(fs,string,(int)strlen(string),&dummy,&dummy,&dummy,&overall);
	XDrawString(dpy,win,gc,xoffset,yoffset+overall.ascent,
		string,(int)strlen(string));
}
/* Derived from xMouseLoc above.  Instead of posting the pick value
the result is returned in PointPick structure. */

SeismicPick MouseLocate(Display *dpy, Window win, XEvent event, int style, 
        int x, int y, int width, int height,
        float x1begb, float x1endb, float x2begb, float x2endb,
        float p2beg, float p2end)
{
	float x1,x2;

	/* convert mouse location to (x1,x2) coordinates */
	if (style==NORMAL) {

		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
		x2 = p2end+x2endb+(p2beg+x2begb-x2endb-p2end)*
			(event.xmotion.y-y)/height;
		//Peng Wang, above commented out stuff is original
		//This is a temperary fix, fabs(x2begb-x2endb)+2 means imagine there are 10 traces,
		//but actually there are 11 spaces, so we can not simply do fabs(x2begb-x2endb)+1
		//which only yields 10, so one will get weird spacing. 
		//The assumption here is x2begb and x2endb are from 1 to nmembers, i.e., auto scaling
		//Most importantly, we need to return something between 1 to nmember.
//		x2 = (fabs(p2beg+x2begb-p2end-x2endb)+1)-(fabs(p2beg+x2begb-p2end-x2endb)+2)*((float)(event.xmotion.y-y))/((float)height);

	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = p2beg+x2begb+(p2end+x2endb-x2begb-p2beg)*
			(event.xmotion.x-x)/width;
//		x2=(fabs(p2end+x2endb-p2beg-x2begb)+1)-((fabs(p2end+x2endb-p2beg-x2begb)+2)*((float)(event.xmotion.x-x))/((float)width));
	}
	return(SeismicPick(x1,x2));

}


void xMousePrint(XEvent event, int style, FILE *mpicksfp,
		 int x, int y, int width, int height,
		 float x1begb, float x1endb, float x2begb, float x2endb,
		 float p2beg, float p2end, stringstream * ss)
{
	float x1,x2;

	/* convert mouse location to (x1,x2) coordinates */
	if (style==NORMAL) {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
		x2 = p2end+x2endb+(p2beg+x2begb-x2endb-p2end)*
			(event.xmotion.y-y)/height;
	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = p2beg+x2begb+(p2end+x2endb-x2begb-p2beg)*
			(event.xmotion.x-x)/width;
	}

	/* write string indicating mouse location */
	if (ss != NULL) (*ss) << x1<<" "<<x2<<endl;
	else 
	    fprintf(mpicksfp, "%0.6g  %0.6g\n", x1, x2);
}
/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */


/*********************** self documentation **********************/
/*****************************************************************************
RFWTVA - Rasterize a Float array as Wiggle-Trace-Variable-Area.

rfwtva	rasterize a float array as wiggle-trace-variable-area.

******************************************************************************
Function Prototype:
void rfwtva (int n, float z[], float zmin, float zmax, float zbase,
	int yzmin, int yzmax, int xfirst, int xlast,
	int wiggle, int nbpr, unsigned char *bits, int endian);

******************************************************************************
Input:
n		number of samples in array to rasterize
z		array[n] to rasterize
zmin		z values below zmin will be clipped
zmax		z values above zmax will be clipped
zbase		z values between zbase and zmax will be filled (see notes)
yzmin		horizontal raster coordinate corresponding to zmin
yzmax		horizontal raster coordinate corresponding to zmax
xfirst		vertical raster coordinate of z[0] (see notes)
xlast		vertical raster coordinate of z[n-1] (see notes)
wiggle		=0 for no wiggle (VA only); =1 for wiggle (with VA)
		wiggle 2<=wiggle<=5 for solid/grey coloring of VA option
                shade of grey: wiggle=2 light grey, wiggle=5 black
nbpr		number of bytes per row of bits
bits		pointer to first (top,left) byte in image
endian		byte order  =1 big endian  =0 little endian 

Output:
bits		pointer to first (top,left) byte in image

******************************************************************************
Notes:
The raster coordinate of the (top,left) bit in the image is (0,0).
In other words, x increases downward and y increases to the right.
Raster scan lines run from left to right, and from top to bottom.
Therefore, xfirst, xlast, yzmin, and yzmax should not be less than 0.
Likewise, yzmin and yzmax should not be greater than nbpr*8-1, and 
care should be taken to ensure that xfirst and xlast do not cause bits 
to be set outside (off the bottom) of the image. 

Variable area fill is performed on the right-hand (increasing y) side
of the wiggle.  If yzmin is greater than yzmax, then z values between
zmin will be plotted to the right of zmax, and z values between zbase
and zmin are filled.  Swapping yzmin and yzmax is an easy way to 
reverse the polarity of a wiggle.

The variable "endian" must have a value of 1 or 0. If this is
not a case an error is returned.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 07/01/89
MODIFIED:  Paul Michaels, Boise State University, 29 December 2000
           Added solid/grey shading scheme, wiggle>=2 option for peaks/troughs
*****************************************************************************/
/**************** end self doc ********************************/


void rfwtva (
	int n, float z[], float zmin, float zmax, float zbase,
	int yzmin, int yzmax, int xfirst, int xlast,
	int wiggle, int nbpr, unsigned char *bits, int endian)
/*****************************************************************************
Rasterize a float array as wiggle-trace-variable-area.
******************************************************************************
Input:
n		number of samples in array to rasterize
z		array[n] to rasterize
zmin		z values below zmin will be clipped
zmax		z values above zmax will be clipped
zbase		z values between zbase and zmax will be filled (see notes)
yzmin		horizontal raster coordinate corresponding to zmin
yzmax		horizontal raster coordinate corresponding to zmax
xfirst		vertical raster coordinate of z[0] (see notes)
xlast		vertical raster coordinate of z[n-1] (see notes)
wiggle		=0 for no wiggle (VA only); =1 for wiggle (with VA)
		wiggle 2<=wiggle<=5 for solid/grey coloring of VA option
                shade of grey: wiggle=2 light grey, wiggle=5 black
nbpr		number of bytes per row of bits
bits		pointer to first (top,left) byte in image

Output:
bits		pointer to first (top,left) byte in image
******************************************************************************
Notes:
The raster coordinate of the (top,left) bit in the image is (0,0).
In other words, x increases downward and y increases to the right.
Raster scan lines run from left to right, and from top to bottom.
Therefore, xfirst, xlast, yzmin, and yzmax should not be less than 0.
Likewise, yzmin and yzmax should not be greater than nbpr*8-1, and 
care should be taken to ensure that xfirst and xlast do not cause bits 
to be set outside (off the bottom) of the image. 

Variable area fill is performed on the right-hand (increasing y) side
of the wiggle.  If yzmin is greater than yzmax, then z values between
zmin will be plotted to the right of zmax, and z values between zbase
and zmin are filled.  Swapping yzmin and yzmax is an easy way to 
reverse the polarity of a wiggle.

The variable "endian" must have a value of 1 or 0. If this is
not a case an error is returned.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 07/01/89
Modified:  Craig Artley, Colorado School of Mines, 04/14/92
           Fixed bug in computing yoffset.  Previously, when zmin==zmax
           the rasterized trace was shifted to the left by one trace.
MODIFIED:  Paul Michaels, Boise State University, 29 December 2000
           Added solid/grey color scheme, wiggle=2 option for peaks/troughs
*****************************************************************************/
{
	int iscale,xscale,dx,dy,i,x,y,
		ymin,ymax,ybase,ythis,ynext,xthis,xnext,xstep;
	int igrey,ideci;
	float yscale,yoffset,zthis,znext;
	register int bit;
	register unsigned char *byte;

	/* if solid/grey coloring desired      */
	if (wiggle>=2)
	{  igrey=abs(wiggle); wiggle=1; }
	else
	{  igrey=0; }

	/* determine min and max y coordinates */
	ymin = (yzmin<yzmax)?yzmin:yzmax;
	ymax = (yzmax>yzmin)?yzmax:yzmin;

	/* restrict min and max y coordinates */
	ymin = (ymin>0)?ymin:0;
	ymax = (ymax<nbpr*8-1)?ymax:nbpr*8-1;
	
	/* determine sample index scale factor */
	iscale = n-1;
	
	/* determine y scale factor and offset */
	yscale = (zmax!=zmin)?(yzmax-yzmin)/(zmax-zmin):1.0;
	yoffset = (zmax!=zmin)?yzmin-zmin*yscale:0.5*(yzmin+yzmax);
	
	/* determine x scale factor and step */
	xscale = (n>1)?xlast-xfirst:0;
	xstep = (xlast>xfirst)?1:-1;
	
	/* determine base y coordinate */
	ybase = static_cast<int>(yoffset+zbase*yscale);
	ybase = (ybase>ymin)?ybase:ymin;
	ybase = (ybase<ymax)?ybase:ymax;
	
	/* initialize next values of x, y, and z */
	znext = *z;
	ynext = static_cast<int>(yoffset+znext*yscale);
	xnext = xfirst;
	
	/* loop over samples */
	for (i=0; i<n; i++,z++) {
		
		/* determine x coordinate for this sample */
		xthis = xnext;
		
		/* determine x coordinate for next sample */
		xnext = (i<iscale)?xfirst+(i+1)*xscale/iscale:xthis+xstep;

		/* skip sample if next sample falls on same x coordinate */
		if (xnext==xthis) continue;
		
		/* determine difference in x coordinates */
		dx = xnext-xthis;
		
		/* determine this sample value */
		zthis = znext;
		
		/* determine next sample value */
		znext = (i<n-1)?*(z+1):zthis;
		
		/* determine y coordinate for this sample */
		ythis = ynext;
		
		/* determine y coordinate for next sample */
		ynext = static_cast<int>(yoffset+znext*yscale);
		
		/* determine difference in y coordinates */
		dy = ynext-ythis;
		
		/* loop over x coordinates */
		for (x=xthis,y=ythis; x!=xnext;
			x+=xstep,y=ythis+(x-xthis)*dy/dx) {
			
			/* apply clip */
			if (y<ymin) y = ymin;
			if (y>ymax) y = ymax;
			
			/* determine the bit and byte */
			/* original: bit = 7-y&7; */
			bit = (7-y)&7;

			byte = bits+x*nbpr+(y>>3);

			/* if wiggle or filling, then set the bit */
			if (wiggle || y>ybase) { 
				if (endian==0) 
					*byte |= 1<<(-bit+7);
				else if (endian==1)
					*byte |= 1<<bit;
				else
					fprintf(stderr,"endian must equal either 0 or 1\n");
			}

			
			/* while y greater than base, set more bits (SOLID FILL PEAKS) */
			while (y>ybase) {
				y-=1;
				bit+=1;
				if (bit>=8) {
					byte--;
					bit = 0;
				}
				if (endian==0)
					*byte |= 1<<(-bit+7);
				else if (endian==1)
					*byte |= 1<<bit;
				else
					fprintf(stderr,"endian must equal either 0 or 1\n");
			}  /* endwhile */

			/* while y less than base, set more bits (GREY FILL TROUGHS) */

			if (igrey>0)
			{
			ideci=6-igrey;
			if (ideci<1) ideci=1;
			
				while (y<ybase) {
					y+=ideci;
					bit-=ideci;
					if (bit<0) {
						byte++;
						bit = 7;
					}
					if (endian==0)
						*byte |= 1<<(-bit+7);
					else if (endian==1)
						*byte |= 1<<bit;
					else
						fprintf(stderr,"endian must equal either 0 or 1\n");
				}  /* endwhile  */
			}  /*  endif igrey   */

		}  /* next x  */
	}   /* next sample  */
}   /* end rfwtva   */
/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
SINC - Return SINC(x) for as floats or as doubles

fsinc		return float value of sinc(x) for x input as a float
dsinc		return double precision sinc(x) for double precision x

******************************************************************************
Function Prototype:
double dsinc (double x);

******************************************************************************
Input:
x		value at which to evaluate sinc(x)

Returned: 	sinc(x)

******************************************************************************
Notes:
    sinc(x) = sin(PI*x)/(PI*x) 

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/


double dsinc (double x)
/*****************************************************************************
Return sinc(x) = sin(PI*x)/(PI*x) (double version)
******************************************************************************
Input:
x		value at which to evaluate sinc(x)

Returned:	sinc(x)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	double pix;

	if (x==0.0) {
		return 1.0;
	} else {
		pix = M_PI*x;
		return sin(pix)/pix;
	}
}
/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
STOEP - Functions to solve a symmetric Toeplitz linear system of equations
	 Rf=g for f

stoepd		solve a symmetric Toeplitz system - doubles
stoepf		solve a symmetric Toeplitz system - floats

******************************************************************************
Function Prototypes:
void stoepd (int n, double r[], double g[], double f[], double a[]);

******************************************************************************
Input:
n		dimension of system
r		array[n] of top row of Toeplitz matrix
g		array[n] of right-hand-side column vector

Output:
f		array[n] of solution (left-hand-side) column vector
a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)

******************************************************************************
Notes:
These routines do NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/


void stoepd (int n, double r[], double g[], double f[], double a[])
/*****************************************************************************
Solve a symmetric Toeplitz linear system of equations Rf=g for f
(double version)
******************************************************************************
Input:
n		dimension of system
r		array[n] of top row of Toeplitz matrix
g		array[n] of right-hand-side column vector

Output:
f		array[n] of solution (left-hand-side) column vector
a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
******************************************************************************
Notes:
This routine does NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int i,j;
	double v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}

/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
MKSINC - Compute least-squares optimal sinc interpolation coefficients.

mksinc		Compute least-squares optimal sinc interpolation coefficients.

******************************************************************************
Function Prototype:
void mksinc (float d, int lsinc, float sinc[]);

******************************************************************************
Input:
d		fractional distance to interpolation point; 0.0<=d<=1.0
lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20

Output:
sinc		array[lsinc] containing interpolation coefficients

******************************************************************************
Notes:
The coefficients are a least-squares-best approximation to the ideal
sinc function for frequencies from zero up to a computed maximum
frequency.  For a given interpolator length, lsinc, mksinc computes
the maximum frequency, fmax (expressed as a fraction of the nyquist
frequency), using the following empirically derived relation (from
a Western Geophysical Technical Memorandum by Ken Larner):

	fmax = min(0.066+0.265*log(lsinc),1.0)

Note that fmax increases as lsinc increases, up to a maximum of 1.0.
Use the coefficients to interpolate a uniformly-sampled function y(i) 
as follows:

            lsinc-1
    y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
              j=0

Interpolation error is greatest for d=0.5, but for frequencies less
than fmax, the error should be less than 1.0 percent.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/


void mksinc (float d, int lsinc, float sinc[])
/*****************************************************************************
Compute least-squares optimal sinc interpolation coefficients.
******************************************************************************
Input:
d		fractional distance to interpolation point; 0.0<=d<=1.0
lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20

Output:
sinc		array[lsinc] containing interpolation coefficients
******************************************************************************
Notes:
The coefficients are a least-squares-best approximation to the ideal
sinc function for frequencies from zero up to a computed maximum
frequency.  For a given interpolator length, lsinc, mksinc computes
the maximum frequency, fmax (expressed as a fraction of the nyquist
frequency), using the following empirically derived relation (from
a Western Geophysical Technical Memorandum by Ken Larner):

	fmax = min(0.066+0.265*log(lsinc),1.0)

Note that fmax increases as lsinc increases, up to a maximum of 1.0.
Use the coefficients to interpolate a uniformly-sampled function y(i) 
as follows:

            lsinc-1
    y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
              j=0

Interpolation error is greatest for d=0.5, but for frequencies less
than fmax, the error should be less than 1.0 percent.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int j;
	double s[20],a[20],c[20],work[20],fmax;

	/* compute auto-correlation and cross-correlation arrays */
	fmax = 0.066+0.265*std::log((double)lsinc);
	fmax = (fmax<1.0)?fmax:1.0;
	for (j=0; j<lsinc; j++) {
		a[j] = dsinc(fmax*j);
		c[j] = dsinc(fmax*(lsinc/2-j-1+d));
	}

	/* solve symmetric Toeplitz system for the sinc approximation */
	stoepd(lsinc,a,c,s,work);
	for (j=0; j<lsinc; j++)
		sinc[j] = s[j];
}
/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
INTTABLE8 -  Interpolation of a uniformly-sampled complex function y(x)
		via a table of 8-coefficient interpolators

intt8r	interpolation of a uniformly-sampled real function y(x) via a
		table of 8-coefficient interpolators

******************************************************************************
Function Prototype:
void intt8r (int ntable, float table[][8],
	int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[]);

******************************************************************************
Input:
ntable		number of tabulated interpolation operators; ntable>=2
table		array of tabulated 8-point interpolation operators
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array of x values at which y(x) is output

Output:
yout		array of output y(x) values:  yout[0] = y(xout[0]), etc.

******************************************************************************
NOTES:
ntable must not be less than 2.

The table of interpolation operators must be as follows:

Let d be the distance, expressed as a fraction of dxin, from a particular
xout value to the sampled location xin just to the left of xout.  Then,
for d = 0.0,

table[0][0:7] = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.
Likewise, for d = 1.0,

table[ntable-1][0:7] = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.  In general,
for d = (float)itable/(float)(ntable-1), table[itable][0:7] are the
weights applied to the 8 input samples nearest xout.  If the actual sample
distance d does not exactly equal one of the values for which interpolators
are tabulated, then the interpolator corresponding to the nearest value of
d is used.

Because extrapolation of the input function y(x) is defined by the left
and right values yinl and yinr, the xout values are not restricted to lie
within the range of sample locations defined by nxin, dxin, and fxin.

******************************************************************************
AUTHOR:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/


void intt8r (int ntable, float table[][8],
	int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
	int nxout, float xout[], float yout[])
{
	int ioutb,nxinm8,ixout,ixoutn,kyin,ktable,itable;
	float xoutb,xoutf,xouts,xoutn,frac,fntablem1,yini,sum,
		*yin0,*table00,*pyin,*ptable;

	/* compute constants */
	ioutb = -3-8;
	xoutf = fxin;
	xouts = 1.0/dxin;
	xoutb = 8.0-xoutf*xouts;
	fntablem1 = (float)(ntable-1);
	nxinm8 = nxin-8;
	yin0 = &yin[0];
	table00 = &table[0][0];

	/* loop over output samples */
	for (ixout=0; ixout<nxout; ixout++) {

		/* determine pointers into table and yin */
		xoutn = xoutb+xout[ixout]*xouts;
		ixoutn = (int)xoutn;
		kyin = ioutb+ixoutn;
		pyin = yin0+kyin;
		frac = xoutn-(float)ixoutn;
		ktable = frac>=0.0?frac*fntablem1+0.5:(frac+1.0)*fntablem1-0.5;
		ptable = table00+ktable*8;
		
		/* if totally within input array, use fast method */
		if (kyin>=0 && kyin<=nxinm8) {
			yout[ixout] = 
				pyin[0]*ptable[0]+
				pyin[1]*ptable[1]+
				pyin[2]*ptable[2]+
				pyin[3]*ptable[3]+
				pyin[4]*ptable[4]+
				pyin[5]*ptable[5]+
				pyin[6]*ptable[6]+
				pyin[7]*ptable[7];
		
		/* else handle end effects with care */
		} else {
	
			/* sum over 8 tabulated coefficients */
			for (itable=0,sum=0.0; itable<8; itable++,kyin++) {
				if (kyin<0)
					yini = yinl;
				else if (kyin>=nxin)
					yini = yinr;
				else
					yini = yin[kyin];
				sum += yini*(*ptable++);
			}
			yout[ixout] = sum;
		}
	}
}
/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
INTSINC8 - Functions to interpolate uniformly-sampled data via 8-coeff. sinc
		approximations:

ints8r	Interpolation of a uniformly-sampled real function y(x) via a
		table of 8-coefficient sinc approximations

******************************************************************************
Function Prototypes:
void ints8r (int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[]);

******************************************************************************
Input:
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array[nxout] of x values at which y(x) is output

Output:
yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.

******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.

The maximum error for frequiencies less than 0.6 nyquist is less than
one percent.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/


/* these are used by ints8r */
#define LTABLE 8
#define NTABLE 513

void ints8r (int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[])
/*****************************************************************************
Interpolation of a uniformly-sampled real function y(x) via a
table of 8-coefficient sinc approximations; maximum error for frequiencies
less than 0.6 nyquist is less than one percent.
******************************************************************************
Input:
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array[nxout] of x values at which y(x) is output

Output:
yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/

{
	static float table[NTABLE][LTABLE];
	static int tabled=0;
	int jtable;
	float frac;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled) {
		for (jtable=1; jtable<NTABLE-1; jtable++) {
			frac = (float)jtable/(float)(NTABLE-1);
			mksinc(frac,LTABLE,&table[jtable][0]);
		}
		for (jtable=0; jtable<LTABLE; jtable++) {
			table[0][jtable] = 0.0;
			table[NTABLE-1][jtable] = 0.0;
		}
		table[0][LTABLE/2-1] = 1.0;
		table[NTABLE-1][LTABLE/2] = 1.0;
		tabled = 1;
	}

	/* interpolate using tabulated coefficients */
	intt8r(NTABLE,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,xout,yout);
}

/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */


/*********************** self documentation **********************/
/*****************************************************************************
RFWTVAINT - Rasterize a Float array as Wiggle-Trace-Variable-Area, with
	    8 point sinc INTerpolation.

rfwtvaint	rasterize a float array as wiggle-trace-variable-area, and
		apply sinc interploation for display purposes.

******************************************************************************
Function Prototype:
void rfwtvaint (int n, float z[], float zmin, float zmax, float zbase,
	int yzmin, int yzmax, int xfirst, int xlast,
	int wiggle, int nbpr, unsigned char *bits, int endian);

******************************************************************************
Input:
n		number of samples in array to rasterize
z		array[n] to rasterize
zmin		z values below zmin will be clipped
zmax		z values above zmax will be clipped
zbase		z values between zbase and zmax will be filled (see notes)
yzmin		horizontal raster coordinate corresponding to zmin
yzmax		horizontal raster coordinate corresponding to zmax
xfirst		vertical raster coordinate of z[0] (see notes)
xlast		vertical raster coordinate of z[n-1] (see notes)
wiggle		=0 for no wiggle (VA only); =1 for wiggle (with VA)
                wiggle 2<=wiggle<=5 for solid/grey coloring of VA option
                shade of grey: wiggle=2 light grey, wiggle=5 black
nbpr		number of bytes per row of bits
bits		pointer to first (top,left) byte in image
endian		byte order  =1 big endian  =0 little endian 

Output:
bits		pointer to first (top,left) byte in image

******************************************************************************
Notes:
The raster coordinate of the (top,left) bit in the image is (0,0).
In other words, x increases downward and y increases to the right.
Raster scan lines run from left to right, and from top to bottom.
Therefore, xfirst, xlast, yzmin, and yzmax should not be less than 0.
Likewise, yzmin and yzmax should not be greater than nbpr*8-1, and 
care should be taken to ensure that xfirst and xlast do not cause bits 
to be set outside (off the bottom) of the image. 

Variable area fill is performed on the right-hand (increasing y) side
of the wiggle.  If yzmin is greater than yzmax, then z values between
zmin will be plotted to the right of zmax, and z values between zbase
and zmin are filled.  Swapping yzmin and yzmax is an easy way to 
reverse the polarity of a wiggle.

The variable "endian" must have a value of 1 or 0. If this is
not a case an error is returned.

The interpolation is by the 8 point sinc interpolation routine s8r.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 07/01/89
	Memorial University of Newfoundland: Tony Kocurko, Sept 1995.
	 Added sinc interpolation.
MODIFIED: Paul Michaels, Boise State University, 29 December 2000
          added solid/grey color scheme for peaks/troughs  wiggle=2 option
*****************************************************************************/
/**************** end self doc ********************************/


void rfwtvaint (
	int n, float z[], float zmin, float zmax, float zbase,
	int yzmin, int yzmax, int xfirst, int xlast,
	int wiggle, int nbpr, unsigned char *bits, int endian)
/*****************************************************************************
Rasterize a float array as wiggle-trace-variable-area.
******************************************************************************
Input:
n		number of samples in array to rasterize
z		array[n] to rasterize
zmin		z values below zmin will be clipped
zmax		z values above zmax will be clipped
zbase		z values between zbase and zmax will be filled (see notes)
yzmin		horizontal raster coordinate corresponding to zmin
yzmax		horizontal raster coordinate corresponding to zmax
xfirst		vertical raster coordinate of z[0] (see notes)
xlast		vertical raster coordinate of z[n-1] (see notes)
wiggle		=0 for no wiggle (VA only); =1 for wiggle (with VA)
                wiggle 2<=wiggle<=5 for solid/grey coloring of VA option
                shade of grey: wiggle=2 light grey, wiggle=5 black
nbpr		number of bytes per row of bits
bits		pointer to first (top,left) byte in image

Output:
bits		pointer to first (top,left) byte in image
******************************************************************************
Notes:
The raster coordinate of the (top,left) bit in the image is (0,0).
In other words, x increases downward and y increases to the right.
Raster scan lines run from left to right, and from top to bottom.
Therefore, xfirst, xlast, yzmin, and yzmax should not be less than 0.
Likewise, yzmin and yzmax should not be greater than nbpr*8-1, and 
care should be taken to ensure that xfirst and xlast do not cause bits 
to be set outside (off the bottom) of the image. 

Variable area fill is performed on the right-hand (increasing y) side
of the wiggle.  If yzmin is greater than yzmax, then z values between
zmin will be plotted to the right of zmax, and z values between zbase
and zmin are filled.  Swapping yzmin and yzmax is an easy way to 
reverse the polarity of a wiggle.

The variable "endian" must have a value of 1 or 0. If this is
not a case an error is returned.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 07/01/89
Modified:  Craig Artley, Colorado School of Mines, 04/14/92
           Fixed bug in computing yoffset.  Previously, when zmin==zmax
           the rasterized trace was shifted to the left by one trace.
MODIFIED: Paul Michaels, Boise State University, 29 December 2000
          added solid/grey color scheme for peaks/troughs  wiggle=2 option
*****************************************************************************/
{
	int i,y,
		ymin,ymax,ybase;
	int igrey,ideci;
	float yscale,yoffset,zthis;
	register int bit;
	register unsigned char *byte;

	static float *xout  , *yout;
	static int    nx = 0;
	float        *yin   ,  fxin , yinl, yinr, dxin, x0, deltax;
	int           nxin  ,  nxout;

	/* if solid/grey coloring desired      */
	if (wiggle>=2)
	{  igrey=abs(wiggle); wiggle=1; }
	else
	{  igrey=0; }

	/* Compute the number of raster scan lines. */
	nxout = ABS(xlast - xfirst + 1);

	/* If the # of scan lines has increased from a previous call,
	   allocate more*/
	if ( nxout > nx ) {

	  /* If a previous call allocated space for output values,
	     free them. */
	  if ( nx > 0 ) {
	    free (xout);
	    free (yout);
	  }

	  /* Allocate space for the scan line x values and interpolated
	     z values. */

		xout = (float *)calloc ((size_t)nxout, sizeof(float));
		yout = (float *)calloc ((size_t)nxout, sizeof(float));
	  nx   = nxout;
	}

	nxin   = n       ; /* There are n z-values.  */
	dxin   = 1.0     ; /* We go from index 0 to index n - 1 in steps
			      of 1.0  */
	fxin   = 0.0     ; /* The first index is 0.  */
	yin    = z       ; /* The input array is the z array.  */
	yinl   = z[0]    ; /* Set the values to the left of the array
			      to z[0]. */
	yinr   = z[n - 1]; /* Set the values to the right of the array
			      to z[n-1].*/

	deltax = (float)(n - 1) / (float)(nxout - 1);
	x0     =  0.0;
	if ( xfirst > xlast ) { /* If the z array is to be output backwards, */
		x0 = (float)(n - 1);  /* Then the first output index is n-1, */
		deltax = -deltax;     /*   and we decrement rather than
					   increment. */
	}
	for (i = 0; i < nxout; i++) /* Load the indices of the output values.*/
		xout[i] = x0 + (float)i * deltax;

	ints8r (nxin, dxin, fxin, yin, yinl, yinr, nxout, xout, yout);
		

	/* determine min and max y coordinates */
	ymin = (yzmin<yzmax)?yzmin:yzmax;
	ymax = (yzmax>yzmin)?yzmax:yzmin;

	/* restrict min and max y coordinates */
	ymin = (ymin>0)?ymin:0;
	ymax = (ymax<nbpr*8-1)?ymax:nbpr*8-1;
	
	/* determine y scale factor and offset */
	yscale = (zmax!=zmin)?(yzmax-yzmin)/(zmax-zmin):1.0;
	yoffset = (zmax!=zmin)?yzmin-zmin*yscale:0.5*(yzmin+yzmax);
	
	/* determine base y coordinate */
	ybase = static_cast<int>(yoffset+zbase*yscale);
	ybase = (ybase>ymin)?ybase:ymin;
	ybase = (ybase<ymax)?ybase:ymax;
	
	/* loop over scan lines */
	for (i = 0; i < nxout; i++) {
		zthis = yout[i];
		y     = static_cast<int>(yoffset + zthis * yscale);

		/* apply clip */
		if (y < ymin) y = ymin;
		if (y > ymax) y = ymax;
			
		/* determine the bit and byte */
		/* original: bit = 7-y&7; */
		bit = (7-y)&7;

		/* Tony Kocurko: Had been "bits+x*nbpr+(y>>3)".*/
		byte = bits+i*nbpr+(y>>3);

		/* if wiggle or filling, then set the bit */
		if (wiggle || y>ybase) {
			if (endian==0)
				*byte |= 1<<(-bit+7);
			else if (endian==1)
				*byte |= 1<<bit;
			else
				fprintf(stderr,"endian must equal either 0 or 1\n");
		}

		
		/* while y greater than base, set more bits (SOLID FILL PEAKS) */
		while (y>ybase) {
			y-=1;
			bit+=1;
			if (bit>=8) {
				byte--;
				bit = 0;
			}
			if (endian==0)
				*byte |= 1<<(-bit+7);
			else if (endian==1)
				*byte |= 1<<bit;
			else
				fprintf(stderr,"endian must equal either 0 or 1\n");
		}  /*  endwhile  */

		/* while y less than base, set more bits (GREY FILL TROUGHS) */
	        if (igrey>0)
	        {
                ideci=6-igrey;
                if (ideci<1) ideci=1;
		
			while (y<ybase) {
				y+=ideci;
				bit-=ideci;
				if (bit<0) {
					byte++;
					bit = 7;
				}
				if (endian==0)
					*byte |= 1<<(-bit+7);
				else if (endian==1)
					*byte |= 1<<bit;
				else
					fprintf(stderr,"endian must equal either 0 or 1\n");
			}  /* endwhile */
		}  /* endif igrey  */

	} /* next scan line  */
}  /*   end rfwtvaint   */

/* return pointer to new image bitmap of rasterized wiggles */
static XImage *newBitmap (Display *dpy, int width, int height,
	int n1, float d1, float f1, int n2, float *x2, float *z,
	float x1beg, float x1end, float x2beg, float x2end,
	float xcur, float clip, int wt, int va,
	float *p2begp, float *p2endp, int endian, int interp,
	int wigclip, int style)
{
	int widthpad,nbpr,i1beg,i1end,if1r,n1r,b1fz,b1lz,i2,i,n2in;
	float x2min,x2max,p2beg,p2end,bscale,boffset,bxcur,bx2;
	unsigned char *bits;
	int scr=DefaultScreen(dpy);
	XImage *image,*image2;
	float	x2margin,clip1,clip2;
	int	bx1max,bx2min,bx2max,b2f,b2l;
	int	width1,height1;


	/* determine bitmap dimensions and allocate space for bitmap */
	width1 =  (style==SEISMIC) ? width : height;
	height1 = (style==SEISMIC) ? height : width;
	widthpad = (1+(width1-1)/(BitmapPad(dpy)/8))*BitmapPad(dpy)/8;
	nbpr = 1+(widthpad-1)/8;
	bits = static_cast<unsigned char *>(calloc(nbpr*height1,sizeof(unsigned char)));
	if(bits==NULL)
		throw SeisppError("SeismicPlot::newBitmap:  allocation failure for bitmap array");
	
	for (i=0; i<nbpr*height1; ++i) bits[i] = 0;

	/* determine number of traces that fall within axis 2 bounds */
	x2min = MIN(x2beg,x2end);
	x2max = MAX(x2beg,x2end);
	for (i2=0,n2in=0; i2<n2; i2++)
		if (x2[i2]>=x2min && x2[i2]<=x2max) n2in++;

	/* determine pads for wiggle excursion along axis 2 */
	xcur = fabs(xcur);
	if (n2in>1) xcur *= (x2max-x2min)/(n2in-1);
	x2margin = (wigclip && n2in>1) ? (x2max-x2min)/(2*(n2in-1)) : xcur;
	p2beg = (x2end>=x2beg)?-x2margin:x2margin;
	p2end = -p2beg;

	bx2min = 0;
	bx2max = width1 - 1;
	bx1max = height1 - 1;

	/* determine scale and offset to map x2 units to bitmap units */
	bscale = bx2max/(x2end+p2end-x2beg-p2beg);
	boffset = -(x2beg+p2beg)*bscale;
	bxcur = xcur*bscale;

	/* adjust x1beg and x1end to fall on sampled values */
	i1beg = nint((x1beg-f1)/d1);
	i1beg = MAX(0,MIN(n1-1,i1beg));
	x1beg = f1+i1beg*d1;
	i1end = nint((x1end-f1)/d1);
	i1end = MAX(0,MIN(n1-1,i1end));
	x1end = f1+i1end*d1;

	/* determine first sample and number of samples to rasterize */
	if1r = MIN(i1beg,i1end);
	n1r = MAX(i1beg,i1end)-if1r+1;

	/* determine bits corresponding to first and last samples */
	b1fz = (x1end > x1beg) ? 0 : bx1max;
	b1lz = (x1end > x1beg) ? bx1max : 0;

	/* rasterize traces */
	for (i2=0; i2<n2; i2++,z+=n1) {

		/* skip traces not in bounds */
		if (x2[i2]<x2min || x2[i2]>x2max) continue;

		/* determine bitmap coordinate of trace */
		bx2 = boffset+x2[i2]*bscale;
		b2f = (int)(bx2-bxcur);
		b2l = (int)(bx2+bxcur);
		clip1 = -clip;
		clip2 = clip;
		if (b2f < bx2min) {
			clip1 *= ((bx2-bx2min) / bxcur);
			b2f = bx2min;
		}
		if (b2l > bx2max) {
			clip2 *= ((bx2max-bx2) / bxcur);
			b2l = bx2max;
		}

		/* rasterize one trace */
		if (interp==0) { /* don't use interpolation */
			rfwtva(n1r,&z[if1r],clip1,clip2,va?0:clip2,
				b2f,b2l,b1fz,b1lz,
				wt,nbpr,bits,endian);
		} else { /* use 8 point sinc interpolation */
			rfwtvaint(n1r,&z[if1r],clip1,clip2,va?0:clip2,
				b2f,b2l,b1fz,b1lz,
				wt,nbpr,bits,endian);
		}
		
	}
	
	/* return axis 2 pads */
	*p2begp = p2beg;  *p2endp = p2end;
	
	/* get pointer to image */
	image = XCreateImage(	(Display *) dpy,
				(Visual *) DefaultVisual(dpy,scr),
				(unsigned int) 1,
				(int) XYBitmap,
				(int) 0,
				(char *) bits,
				(unsigned int) widthpad,
				(unsigned int) height1,
				(int) BitmapPad(dpy),
				(int) nbpr);

	if (style == NORMAL) {
		image2 = RotImage90(dpy,image);
		XDestroyImage(image);
		image = image2;
	}

	return image;
}	


/* This is a relation function with the GUI, this is in case there is a scroll
   bar in GUI, we have to automatic slide the scroll bar is the area we want
   to go to is not in the visible area, that is why we have this function.
   This is also the only place in libseisppplot that we uses Motif instead of
   plain old Xlib */

/* Peng, I had to remove this for now.  Not clear how it should enter
in the interface */
/****************************************************************
void adjust_scrollbar(Display * dpy, Window win, Widget wdgt, int bx, int by, SeismicPlot * sptr)
{
            //scroll bars
            int xval, xincr, xpg_incr, xsz, xmax, xmin, ymax, ymin;
            int yval, ysz, yincr, ypg_incr;
	    int n1;
            Arg args1[5];
            Widget hsb, vsb;
	    Dimension vwidth, vheight;

            n1 = 0;
            XtSetArg (args1[n1], XmNhorizontalScrollBar, &hsb);  n1++;
            XtSetArg (args1[n1], XmNverticalScrollBar, &vsb);  n1++;
            XtGetValues (wdgt, args1, n1);
            n1 = 0;
            XtSetArg (args1[n1], XmNwidth, &vwidth);  n1++;
            XtSetArg (args1[n1], XmNminimum, &xmin); n1++;
            XtSetArg (args1[n1], XmNmaximum, &xmax); n1++;
            XtGetValues (hsb, args1, n1);
            n1 = 0;
            XtSetArg (args1[n1], XmNheight, &vheight);  n1++;
            XtSetArg (args1[n1], XmNminimum, &ymin); n1++;
            XtSetArg (args1[n1], XmNmaximum, &ymax); n1++;
            XtGetValues (vsb, args1, n1);

            XmScrollBarGetValues(hsb,&xval,&xsz,&xincr,&xpg_incr);
            XmScrollBarGetValues(vsb,&yval,&ysz,&yincr,&ypg_incr);

//	    cout << bx <<" "<<by<<" "<<xval<<" "<<vwidth<<endl;

                //Oddly, the value of the scrollbar appears to be the
                //left/top starting location of the viewport.

            if (bx >= xval+vwidth) {
                xval=xval+xincr;
                if (xval >= xmax-xsz) xval=xmax-xsz;
                XmScrollBarSetValues(hsb,xval,xsz,xincr,xpg_incr,1);
//		XClearArea(dpy,win,0,0,0,0,True);
		sptr->render_image(dpy,win);
            } else if (bx <= xval) {
                xval=xval-xincr;
                if (xval <= xmin) xval=xmin;
                XmScrollBarSetValues(hsb,xval,xsz,xincr,xpg_incr,1);
//		XClearArea(dpy,win,0,0,0,0,True);
		sptr->render_image(dpy,win);
            } else if (by >= yval+vheight) {
                yval=yval+yincr;
                if (yval > ymax-ysz) yval=ymax-ysz;
                XmScrollBarSetValues(vsb,yval,ysz,yincr,ypg_incr,1);
//		XClearArea(dpy,win,0,0,0,0,True);
		sptr->render_image(dpy,win);
            } else if (by <= yval) {
                yval=yval-yincr;
                if (yval < ymin) yval=ymin;
                XmScrollBarSetValues(vsb,yval,ysz,yincr,ypg_incr,1);
//		XClearArea(dpy,win,0,0,0,0,True);
		sptr->render_image(dpy,win);
            }

}
****************************************************************/


/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
RUBBERBOX -  Function to draw a rubberband box in X-windows plots

xRubberBox	Track pointer with rubberband box

******************************************************************************
Function Prototype:
void xRubberBox (Display *dpy, Window win, XEvent event,
	int *x, int *y, int *width, int *height);

******************************************************************************
Input:
dpy		display pointer
win		window ID
event		event of type ButtonPress

Output:
x		x of upper left hand corner of box in pixels
y		y of upper left hand corner of box in pixels
width		width of box in pixels
height		height of box in pixels

******************************************************************************
Notes:
xRubberBox assumes that event is a ButtonPress event for the 1st button;
i.e., it tracks motion of the pointer while the 1st button is down, and
it sets x, y, w, and h and returns after a ButtonRelease event for the
1st button.

Before calling xRubberBox, both ButtonRelease and Button1Motion events 
must be enabled.

This is the same rubberbox.c as in Xtcwp/lib, only difference is
that xRubberBox here is XtcwpRubberBox there, and a shift has been
added to make the rubberbox more visible.

******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
/**************** end self doc ********************************/


void 
xRubberBox (Display *dpy, Window win, XEvent event,
	int *x, int *y, int *width, int *height)
/*****************************************************************************
Track pointer with rubber box
******************************************************************************
Input:
dpy		display pointer
win		window ID
event		event of type ButtonPress

Output:
x		x of upper left hand corner of box in pixels
y		y of upper left hand corner of box in pixels
width		width of box in pixels
height		height of box in pixels
******************************************************************************
Notes:
xRubberBox assumes that event is a ButtonPress event for the 1st button;
i.e., it tracks motion of the pointer while the 1st button is down, and
it sets x, y, w, and h and returns after a ButtonRelease event for the
1st button.

Before calling xRubberBox, both ButtonRelease and Button1Motion events 
must be enabled.
******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
{
	GC gc;
	XGCValues *values=NULL;
	XEvent eventb;
	XStandardColormap scmap;
	int scr=DefaultScreen(dpy);
	int xb,yb,w,h,x1,x2,y1,y2,xorig,yorig,xold,yold;
	unsigned long background;

	/* determine typical background color */
	/* +1 added by John Stockwell 23 Jun 1993 */
	/* to shift xwigb rubberbox from light green to red */
	if (xCreateRGBDefaultMap(dpy,&scmap))
		background = (xGetFirstPixel(dpy)+xGetLastPixel(dpy) + 1)/2;
	else
		background = WhitePixel(dpy,scr);


	/* make graphics context */
	gc = XCreateGC(dpy,win,0,values);
  	XSetFunction(dpy,gc,GXxor);
  	XSetForeground(dpy,gc,BlackPixel(dpy,scr)^background);

	/* track pointer */
	xorig = event.xbutton.x;
	yorig = event.xbutton.y;
	xold = xorig;
	yold = yorig;
	x1 = xorig;
	y1 = yorig;
	w = 0;
	h = 0;
	while(h|(~h)/*True*/) {
		XNextEvent(dpy,&eventb);
		if (eventb.type==ButtonRelease) {
			xb = eventb.xbutton.x;
			yb = eventb.xbutton.y;
			break;
		} else if (eventb.type==MotionNotify) {
			xb = eventb.xmotion.x;
			yb = eventb.xmotion.y;

			/* if box is the same, continue */
			if (xb==xold && yb==yold) 
				continue;

			/* erase old box */
			x1 = (xold<xorig)?xold:xorig;
			y1 = (yold<yorig)?yold:yorig;
			x2 = (xold>xorig)?xold:xorig;
			y2 = (yold>yorig)?yold:yorig;
			w = x2-x1;
			h = y2-y1;
			XDrawRectangle(dpy,win,gc,x1,y1,w,h);

			/* draw current box */
			x1 = (xb<xorig)?xb:xorig;
			y1 = (yb<yorig)?yb:yorig;
			x2 = (xb>xorig)?xb:xorig;
			y2 = (yb>yorig)?yb:yorig;
			w = x2-x1;
			h = y2-y1;
			XDrawRectangle(dpy,win,gc,x1,y1,w,h);

			/* remember current pointer position */
			xold = xb;
			yold = yb;
		}
	}

	/* erase rubber box */
	XDrawRectangle(dpy,win,gc,x1,y1,w,h);

	/* free graphics context */
	XFreeGC(dpy,gc);

	/* set output parameters */
	*x = x1;
	*y = y1;
	*width = w;
	*height = h;
}

/* This function picks an interval and returns the screen coordinates
of the interval in xleft and xright.  Note that xleft and xright
are screen x in NORMAL style and screen y in SEISMIC style.
Caller must be conscious of this.

The procedure was derived from the xRubberBox function immediately
above.  The main difference is that instead of drawing a zoom box
from click to releas, the box is expanded to the limits of the plot.
All the plotting stuff is identical to xRubberBox except it contains
a conditional for SEISMIC versus NORMAL plot style.

Note there is an error in the comments above about xRubberBox.
This function works with release of any mouse button.  It is
not specific to button 1, and in fact at the moment button3
is used for this function.

Author:  Gary Pavlis and Peng Wang
*/
void pick_x1_interval(Display *dpy,Window win,XEvent event,
	int style, int x, int y, int width, int height,
	int *xleft, int *xright, Widget * wdgt, SeismicPlot * sptr)

{
	GC gc;
	XGCValues *values=NULL;
	XEvent eventb;
	XStandardColormap scmap;
	int scr=DefaultScreen(dpy);
	int xb,yb,w,h,x1,x2,y1,y2,xorig,yorig,xold,yold;
	unsigned long background;

	/* determine typical background color */
	/* +1 added by John Stockwell 23 Jun 1993 */
	/* to shift xwigb rubberbox from light green to red */
	if (xCreateRGBDefaultMap(dpy,&scmap))
		background = (xGetFirstPixel(dpy)+xGetLastPixel(dpy) + 1)/2;
	else
		background = WhitePixel(dpy,scr);


	/* make graphics context */
	gc = XCreateGC(dpy,win,0,values);
  	XSetFunction(dpy,gc,GXxor);
  	XSetForeground(dpy,gc,BlackPixel(dpy,scr)^background);
  	if(style==SEISMIC)
  	{

	    /* track pointer */
	    xorig = x;
	    yorig = event.xbutton.y;
	    xold = xorig;
	    yold = yorig;
	    x1 = xorig;
	    y1 = yorig;
	    w = width;
	    h = 0;
	    while(h|(~h)/*True*/) {
		XNextEvent(dpy,&eventb);
		if (eventb.type==ButtonRelease) {
			xb = eventb.xbutton.x;
			yb = eventb.xbutton.y;
			break;
		} else if (eventb.type==MotionNotify) {
			xb = eventb.xmotion.x;
			yb = eventb.xmotion.y;

			/* if box is the same, continue */
			if (xb==xold && yb==yold) 
				continue;

			/* erase old box */
			x1 = (xold<xorig)?xold:xorig;
			y1 = (yold<yorig)?yold:yorig;
			x2 = (xold>xorig)?xold:xorig;
			y2 = (yold>yorig)?yold:yorig;
			w = width;
			h = y2-y1;
			XDrawRectangle(dpy,win,gc,x,y1,w,h);

			/* draw current box */
			x1 = (xb<xorig)?xb:xorig;
			y1 = (yb<yorig)?yb:yorig;
			x2 = (xb>xorig)?xb:xorig;
			y2 = (yb>yorig)?yb:yorig;
			w = width;
			h = y2-y1;
			XDrawRectangle(dpy,win,gc,x,y1,w,h);

                        //only relation with gui here
			/*  Removed for now (12/30/2005) by GLP
                        if (wdgt != NULL && wdgt != 0)
                            adjust_scrollbar(dpy, win, *wdgt, xb, yb, sptr);
			****************************************/

			/* remember current pointer position */
			xold = x;
			yold = yb;
		}
	    }
	}
	else
	{
	    /* track pointer */
	    xorig = event.xbutton.x;
	    yorig = y;
	    xold = xorig;
	    yold = yorig;
	    x1 = xorig;
	    y1 = yorig;
	    w = 0;
	    h = height;
	    while(h|(~h)/*True*/) {
		XNextEvent(dpy,&eventb);
		if (eventb.type==ButtonRelease) {
			xb = eventb.xbutton.x;
			yb = eventb.xbutton.y;
			break;
		} else if (eventb.type==MotionNotify) {
			xb = eventb.xmotion.x;
			yb = eventb.xmotion.y;

			/* if box is the same, continue */
			if (xb==xold && yb==yold) 
				continue;

			/* erase old box */
			x1 = (xold<xorig)?xold:xorig;
			y1 = (yold<yorig)?yold:yorig;
			x2 = (xold>xorig)?xold:xorig;
			y2 = (yold>yorig)?yold:yorig;
			w = x2-x1;
			h = height;
			XDrawRectangle(dpy,win,gc,x1,y,w,h);

			/* draw current box */
			x1 = (xb<xorig)?xb:xorig;
			y1 = (yb<yorig)?yb:yorig;
			x2 = (xb>xorig)?xb:xorig;
			y2 = (yb>yorig)?yb:yorig;
			w = x2-x1;
			h = height;
			XDrawRectangle(dpy,win,gc,x1,y,w,h);

                        //added by Peng Wang, only relation with gui here
			/*  Removed for now (12/30/2005) by GLP
                        if (wdgt != NULL && wdgt != 0)
                            adjust_scrollbar(dpy, win, *wdgt, xb, yb, sptr);
			**********************************************/

			/* remember current pointer position */
			xold = xb;
			yold = y;
		}
	    }	
	}
	/* erase rubber box */
	if(style==SEISMIC)
		XDrawRectangle(dpy,win,gc,x,y1,w,h);
	else
		XDrawRectangle(dpy,win,gc,x1,y,w,h);		

	/* free graphics context */
	XFreeGC(dpy,gc);

	/* set output parameters */
	if(style==SEISMIC)
	{
		*xleft=y1;
		*xright=y1+h;
	}
	else
	{
		*xleft=x1;
		*xright=x1+w;
	}
}
/*  Formula convert screen coordinates to plot units.  
Formula is independent of plot style because in X 
the upper left corner of the screen is the origin.
See comment below about offset argument.

Arguments:
	pick - int screen value to convert 
	offset - screen coordinate of left (NORMAL) or top (SEISMIC) 
		of plot box.
	plot_width - width of plotting box in screen coordinates
	x1begb - x1 value of left (top) of plotting box
	x1endb - x1 value of right (bottom) of plotting box
*/
double screen_to_time(int pick,
		int offset,
			int plot_width,
				float x1begb,
					float x1endb)

{
	double time;
	time=static_cast<double>(x1begb)
		+ static_cast<double>(x1endb-x1begb)
		* static_cast<double>((pick-offset))
		/ static_cast<double>(plot_width);
	return(time);
}

/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
SCAXIS - compute a readable scale for use in plotting axes

scaxis		compute a readable scale for use in plotting axes

******************************************************************************
Function Prototype:
void scaxis (float x1, float x2, int *nxnum, float *dxnum, float *fxnum);

******************************************************************************
Input:
x1		first x value
x2		second x value
nxnum		desired number of numbered values

Output:
nxnum		number of numbered values
dxnum		increment between numbered values (dxnum>0.0)
fxnum		first numbered value

******************************************************************************
Notes:
scaxis attempts to honor the user-specified nxnum.  However, nxnum
will be modified if necessary for readability.  Also, fxnum and nxnum
will be adjusted to compensate for roundoff error; in particular, 
fxnum will not be less than xmin-eps, and fxnum+(nxnum-1)*dxnum 
will not be greater than xmax+eps, where eps = 0.0001*(xmax-xmin).
xmin is the minimum of x1 and x2.  xmax is the maximum of x1 and x2.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
/**************** end self doc ********************************/



void scaxis (float x1, float x2, int *nxnum, float *dxnum, float *fxnum)
/*****************************************************************************
compute a readable scale for use in plotting axes
******************************************************************************
Input:
x1		first x value
x2		second x value
nxnum		desired number of numbered values

Output:
nxnum		number of numbered values
dxnum		increment between numbered values (dxnum>0.0)
fxnum		first numbered value
******************************************************************************
Notes:
scaxis attempts to honor the user-specified nxnum.  However, nxnum
will be modified if necessary for readability.  Also, fxnum and nxnum
will be adjusted to compensate for roundoff error; in particular, 
fxnum will not be less than xmin-eps, and fxnum+(nxnum-1)*dxnum 
will not be greater than xmax+eps, where eps = 0.0001*(xmax-xmin).
xmin is the minimum of x1 and x2.  xmax is the maximum of x1 and x2.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int n,i,iloga;
	float d,f,rdint[4],eps,a,b,xmin,xmax;

	/* set readable intervals */
	rdint[0] = 1.0;  rdint[1] = 2.0;  rdint[2] = 5.0;  rdint[3] = 10.0;

	/* handle x1==x2 as a special case */
	if  (x1==x2) {
		*nxnum = 1;
		*dxnum = 1.0;
		*fxnum = x1;
		return;
	}

	/* determine minimum and maximum x */
	xmin = (x1<x2)?x1:x2;
	xmax = (x1>x2)?x1:x2;
	
	/* get desired number of numbered values */
	n = *nxnum;
	n = (2>n)?2:n;
	
	/* determine output parameters, adjusted for roundoff */
	a = (xmax-xmin)/(float)(n-1);
	iloga = (int)log10(a);
	if (a<1.0) iloga = iloga - 1;
	b = a/pow(10.0,(double)iloga);
	for (i=0; i<3 && b>=sqrt(rdint[i]*rdint[i+1]); i++);
	d = rdint[i]*static_cast<float>(pow(10.0,(double)iloga));
	f = ((int)(xmin/d))*d-d;
	eps = 0.0001*(xmax-xmin);
	while(f<(xmin-eps))
		 f = f+d;
	n = 1+(int)((xmax+eps-f)/d); 
        
	/* set output parameters before returning */
	*nxnum = n;
	*dxnum = d;
	*fxnum = f;
}
/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/* AXESBOX: $Revision: 1.2 $ ; $Date: 2006/02/07 14:37:37 $	*/

/*********************** self documentation **********************/
/*****************************************************************************
AXESBOX - Functions to draw axes in X-windows graphics

xDrawAxesBox	draw a labeled axes box
xSizeAxesBox	determine optimal origin and size for a labeled axes box

*****************************************************************************
Function Prototypes:
void xDrawAxesBox (Display *dpy, Window win,
	int x, int y, int width, int height,
	float x1beg, float x1end, float p1beg, float p1end,
	float d1num, float f1num, int n1tic, int grid1, char *label1,
	float x2beg, float x2end, float p2beg, float p2end,
	float d2num, float f2num, int n2tic, int grid2, char *label2,
	char *labelfont, char *title, char *titlefont, 
	char *axescolor, char *titlecolor, char *gridcolor,
	int style);
void xSizeAxesBox (Display *dpy, Window win, 
	char *labelfont, char *titlefont, int style,
	int *x, int *y, int *width, int *height);

*****************************************************************************
xDrawAxesBox:
Input:
dpy		display pointer
win		window
x		x coordinate of upper left corner of box
y		y coordinate of upper left corner of box
width		width of box
height		height of box
x1beg		axis value at beginning of axis 1
x1end		axis value at end of axis 1
p1beg		pad value at beginning of axis 1
p1end		pad value at end of axis 1
d1num		numbered tic increment for axis 1 (0.0 for automatic)
f1num		first numbered tic for axis 1
n1tic		number of tics per numbered tic for axis 1
grid1		grid code for axis 1:  NONE, DOT, DASH, or SOLID
label1		label for axis 1
x2beg		axis value at beginning of axis 2
x2end		axis value at end of axis 2
p2beg		pad value at beginning of axis 2
p2end		pad value at end of axis 2
d2num		numbered tic increment for axis 2 (0.0 for automatic)
f2num		first numbered tic for axis 2
n2tic		number of tics per numbered tic for axis 2
grid2		grid code for axis 2:  NONE, DOT, DASH, or SOLID
label2		label for axis 2
labelfont	name of font to use for axes labels
title		axes box title
titlefont	name of font to use for title
axescolor	name of color to use for axes
titlecolor	name of color to use for title
gridcolor	name of color to use for grid
int style	NORMAL (axis 1 on bottom, axis 2 on left)
		SEISMIC (axis 1 on left, axis 2 on top)

******************************************************************************
xSizeAxesBox:
Input:
dpy		display pointer
win		window
labelfont	name of font to use for axes labels
titlefont	name of font to use for title
int style	NORMAL (axis 1 on bottom, axis 2 on left)
		SEISMIC (axis 1 on left, axis 2 on top)

Output:
x		x coordinate of upper left corner of box
y		y coordinate of upper left corner of box
width		width of box
height		height of box
******************************************************************************
{
	XFontStruct *fa,*ft;
******************************************************************************
Notes:
xDrawAxesBox:
will determine the numbered tic incremenet and first
numbered tic automatically, if the specified increment is zero.

Pad values must be specified in the same units as the corresponding
axes values.  These pads are useful when the contents of the axes box
requires more space than implied by the axes values.  For example,
the first and last seismic wiggle traces plotted inside an axes box
will typically extend beyond the axes values corresponding to the
first and last traces.  However, all tics will lie within the limits
specified in the axes values (x1beg, x1end, x2beg, x2end).

xSizeAxesBox:
is intended to be used prior to xDrawAxesBox.

An "optimal" axes box is one that more or less fills the window, 
with little wasted space around the edges of the window.

******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
/**************** end self doc ********************************/


void
xDrawAxesBox (Display *dpy, Window win,
	int x, int y, int width, int height,
	float x1beg, float x1end, float p1beg, float p1end,
	float d1num, float f1num, int n1tic, int grid1, char *label1,
	float x2beg, float x2end, float p2beg, float p2end,
	float d2num, float f2num, int n2tic, int grid2, char *label2,
	char *labelfont, char *title, char *titlefont, 
	char *axescolor, char *titlecolor, char *gridcolor,
	int style)
/*****************************************************************************
draw a labeled axes box
******************************************************************************
Input:
dpy		display pointer
win		window
x		x coordinate of upper left corner of box
y		y coordinate of upper left corner of box
width		width of box
height		height of box
x1beg		axis value at beginning of axis 1
x1end		axis value at end of axis 1
p1beg		pad value at beginning of axis 1
p1end		pad value at end of axis 1
d1num		numbered tic increment for axis 1 (0.0 for automatic)
f1num		first numbered tic for axis 1
n1tic		number of tics per numbered tic for axis 1
grid1		grid code for axis 1:  NONE, DOT, DASH, or SOLID
label1		label for axis 1
x2beg		axis value at beginning of axis 2
x2end		axis value at end of axis 2
p2beg		pad value at beginning of axis 2
p2end		pad value at end of axis 2
d2num		numbered tic increment for axis 2 (0.0 for automatic)
f2num		first numbered tic for axis 2
n2tic		number of tics per numbered tic for axis 2
grid2		grid code for axis 2:  NONE, DOT, DASH, or SOLID
label2		label for axis 2
labelfont	name of font to use for axes labels
title		axes box title
titlefont	name of font to use for title
axescolor	name of color to use for axes
titlecolor	name of color to use for title
gridcolor	name of color to use for grid
int style	NORMAL (axis 1 on bottom, axis 2 on left)
		SEISMIC (axis 1 on left, axis 2 on top)
******************************************************************************
Notes:
xDrawAxesBox will determine the numbered tic incremenet and first
numbered tic automatically, if the specified increment is zero.

Pad values must be specified in the same units as the corresponding
axes values.  These pads are useful when the contents of the axes box
requires more space than implied by the axes values.  For example,
the first and last seismic wiggle traces plotted inside an axes box
will typically extend beyond the axes values corresponding to the
first and last traces.  However, all tics will lie within the limits
specified in the axes values (x1beg, x1end, x2beg, x2end).
******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
{
	GC gca,gct,gcg;
	XGCValues *values=NULL;
	XColor scolor,ecolor;
	XFontStruct *fa,*ft;
	XWindowAttributes wa;
	Colormap cmap;
	int labelca,labelcd,labelch,labelcw,titleca,
		ntic,xa,ya,tw,ticsize,ticb,numb,labelb,lstr,grided,grid,
		n1num,n2num;
	float dnum,fnum,dtic,amin,amax,base,scale,anum,atic,azero;
	char str[256],dash[2],*label;

	/* create graphics contexts */
	gca = XCreateGC(dpy,win,0,values);
	gct = XCreateGC(dpy,win,0,values);
	gcg = XCreateGC(dpy,win,0,values);


	/* get and set fonts and determine character dimensions */
	fa = XLoadQueryFont(dpy,labelfont);
	if (fa==NULL) fa = XLoadQueryFont(dpy,"fixed");
	if (fa==NULL) {
		fprintf(stderr,"Cannot load/query labelfont=%s\n",labelfont);
		exit(-1);
	}
	XSetFont(dpy,gca,fa->fid);
	labelca = fa->max_bounds.ascent;
	labelcd = fa->max_bounds.descent;
	labelch = fa->max_bounds.ascent+fa->max_bounds.descent;
	labelcw = fa->max_bounds.lbearing+fa->max_bounds.rbearing;
	ft = XLoadQueryFont(dpy,titlefont);
	if (ft==NULL) ft = XLoadQueryFont(dpy,"fixed");
	if (ft==NULL) {
		fprintf(stderr,"Cannot load/query titlefont=%s\n",titlefont);
		exit(-1);
	}
	XSetFont(dpy,gct,ft->fid);
	titleca = ft->max_bounds.ascent;

	/* determine window's current colormap */
	XGetWindowAttributes(dpy,win,&wa);
	cmap = wa.colormap;

	/* get and set colors */
	if (XAllocNamedColor(dpy,cmap,axescolor,&scolor,&ecolor))
		XSetForeground(dpy,gca,ecolor.pixel);
	else
		XSetForeground(dpy,gca,1L);
	if (XAllocNamedColor(dpy,cmap,titlecolor,&scolor,&ecolor))
		XSetForeground(dpy,gct,ecolor.pixel);
	else
		XSetForeground(dpy,gct,1L);
	if (XAllocNamedColor(dpy,cmap,gridcolor,&scolor,&ecolor))
		XSetForeground(dpy,gcg,ecolor.pixel);
	else
		XSetForeground(dpy,gcg,1L);

	/* determine tic size */
	ticsize = labelcw;

	/* determine numbered tic intervals */
	if (d1num==0.0) {
		n1num = (style==NORMAL ? width : height)/(8*labelcw);
		scaxis(x1beg,x1end,&n1num,&d1num,&f1num);
	}
	if (d2num==0.0) {
		n2num = (style==NORMAL ? height : width)/(8*labelcw);
		scaxis(x2beg,x2end,&n2num,&d2num,&f2num);
	}

	/* draw horizontal axis */
	if (style==NORMAL) {
		amin = (x1beg<x1end)?x1beg:x1end;
		amax = (x1beg>x1end)?x1beg:x1end;
		dnum = d1num;  fnum = f1num;  ntic = n1tic;
		scale = width/(x1end+p1end-x1beg-p1beg);
		base = x-scale*(x1beg+p1beg);
		ya = y+height;
		ticb = ticsize;
		numb = ticb+labelca;
		labelb = numb+labelch;
		grid = grid1;
		label = label1;
	} else {
		amin = (x2beg<x2end)?x2beg:x2end;
		amax = (x2beg>x2end)?x2beg:x2end;
		dnum = d2num;  fnum = f2num;  ntic = n2tic;
		scale = width/(x2end+p2end-x2beg-p2beg);
		base = x-scale*(x2beg+p2beg);
		ya = y;
		ticb = -ticsize;
		numb = ticb-labelcd;
		labelb = numb-labelch;
		grid = grid2;
		label = label2;
	}
	if (grid==SOLID)
		grided = True;
	else if (grid==DASH) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 8;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else if (grid==DOT) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 1;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else
		grided = False;
	azero = 0.0001*(amax-amin);
	for (anum=fnum; anum<=amax; anum+=dnum) {
		if (anum<amin) continue;
		xa = base+scale*anum;
		if (grided) XDrawLine(dpy,win,gcg,xa,y,xa,y+height);
		XDrawLine(dpy,win,gca,xa,ya,xa,ya+ticb);
		if (anum>-azero && anum<azero)
			sprintf(str,"%1.5g",0.0);
		else
			sprintf(str,"%1.5g",anum);
		lstr = (int) strlen(str);
		tw = XTextWidth(fa,str,lstr);
		XDrawString(dpy,win,gca,xa-tw/2,ya+numb,str,lstr);
	}
	dtic = dnum/ntic;
	for (atic=fnum-ntic*dtic-dtic; atic<=amax; atic+=dtic) {
		if (atic<amin) continue;
		xa = base+scale*atic;
		XDrawLine(dpy,win,gca,xa,ya,xa,ya+ticb/2);
	}
	lstr = (int) strlen(label);
	tw = XTextWidth(fa,label,lstr);
	XDrawString(dpy,win,gca,x+width-tw,ya+labelb,label,lstr);

	/* draw vertical axis */
	if (style==NORMAL) {
		amin = (x2beg<x2end)?x2beg:x2end;
		amax = (x2beg>x2end)?x2beg:x2end;
		dnum = d2num;  fnum = f2num;  ntic = n2tic;
		scale = -height/(x2end+p2end-x2beg-p2beg);
		base = y+height-scale*(x2beg+p2beg);
		grid = grid2;
		label = label2;
	} else {
		amin = (x1beg<x1end)?x1beg:x1end;
		amax = (x1beg>x1end)?x1beg:x1end;
		dnum = d1num;  fnum = f1num;  ntic = n1tic;
		scale = height/(x1end+p1end-x1beg-p1beg);
		base = y-scale*(x1beg+p1beg);
		grid = grid1;
		label = label1;
	}
	xa = x;
	ticb = -ticsize;
	numb = ticb-ticsize/4;
	if (grid==SOLID)
		grided = True;
	else if (grid==DASH) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 8;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else if (grid==DOT) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 1;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else
		grided = False;
	azero = 0.0001*(amax-amin);
	for (anum=fnum; anum<=amax; anum+=dnum) {
		if (anum<amin) continue;
		ya = base+scale*anum;
		if (grided) XDrawLine(dpy,win,gcg,x,ya,x+width,ya);
		XDrawLine(dpy,win,gca,xa,ya,xa+ticb,ya);
		if (anum>-azero && anum<azero)
			sprintf(str,"%1.5g",0.0);
		else
			sprintf(str,"%1.5g",anum);
		lstr = (int) strlen(str);
		tw = XTextWidth(fa,str,lstr);
		XDrawString(dpy,win,gca,xa+numb-tw,ya+labelca/4,str,lstr);
	}
	dtic = dnum/ntic;
	for (atic=fnum-ntic*dtic-dtic; atic<=amax; atic+=dtic) {
		if (atic<amin) continue;
		ya = base+scale*atic;
		XDrawLine(dpy,win,gca,xa,ya,xa+ticb/2,ya);
	}
	lstr = (int) strlen(label);
	if (style==NORMAL)
		XDrawString(dpy,win,gca,
			x+ticb-9*labelcw,
			y+labelca/4-labelch,label,lstr);
	else
		XDrawString(dpy,win,gca,
			x+ticb-9*labelcw,
			y+height+labelca/4+labelch,label,lstr);
	
	/* draw title */
	lstr = (int) strlen(title);
	tw = XTextWidth(ft,title,lstr);
	if (style==NORMAL)
		XDrawString(dpy,win,gct,
			x+width/2-tw/2,
			y+labelca/4-labelch-labelch,title,lstr);
	else
		XDrawString(dpy,win,gct,
			x+width/2-tw/2,
			y+height+labelca/4+labelch+titleca,title,lstr);

	/* draw axes box */
	XDrawRectangle(dpy,win,gca,x,y,width,height);

	/* free resources before returning */
	XFreeGC(dpy,gca);
	XFreeGC(dpy,gct);
	XFreeGC(dpy,gcg);
	XFreeFont(dpy,fa);
	XFreeFont(dpy,ft);
}

void
xSizeAxesBox (Display *dpy, Window win, 
	char *labelfont, char *titlefont, int style,
	int *x, int *y, int *width, int *height)
/*****************************************************************************
determine optimal origin and size for a labeled axes box
******************************************************************************
Input:
dpy		display pointer
win		window
labelfont	name of font to use for axes labels
titlefont	name of font to use for title
int style	NORMAL (axis 1 on bottom, axis 2 on left)
		SEISMIC (axis 1 on left, axis 2 on top)

Output:
x		x coordinate of upper left corner of box
y		y coordinate of upper left corner of box
width		width of box
height		height of box
******************************************************************************
Notes:
xSizeAxesBox is intended to be used prior to xDrawAxesBox.

An "optimal" axes box is one that more or less fills the window, 
with little wasted space around the edges of the window.
******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
{
	XFontStruct *fa,*ft;
	XWindowAttributes attr;
	int labelch,labelcw,titlech,bl,bt,br,bb;

	/* get fonts and determine character dimensions */
	fa = XLoadQueryFont(dpy,labelfont);
	if (fa==NULL) fa = XLoadQueryFont(dpy,"fixed");
	if (fa==NULL) {
		fprintf(stderr,"Cannot load/query labelfont=%s\n",labelfont);
		exit(-1);
	}
	labelch = fa->max_bounds.ascent+fa->max_bounds.descent;
	labelcw = fa->max_bounds.lbearing+fa->max_bounds.rbearing;
	ft = XLoadQueryFont(dpy,titlefont);
	if (ft==NULL) ft = XLoadQueryFont(dpy,"fixed");
	if (ft==NULL) {
		fprintf(stderr,"Cannot load/query titlefont=%s\n",titlefont);
		exit(-1);
	}
	titlech = ft->max_bounds.ascent+ft->max_bounds.descent;

	/* determine axes box origin and size */
	XGetWindowAttributes(dpy,win,&attr);
	bl = 10*labelcw;
	br = attr.width-5*labelcw;
	while (br<=bl) {
		br += labelcw;
		bl -= labelcw;
	}
	if (bl<0) bl = 0;
	if (br>attr.width) br = attr.width;
	if (style==NORMAL) {
		bt = labelch+labelch/2+titlech;
		bb = attr.height-3*labelch;
	} else {
		bt = 3*labelch;
		bb = attr.height-labelch-labelch/2-titlech;
	}
	while (bb<=bt) {
		bb += labelch;
		bt -= labelch;
	}
	if (bt<0) bt = 0;
	if (bb>attr.height) bb = attr.height;
	
	*x = bl;
	*y = bt;
	*width = br-bl;
	*height = bb-bt;

	XFreeFont(dpy,fa);
	XFreeFont(dpy,ft);
}

//
// Helper for SetParameters method below
//
int parse_grid_style(string keyword)
{
	if(keyword=="none")
		return NONE;
	else if(keyword=="dash")
		return DASH;
	else if(keyword=="solid")
		return SOLID;
	else
		return NONE;
}
void SeismicPlot::SetParameters(Metadata md)


{
	// These are temporaries needed for type conversion
	string s;
	bool bval;

	try {
		verbose=md.get_bool("verbose");
		bval=md.get_bool("WiggleTrace");
		if(bval)
		{
			wt=md.get_int("SUVariableArea_grey_value");
			if(wt<0) wt=0;
			if(wt>5) wt=5;
		}
		else
			wt=0;
		bval=md.get_bool("VariableArea");
		if(bval)
		{
			va=md.get_int("SUVariableArea_grey_value");
			if(va<0) va=0;
			if(va>5) va=5;
		}
		else
			va=0;
		
		clip_data=md.get_bool("clip_data");
		if(clip_data)
			perc=static_cast<float>(md.get_double("clip_percent"));
		bval=md.get_bool("clip_wiggle_traces");
		if(bval)
			wigclip=1;
		else
			wigclip=0;
		// wiggle excursion in traces for clip
		// cryptic name, but no obvious longer name that is self-documenting
		xcur=static_cast<float>(md.get_double("xcur"));
		
		xbox=md.get_int("xbox");
		ybox=md.get_int("ybox");
		wbox=md.get_int("wbox");
		hbox=md.get_int("hbox");
		s=md.get_string("style");
		d1num=static_cast<float>(md.get_double("d1num"));
		d2num=static_cast<float>(md.get_double("d2num"));
		f1num=static_cast<float>(md.get_double("f1num"));
		f2num=static_cast<float>(md.get_double("f2num"));
		n1tic=md.get_int("n1tic");
		n2tic=md.get_int("n2tic");
		s=md.get_string("label1");
		label1=strdup(s.c_str());
		s=md.get_string("blabel2");
		label2=strdup(s.c_str());
		s=md.get_string("title");
		title=strdup(s.c_str());
		s=md.get_string("windowtitle");
		windowtitle=strdup(s.c_str());
		s=md.get_string("labelfont");
		labelfont=strdup(s.c_str());
		s=md.get_string("titlefont");
		titlefont=strdup(s.c_str());
		s=md.get_string("style");
		if(s=="seismic" || s=="SEISMIC")
			style=SEISMIC;
		else
			style=NORMAL;
		s=md.get_string("time_axis_grid_type");
		grid1=parse_grid_style(s);
		s=md.get_string("trace_axis_grid_type");
		grid2=parse_grid_style(s);

		s=md.get_string("labelcolor");
		labelcolor=strdup(s.c_str());
		s=md.get_string("titlecolor");
		titlecolor=strdup(s.c_str());
		s=md.get_string("gridcolor");
		gridcolor=strdup(s.c_str());
		grid1=md.get_int("grid1");
		grid2=md.get_int("grid2");
		labelsize=static_cast<float>(md.get_double("labelsize"));
		titlesize=static_cast<float>(md.get_double("titlesize"));
		trace_spacing=md.get_double("trace_spacing");
		first_trace_offset=md.get_double("first_trace_offset");		
		time_scaling=md.get_string("time_scaling");
		x1beg=static_cast<float>(md.get_double("x1beg"));
		x1end=static_cast<float>(md.get_double("x1end"));
		trace_axis_scaling=md.get_string("trace_axis_scaling");
		x2beg=static_cast<float>(md.get_double("x2beg"));
		x2end=static_cast<float>(md.get_double("x2end"));
		bval=md.get_bool("interpolate");
		if(bval)
			interp=1;
		else
			interp=0;
	
		default_curve_color=md.get_string("default_curve_color");
		plotfile=md.get_string("plot_file_name");
		// variable trace spacing option
		use_variable_trace_spacing=md.get_bool("use_variable_trace_spacing");
	}
	catch (MetadataGetError mde)
	{
		cerr << "SeismicPlot constructor:  input parameter error"<<endl;
		mde.log_error();
		throw SeisppError("SeismicPlot constructor failure");
	}
}
SeismicPlot::~SeismicPlot()
{

	/* close connection to X server */
	XCloseDisplay(dpy);
	// free char * variables 
	free(labelfont);
	free(titlefont);
	free(labelcolor);
	free(titlecolor);
	free(label1);
	free(label2);
	free(title);
	free(windowtitle);
	free(gridcolor);
	if(x2!=NULL) delete [] x2;
	if(z!=NULL) delete [] z;
}


void SeismicPlot::X_init()
{
        XSizeHints size_hints;
        XWMHints wm_hints;

	/**************************************************************
	// Peng:  Previous version I'm pretty sure got in trouble
	// from violating the rule that in OOP design of a resource
	// management class.  The SeismicPlot is a resource management
	// object with the resource being the display and window
	// combination.  To use this class correctly in your gui
	// you must never mess the display or Window variables.  
	// For the gui with Display derived from a widget, 
	// this resource must be passed into this object from 
	// outside.  You'll need a method to do this and I suggest
	// you use the display_name here to pass it.  You could
	// change this block to trigger different behaviour to handle
	// getting Display from a widget as you do in gui.cc
	// with the call to XtDisplay.  The Screen and Window 
	// variables should be handled the same way.  For now
	// I'm reverting to a simple window version that I have
	// some hope of understanding.  glp
	**************************************************************/
	try {
		string display_name=this->get_string("display_name");
		dpy=XOpenDisplay(display_name.c_str());
	} catch (MetadataError mderr)
	{
		// Silently try this if display_name is not defined
		dpy=XOpenDisplay(NULL);
	}
	if(dpy==NULL)
		throw SeisppError("SeismicPlot::X_init:  XOpenDisplay failed");

	// simplest startup solution.  Needs to be generalized
        scr = DefaultScreen(dpy);
	win=RootWindow(dpy,scr);
        unsigned long black, white;
        black = BlackPixel(dpy,scr);
        white = WhitePixel(dpy,scr);
	// Simple window launch with SU procedure
	win=xNewWindow(dpy,xbox,ybox,wbox,hbox,
		(int) black,(int) white,windowtitle);

        /* set endian for display.  Original a bit more complex and involved
        an apparent default passed as a parameter.  This allowed bypass of
        BitmapBitOrder call below.*/
        if(BitmapBitOrder(dpy)==LSBFirst)
                endian=0;
        else if(BitmapBitOrder(dpy)==MSBFirst)
                endian=1;

        /* make GC for image */
        gci = XCreateGC(dpy,win,0,NULL);

        /* make sure foreground/background are black/white */
        XSetForeground(dpy,gci,black);
        XSetBackground(dpy,gci,white);

        /* set normal event mask */
        XSelectInput(dpy,win,
                StructureNotifyMask |
                ExposureMask |
                KeyPressMask |
                PointerMotionMask |
                ButtonPressMask |
                ButtonReleaseMask |
                Button1MotionMask |
                Button2MotionMask);

        /* map window */
        XMapWindow(dpy,win);

        /* clear the window */
        XClearWindow(dpy,win);
}

SeismicPlot::SeismicPlot(TimeSeriesEnsemble& tse, Metadata md)  : TimeSeriesEnsemble(tse)
{
	int i;
	// first thing to do is delete data marked dead from the ensemble.  
	// for this implementation we will be dogmatic about this and demand it be
	// so. There probably should be an option to zero dead traces and leave them
	// in the ensemble, but I prefer this approach for now.
	vector<TimeSeries>::iterator memiter;
	for(memiter=member.begin();memiter!=member.end();++memiter)
		if(!memiter->live) member.erase(memiter);
	nmember=member.size();
	if(nmember<=0) throw SeisppError("Seisplot constructor: no data to plot\n");
	this->SetParameters(md);
	for(i=0;i<nmember;++i)
		curvecolor.push_back(default_curve_color);

	x2=new float[nmember];
	if(use_variable_trace_spacing)
	{
		try {
			for(i=0;i<nmember;++i)
				x2[i]=static_cast<double>(member[i].get_double(trace_axis_attribute));
		} catch (MetadataGetError mde)
		{
			mde.log_error();
			cerr << "Reverting to equal space tracing" << endl;
			for(i=0;i<nmember;++i) x2[i]=static_cast<float>(i+1);
		}
	}
	else
	{
		 for(i=0;i<nmember;++i) x2[i]=static_cast<float>(i+1);
	}
	// explicit initialization
	z=NULL;

        //initialize command state to NoCmd
        scmd=NoCmd;
	// This is standard OOP style that construction is 
	// acquisition of a resource.  In this case, the resource
	// is the X server connect.  The interface hides the details
	// and insulates this API from an X implementation
	this->X_init();

}


SeismicPlot::SeismicPlot(ThreeComponentEnsemble& tce, int comp, Metadata md) 
{
	int i;
	nmember=tce.member.size();
	if(nmember<=0) throw SeisppError("Seisplot constructor: no data to plot\n");
	this->SetParameters(md);
	member.reserve(nmember);
	// this loads component comp as the member to be plotted.
	for(int i=0;i<nmember;++i)
	{
		TimeSeries *c;
		if(tce.member[i].live)
		{
			c=ExtractComponent(tce.member[i],comp);
			member.push_back(*c);
			delete c;
			curvecolor.push_back(default_curve_color);
		}
	}
	nmember=member.size();
	x2=new float[nmember];
	if(use_variable_trace_spacing)
	{
		try {
			for(i=0;i<nmember;++i)
			{
				x2[i]=static_cast<double>(member[i].get_double(trace_axis_attribute));
			}
		} catch (MetadataGetError mde)
		{
			mde.log_error();
			cerr << "Reverting to equal space tracing" << endl;
			for(i=0;i<nmember;++i) x2[i]=static_cast<float>(i);
		}
	}
	else
	{
		 for(i=0;i<nmember;++i) x2[i]=static_cast<float>(i);
	}
	// explicit initialization
	z=NULL;

	//initialize command state to NoCmd
	scmd=NoCmd;
	this->X_init();
}

void SeismicPlot::LoadNewData(ThreeComponentEnsemble& tse,int comp) 
{
	nmember=tse.member.size();
	if(nmember<=0) throw SeisppError("Seisplot::LoadNewData: no data to plot\n");
	member.reserve(nmember);
	for(int i=0;i<nmember;++i)
	{
		TimeSeries *c;
		if(tse.member[i].live)
		{
			c=ExtractComponent(tse.member[i],comp);
			member[i]=(*c);
			delete c;
		}
	}
}

void SeismicPlot::LoadNewData(TimeSeriesEnsemble& tse)
{
	nmember=tse.member.size();
	if(nmember<=0) throw SeisppError("Seisplot::LoadNewData: no data to plot\n");
	member.resize(nmember);
	for(int i=0;i<nmember;++i)
		if(tse.member[i].live) member[i]=tse.member[i];
}
/* Internal function need to rasterize traces.  Loads a float fortran-like matrix of
data from an ensemble.   We pass the matrix to avoid need for exception handling
in creation of z, which is not trivial because it interacts with data concepts. 
In particular, one could easily request an absurd plot window that would make the
buffer huge or zero length.

The algorithm used here is most efficient if the window defined by z is less than 
or equal to the time span of the data because it hits every sample in z.


args:
	data - ensemble of data 
	z - float buffer of size n1*n2
	n1 - first dimension used to build z as n1xn2 fortran-like matrix
	n2 - second dimension
	t0 - start time to use for each column in z. 
*/
void load_z_matrix(vector<TimeSeries> data,
		float *z,
			int n1,
				int n2,
					double t0,
						double dt)
{
	int i,j,iz,id;
	// could initialize z to zero, but algorithm used hits every
	// sample so this isn't necessary.  The first loop assumes n2=data.size()
	for(j=0,iz=0;j<n2;++j)
	{
		double t;
		// avoid main loop below if this trace has no data in this time window
		// Test for dead data is redundant for now, but better safe than sorry.
		if( (data[j].endtime() < t0)
			|| !data[j].live )
		{
			for(i=0;i<n1;++i,++iz) z[iz]=0.0;
		}
		else
		{
			for(i=0;i<n1;++i,++iz)
			{
				t=t0+dt*i;
				if(data[j].is_gap(t))
					z[iz]=0.0;
				else if( (data[j].t0>t) && (fabs(t-data[j].t0)>dt)) 
					z[iz]=0.0;
				else
				{
					id=data[j].sample_number(t);
					if( (id>=0) && (id<data[j].ns) ) {
					     if (data[j].s.size() > 0) {   
						z[iz]=static_cast<float>(data[j].s[i]);
					     } else {
                				char message[256];
                				sprintf(message,"SeismicPlot::load_z_matrix():  "
							"time series index %d is empty\n",j);
                				throw SeisppError(string(message));
					     }
					}
				}
cerr << "DEBUG:  i, j, z[iz]: "<<i<<", "<<j<<", "<<z[iz]<<endl;
			}
		}
	}
}
void SeismicPlot::draw() 
{
	int i,j;
	int i1,i2;
	// First compute range of x1 and x2
	float x1min,x1max,x2min,x2max;

	//set mvefac=1.
	mvefac=1.0;

	x1min=member[0].t0;
	x1max=member[0].endtime();
	for(i=0;i<nmember;++i)
	{
		x1min=MIN(member[i].t0,x1min);
		x1max=MAX(member[i].endtime(),x1max);
	}

	if(use_variable_trace_spacing)
	{
		x2min=x2[0];
		x2max=x2[0];
		for(i=1;i<nmember;++i)
		{
			x2min=min(x2min,x2[i]);
			x2max=max(x2max,x2[i]);
		}
	}
	else
	{
		// this is for equal spacing.  then we know
		// data are ordered 1 to nmember
		x2min=x2[0];
		x2max=x2[nmember-1];
	}
  
	//
	// Set x1beg and x1end for auto scaling.  Not else currently
	// because in that case assume x1beg and x1end are set in construction
	//
	if(time_scaling=="auto")
	{
		x1begb=x1min;
		x1endb=x1max;
		x1beg=x1min;
		x1end=x1max;
	}
	else
	{
		x1begb=x1beg;
		x1endb=x1end;
	}
	if(trace_axis_scaling=="auto")
	{
		x2begb=x2min;
		x2endb=x2max;
		x2beg=x2min;
		x2end=x2max;
	}
	else
	{
		x2begb=x2beg;
		x2endb=x2end;
	}
	f1=x1beg;
	x2endb = static_cast<float>(x2end);

					
	/* determine good size for axes box */
	xSizeAxesBox(dpy,win,
		labelfont,titlefont,style,
		&x,&y,&width,&height);
	imageOutOfDate = 1;

	int width1=width;
	int height1=height;
	// Determine size.  Depends on setting of time_scaling.  If auto get it fromt he
	// data, otherwise use the internal parameters.
	//
	d1=static_cast<float>(member[0].dt);
	//
	// Get the number of samples from x1beg and x1end.  
	// With above logic this should work if time_scaling is auto or manual
	//
	n1=nint((x1endb-x1begb)/d1);
        n2=nmember;

	// Sanity checks on n1
	if(n1<2)
	{
		char message[256];
		sprintf(message,"SeismicPlot::refresh():  Plot range error.\n%d samples in requested time range of %lf to %lf\n",
				n1,x1begb,x1endb);
		throw SeisppError(string(message));
	}
	else if(n1>MAXPLOTSAMPLES)
	{
		char message[512];
		sprintf(message,"SeismicPlot::refresh(): Plot range request too large\nTime range of %lf to %lf requires %n samples.  Sanity check overrides request.  Check input parameters\n",
			x1begb,x1endb,n1);
		throw SeisppError(string(message));
	}
	// d2 and f2 are defined in this object by less obscure names
	d2=static_cast<float>(this->trace_spacing);
	f2=static_cast<float>(this->first_trace_offset);
	// This buffer is needed in plot loop below.  We construct a fortran-type matrix
	// of floats with the data aligned by member[i].t0 values.  Traces are in columns
	// of the matrix.  Depend on new throwing an exception if alloc fails.
	//
	//float *z=new float[n1*n2];
	if(z!=NULL)
		delete [] z;
	z=new float[n1*n2];
	// internal function used here
	int iz,nz=n1*n2;
	load_z_matrix(this->member,z,n1,n2,x1begb,d1);

	// handle clip stuff.
	// Modified from xwigb to set clip levels
	// Percentage determines percent of samples to be left unclipped.
	// This is as in xwigb.  Difference here is we use STL vector
	// and standard STL nth_element algorithm
	// instead of SU's internal quick sort routine.
	vector<float> temp;
	temp.reserve(nz);
	if(!clip_data) perc=100.0;
        for (iz=0; iz<nz; iz++)temp.push_back(fabs(z[iz]));
	vector<float>::iterator iziter;
        iz = (nz*perc/100.0);
        if (iz<0) iz = 0;
        if (iz>nz-1) iz = nz-1;
	iziter=temp.begin()+iz;
	nth_element(temp.begin(),iziter,temp.end());
        clip = *iziter;

	/* main event loop */
	p2beg=0.0,p2end=0.0;
	XEvent event;
	if(image==NULL)
		image=NULL;
	else {
		// Not sure why, but this generates a seg
		// fault.  Turn it off now recognizing it will
		// cause a memory leak that needs to eventually
		// be plugged:  glp, 1/5/2006
		//XDestroyImage(image);
		image=NULL;
	}


	// This seens to force an initialization on first pass of these variables
	int winwidth=-1;  
	int winheight=-1;
	showloc=0;

	int num_iter=0;
	imageOutOfDate = 1;
	render_image();
}


void SeismicPlot::render_image()
{

	XEvent event;

	try {
		//dpy=dpy_in; win=win_in;
// This is changed from XCheckTypedEvent, because otherwise, it would clear the menubar and scrollbar 
// events as well if we only specify the display. Obviously, menubar and scrollbars are on the same 
// display, but different windows.
		while (XCheckTypedWindowEvent(dpy,win,Expose,&event));

                /* if necessary, make new image */
                if (imageOutOfDate) {

                    if (image!=NULL) {
                        XDestroyImage(image);
                    }
                    image = newBitmap(dpy,width,height,
                            n1,d1,f1,n2,x2,z,
                            x1begb,x1endb,x2begb,x2endb,
                            xcur,clip,wt,va,
                            &p2beg,&p2end,endian,interp,
                            wigclip,style);

                    imageOutOfDate = 0;
                 }

                 /* draw image (before axes so grid lines visible) */
                 XPutImage(dpy,win,gci,image,0,0,x,y,
                        image->width,image->height);

                 /* draw axes on top of image */
                 xDrawAxesBox(dpy,win,
                          x,y,width,height,
                          x1begb,x1endb,0.0,0.0,
                          d1num,f1num,n1tic,grid1,label1,
                          x2begb,x2endb,p2beg,p2end,
                          d2num,f2num,n2tic,grid2,label2,
                          labelfont,title,titlefont,
                          labelcolor,titlecolor,gridcolor,
                          style);
	} catch (SeisppError) {
	    cerr << "Image Rendering Error in render_image"<<endl;
	    exit(-2);
	}
}

void SeismicPlot::process_commands(SplotCommandType cmd, stringstream * ss)
{
    scmd=cmd;
}
	// PENG:  NOTE FROM GLP:  This is dangerous and should
	// be changed.  dpy and win are acquired resources associated
	// with this object.  Once acquired these are items that should
	// not change should they?  I'm showing how I changed this below.

/* Peng Wang's version:
SeismicPick SeismicPlot::process_callbacks(Display * dpy_in, Window win_in, SplotEventType event_in, XEvent * event_detail, stringstream * ss, Widget * w)
**************************************/
// GLP replacement.
SeismicPick SeismicPlot::process_callbacks(SplotEventType event_in,
	XEvent * event_detail, stringstream * ss)
{
    int xb,yb,wb,hb;
    XEvent event=*event_detail;
    TimeWindow temp_tw;


	// Peng had this line, removed by glp
    //dpy=dpy_in; win=win_in;

// Find the Widget w symbol is difficult.  I think Widget wdgt 
// should be a private data member of the SeismicPlot object.
// The compiler should catch any errors so below w becomes &wdgt.
    try {

    switch(event_in) {
	case B1Down:

            /* track pointer and get new box */
            xRubberBox(dpy,win,event,&xb,&yb,&wb,&hb);

            /* if new box has tiny width or height */
            if (wb<4 || hb<4) {

            /* reset box to initial values */
            x1begb = x1beg;
            x1endb = x1end;
            x2begb = x2beg;
            x2endb = x2end;

            /* else, if new box has non-zero width */
            /* if new box has zero width or height */
            } else {

                /* calculate new box parameters */
                zoomBox(x,y,width,height,
                        xb,yb,wb,hb,
                        x2begb,x2endb,
                        x1begb,x1endb,
                        &x2begb,&x2endb,
                        &x1begb,&x1endb,
                        style);
            }

            /* clear area and force an expose event */
            XClearArea(dpy,win,0,0,0,0,True);

            /* note that image is out of date */
            imageOutOfDate = 1;
	    break;
	case B2Down:
	    if (scmd==SinglePick || scmd==MultiplePicks) {

                SeismicPick ppick=MouseLocate(dpy,win,
                                  event,style,
                                  x,y,width,height,x1begb,x1endb,
                                  x2begb,x2endb,p2beg,p2end);
                ppick.trace_number=
                         this->get_trace_number(ppick.x2);
                int isample=this->member[ppick
                            .trace_number]
                            .sample_number((double)ppick.x1);
                ppick.amplitude=this->member[ppick.trace_number].s[isample];
		if (scmd == SinglePick) scmd=NoCmd;
                return(ppick);

	    } else {
                showloc = 1;
                xMouseLoc(dpy,win,event,style,showloc,
                      x,y,width,height,x1begb,x1endb,
                      x2begb,x2endb,p2beg,p2end);
	    }
	    break;
 	case B3Down:
            int ixb,ixe;
            int wnow;
            pick_x1_interval(dpy,win,event,style,
                     x,y,width,height,
                     &ixb,&ixe, &wdgt, this);
            // Not correct, quick and dirty for testing
            if(style==SEISMIC)
            {
                temp_tw.start=screen_to_time(ixb,
                      y,height,x1begb,x1endb);
                temp_tw.end=screen_to_time(ixe,
                      y,height,x1begb,x1endb);
            }
            else
            {
                temp_tw.start=screen_to_time(ixb,
                         x,width,x1begb,x1endb);
                temp_tw.end=screen_to_time(ixe,
                         x,width,x1begb,x1endb);
            }
	    scmd=NoCmd;
            return(SeismicPick(temp_tw));
	    break;
	case B2Up:
	    showloc = 0;
	    break;
	case B2Motion:
            if (showloc)
               xMouseLoc(dpy,win,event,style,True,
                         x,y,width,height,x1begb,x1endb,
                         x2begb,x2endb,p2beg,p2end);
	    break;
	case KDown:
            char keybuf[256];
            XLookupString(&(event.xkey),keybuf,0,&keysym,&keystat);

	    //another big switch statement for keysym
	    switch (keysym) {
		case XK_s:
		    xMousePrint(event,style, stdout,
                                x,y,width,height,
                                x1begb,x1endb,x2begb,x2endb,
                                p2beg, p2end, ss);
		    break;
		case XK_l:
		    lock=1; 
		    if (verbose) {
 			if (ss != NULL) (*ss) << "zoom lock set "<<lock<<endl;
			else cout << "zoom lock set "<<lock<<endl;
		    }
		    break;
		case XK_u:
		    lock=0;
		    if (verbose) {
			if (ss != NULL) (*ss) << "zoom lock released "<<lock<<endl;
			else cout <<"zoom lock released "<<lock<<endl;
		    }
		    break;	
		case XK_KP_1:
		case XK_1:
		    handle_zoom_keys(1.0,ss);
		    break;
 		case XK_KP_2:
                case XK_2:
                    handle_zoom_keys(2.0,ss);
                    break;
		case XK_KP_3:
                case XK_3:
                    handle_zoom_keys(3.0,ss);
                    break;
		case XK_KP_4:
                case XK_4:
                    handle_zoom_keys(4.0,ss);
                    break;
		case XK_KP_5:
                case XK_5:
                    handle_zoom_keys(5.0,ss);
                    break;
		case XK_KP_6:
                case XK_6:
                    handle_zoom_keys(6.0,ss);
                    break;
		case XK_KP_7:
                case XK_7:
                    handle_zoom_keys(7.0,ss);
                    break;
		case XK_KP_8:
                case XK_8:
                    handle_zoom_keys(8.0,ss);
                    break;
		case XK_KP_9:
                case XK_9:
                    handle_zoom_keys(9.0,ss);
                    break;
		case XK_KP_Left:
		case XK_Left:
		case XK_f:
		    handle_zoom_move(ZLeft,ss);
		    break;
                case XK_KP_Right:
                case XK_Right:
		case XK_h:
                    handle_zoom_move(ZRight,ss);
                    break;
                case XK_KP_Up:
                case XK_Up:
		case XK_t:
                    handle_zoom_move(ZUp,ss);
                    break;
                case XK_KP_Down:
                case XK_Down:
		case XK_b:
                    handle_zoom_move(ZDown,ss);
                    break;
	 	case XK_KP_Subtract:
		case XK_o:
		    handle_zoom_move(ZOut,ss);
		    break;
		case XK_KP_Add:
		case XK_i:
		    handle_zoom_move(ZIn,ss);
		    break; 
		case XK_Page_Down:
		case XK_c:
                    /* Change clip for image */
                    clip += clip/10. ;
                    /* note that image is out of date */
                    imageOutOfDate = 1;
		    break;
		case XK_Page_Up:
		case XK_a:
                    /* Change clip for image */
                    clip -= clip/10. ;
                    /* note that image is out of date */
                    imageOutOfDate = 1;
		    break;
		default:
		    break;
	    }
	    break;
	default:
	    break;
    }

    } catch (SeisppError) {
        if (ss != NULL) (*ss) << "Image Rendering Error in render_image"<<endl;
        else cerr << "Image Rendering Error in render_image"<<endl;
        exit(-3);
    }

}


void SeismicPlot::handle_zoom_keys(float zoom_level, stringstream * ss)
{
    mvefac=zoom_level;
    if (ss != NULL) (*ss)<<"Zoom/Move factor ="<<mvefac<<endl;
    else fprintf(stderr,"Zoom/Move factor = %f\n",mvefac);
}

void SeismicPlot::handle_zoom_move(SplotDirectionType direction, stringstream * ss)
{
    float mve;
    char *msg;

    switch (direction) {
	case ZDown:
            mve = (x2endb - x2begb)/mvefac ;
            x2begb = x2begb - mve ;
            x2endb = x2endb - mve ;
            msg="move ";
            /* check for bounds of full window */
            if (x2begb < x2beg) {
                if ( lock ) { 
		    x2begb = x2begb + mve ;
                    x2endb = x2endb + mve ;
                    msg="limit ";
                    mve=0;
                } else { 
		    x2begb = x2beg ;
		}
            }
	    break;
	case ZUp: 
            mve = (x2endb - x2begb)/mvefac ;
            x2begb = x2begb + mve ;
            x2endb = x2endb + mve ;
            msg="move ";
            /* check for bounds of full window */
            if (x2endb > x2end) {
                if ( lock ) { 
		    x2begb = x2begb - mve ;
                    x2endb = x2endb - mve ;
                    msg="limit ";
                    mve=0;
                } else { 
		    x2endb = x2end ;
		}
            }
	    break;
	case ZRight:
            mve = (x1endb - x1begb)/mvefac ;
            x1begb = x1begb + mve ;
            x1endb = x1endb + mve ;
            msg="move ";
            /* check for bounds of full window */
            if (x1endb > x1end) {
                if ( lock ) {  
		    x1begb = x1begb - mve ;
                    x1endb = x1endb - mve ;
                    msg="limit ";
                    mve=0;
                } else { 
		    x1endb = x1end ;
		}
            }
	    break;
	case ZLeft:
            /* move zoom box down by half window height */
            mve = (x1endb - x1begb)/mvefac ;
            x1begb = x1begb - mve ;
            x1endb = x1endb - mve ;
            msg="move ";
            /* check for bounds of full window */
            if (x1begb < x1beg) {
                if ( lock ) { 
		    x1begb = x1begb + mve ;
                    x1endb = x1endb + mve ;
                    msg="limit ";
                    mve=0;
                } else { x1begb = x1beg ;}
            }
	    break;
	case ZOut:
            /* zoom out .... vertical*/
            mve = (x1endb - x1begb)/mvefac ;
            x1begb = x1begb - mve ;
            x1endb = x1endb + mve ;
            /* check for bounds of full window */
            if (x1begb < x1beg) x1begb = x1beg ;
            if (x1endb > x1end) x1endb = x1end ;
            /*   .... and horizontal */
            mve = (x2endb - x2begb)/mvefac ;
            x2begb = x2begb - mve ;
            x2endb = x2endb + mve ;
            /* check bounds of original image */
            if (x2begb < x2beg) x2begb = x2beg ;
            if (x2endb > x2end) x2endb = x2end ;
	    msg="zoom out ";
	    break;
	case ZIn:
            /* zoom in .... vertical*/
            mve = (x1endb - x1begb)/(2.*mvefac) ;
            x1begb = x1begb + mve ;
            x1endb = x1endb - mve ;
            /*   .... and horizontal */
            mve = (x2endb - x2begb)/(2.*mvefac) ;
            x2begb = x2begb + mve ;
            x2endb = x2endb - mve ;
	    msg="zoom in ";
	    break;
	default:
	    if (ss != NULL) (*ss)<<"Wrong direction specification "<<direction<<endl;
	    else cerr<<"Wrong direction specification "<<direction<<endl;
	    break;
    }

    if (verbose) {
	if (ss!=NULL) (*ss)<<msg<<" "<<mve<<endl;
        else fprintf(stderr,"%s %g\n",msg,mve);
    }

    /* clear area and force an expose event */
    XClearArea(dpy,win,0,0,0,0,True);
    /* note that image is out of date */
    imageOutOfDate = 1;

}

PointPick SeismicPlot::pick_point()
{
	SeismicPick p;
	p=this->pick();
	while(p.type != POINT)
	{
		cerr << "Pick error:  use mouse button 2 to pick a point"<<endl;
		p=this->pick();
	}
	return(p.get_point());
}
TimeWindow SeismicPlot::pick_time_window()
{
	SeismicPick p;
	p=this->pick();
	while(p.type != WINDOW)
	{
		cerr << "Pick error:  click and drag mouse button 3 to pick a time window"<<endl;
		p=this->pick();
	}
	return(p.get_window());
}
int SeismicPlot::select_member()
{
	SeismicPick p;
	p=this->pick();
	while(p.type != POINT)
	{
		cerr << "Pick error:  use mouse button 2 to select a trace"<<endl;
		p=this->pick();
	}
	return(p.get_trace_number());
}

vector<int> SeismicPlot::select_members()
{
	SeismicPick p;
	vector<int> v;
	
	p=pick();
	while (p.type != PICKEND) {
	    if (p.type != POINT) {
		cerr << "Pick error: please use middle mouse button to select a trace"<<endl;
	    } else {
		v.push_back(p.get_trace_number());
	    }
	    p=pick();
	}
	
	return v;
}

int SeismicPlot::get_trace_number(float x2raw)
{
	int i;
	int result;
	float dx2min;


	if(use_variable_trace_spacing)
	{
		// In this case search x2 array for closest match
		// to x2raw
		result=0;
		dx2min=static_cast<float>(fabs(static_cast<double>(x2[0]-x2raw)));

		for(i=1;i<nmember;++i)
		{
			float dx2;
			dx2=static_cast<float>(fabs(static_cast<double>(x2[i]-x2raw)));
			if(dx2<dx2min)
			{
				result=i;
				dx2min=dx2;
			}
		}
	}
	else
	{
		// Block for regular spacing.  Simple formula will
		// do
		result=nint(x2raw/d2);
	}
	// Return result in C form, not FORTRAN indexing.  
	// Plot uses 1,2, ..., but this will be less confusing in C++
	// silently use edges if outside bounds
	if(result<=0)
		return(0);
	else if(result>=nmember)
		return(nmember-1);
	else
		return(result-1);
}
SeismicPick SeismicPlot::pick()
{
	int i,j;
	int i1,i2;
	float mve;
	char *msg;


	/* main event loop */
	XEvent event;
	if(image==NULL)
		image=NULL;
	else if(imageOutOfDate)
		XDestroyImage(image);
		
	// This seens to force an initialization on first pass of these variables
	int winwidth=-1;  
	int winheight=-1;

	//get an EXPOSE event in first
	XClearArea(dpy,win,0,0,0,0,True);

	while(imageOutOfDate|(~imageOutOfDate)/*True*/) {
		XNextEvent(dpy,&event);

		/* if window was resized */
		if (event.type==ConfigureNotify &&
			(event.xconfigure.width!=winwidth ||
			 event.xconfigure.height!=winheight)) {
			winwidth = event.xconfigure.width;
			winheight = event.xconfigure.height;
							
			/* determine good size for axes box */
			xSizeAxesBox(dpy,win,
				labelfont,titlefont,style,
				&x,&y,&width,&height);
	
			/* clear the window */
			XClearWindow(dpy,win);

                        /* force an expose event */
                        XClearArea(dpy,win,0,0,0,0,True);
			
			/* note that image is out of date */
			imageOutOfDate = 1;

		} 
		else if (event.type==Expose) 
		{
		  
			/* clear all expose events from queue */
			while (XCheckTypedEvent(dpy,Expose,&event));
			
			/* if necessary, make new image */
			if (imageOutOfDate) {
				if (image!=NULL) {
				/*	free1(image->data); */
					XDestroyImage(image);
				}
				image = newBitmap(dpy,width,height,
					n1,d1,f1,n2,x2,z,
					x1begb,x1endb,x2begb,x2endb,
					xcur,clip,wt,va,
					&p2beg,&p2end,endian,interp,
					wigclip,style);
				imageOutOfDate = 0;
			}
	
			/* draw image (before axes so grid lines visible) */
			XPutImage(dpy,win,gci,image,0,0,x,y,
				image->width,image->height);


			/* draw axes on top of image */
			xDrawAxesBox(dpy,win,
				x,y,width,height,
				x1begb,x1endb,0.0,0.0,
				d1num,f1num,n1tic,grid1,label1,
				x2begb,x2endb,p2beg,p2end,
				d2num,f2num,n2tic,grid2,label2,
				labelfont,title,titlefont,
				labelcolor,titlecolor,gridcolor,
				style);

		  /* else if key down */
		} else if (event.type==KeyPress) {
		  
			char keybuf[256];
			XLookupString(&(event.xkey),keybuf,0,&keysym,&keystat);
			if (keysym==XK_s) {
			  xMousePrint(event,style, stdout,
				      x,y,width,height,
				      x1begb,x1endb,x2begb,x2endb,
				      p2beg, p2end,NULL);
			  
				
			} else if (keysym==XK_l ) {
				/* set lock */		  
			  lock = 1 ;
			  if (verbose) fprintf(stdout,"zoom lock set  %d\n",lock);

 			} else if (keysym==XK_u ) {
				/* unset lock */		  
			  lock = 0 ;
			  if (verbose) fprintf(stdout,"zoom lock released %d\n",lock);
			/* added by Peng Wang, return an end SeismicPick type so that
			   we could pick consequetively some members to delete */
			} else if (keysym==XK_t) {
			   SeismicPick spick;
			   spick.type=PICKEND;
			   return spick;
			/* end changes by Peng Wang */
 			} else if (keysym==XK_Shift_L ) { 
			  /* if (verbose) 
			     fprintf(stderr,"Shift Left pressed \n");*/
			} else if (keysym==XK_KP_1 || keysym==XK_1 ) { 
			  mvefac=1.;
			  fprintf(stderr,"Zoom/Move factor = 1 \n");
			} else if (keysym==XK_KP_2 || keysym==XK_2 ) { 
			  mvefac=2.;
			  fprintf(stderr,"Zoom/Move factor = 2 \n");
			} else if (keysym==XK_KP_3 || keysym==XK_3 ) { 
			  mvefac=3.;
			  if (verbose) 
			    fprintf(stderr,"Zoom/Move factor = 3 \n");
			} else if (keysym==XK_KP_4 || keysym==XK_4 ) { 
			     mvefac=4.;
			     if (verbose) 
			       fprintf(stderr,"Zoom/Move factor = 4 \n");
			} else if (keysym==XK_KP_5 || keysym==XK_5 ) { 
			  mvefac=5.;
			  if (verbose) 
			    fprintf(stderr,"Zoom/Move factor = 5 \n");
			} else if (keysym==XK_KP_6 || keysym==XK_6 ) { 
			  mvefac=6.;
			     if (verbose) 
			       fprintf(stderr,"Zoom/Move factor = 6 \n");
			} else if (keysym==XK_KP_7 || keysym==XK_7 ) { 
			  mvefac=7.;
			  if (verbose) 
			    fprintf(stderr,"Zoom/Move factor = 7 \n");
			} else if (keysym==XK_KP_8 || keysym==XK_8 ) { 
			  mvefac=8.;
			  if (verbose) 
			    fprintf(stderr,"Zoom/Move factor = 8\n");
			} else if (keysym==XK_KP_9 || keysym==XK_9 ) { 
			  mvefac=9.;
			  if (verbose) 
			    fprintf(stderr,"Zoom/Move factor = 9\n");
			} else if (keysym==XK_Left || keysym==XK_KP_Left ) {
			  
 			  /* move zoom box to left by half window width */
			  mve = (x2endb - x2begb)/mvefac ;
			  x2begb = x2begb - mve ;
			  x2endb = x2endb - mve ;
			  msg="move "; 
			  /* check for bounds of full window */
			  if (x2begb < x2beg) {
			    if ( lock ) { x2begb = x2begb + mve ;
			    x2endb = x2endb + mve ;
			    msg="limit ";
			    mve=0;
			    } else { x2begb = x2beg ;}
			  }
			  
			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);
			  
			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);
			  /* note that image is out of date */
			  imageOutOfDate = 1;
								
			} else if (keysym==XK_Right || keysym==XK_KP_Right ) {
			  /* move zoom box to right by half window width*/
			  mve = (x2endb - x2begb)/mvefac ;
			  x2begb = x2begb + mve ;
			  x2endb = x2endb + mve ;
			  msg="move "; 
			  /* check for bounds of full window */
			  if (x2endb > x2end) {
			    if ( lock ) { x2begb = x2begb - mve ;
			                  x2endb = x2endb - mve ;
					  msg="limit ";
					  mve=0;
			    } else { x2endb = x2end ;}
			  }
			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);

			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);
			  /* note that image is out of date */
			  imageOutOfDate = 1;
								
			} else if (keysym==XK_Down || keysym==XK_KP_Down  ) {
			  /* move zoom box down by half window height */
			  mve = (x1endb - x1begb)/mvefac ;
			  x1begb = x1begb + mve ;
			  x1endb = x1endb + mve ;
			  msg="move "; 
			  /* check for bounds of full window */
			  if (x1endb > x1end) {
			    if ( lock ) {  x1begb = x1begb - mve ;
			                   x1endb = x1endb - mve ;
					  msg="limit ";
					  mve=0;
			    } else { x1endb = x1end ;}
			  }
			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);

			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);
			  /* note that image is out of date */
			  imageOutOfDate = 1;
								
			} else if (keysym==XK_Up || keysym==XK_KP_Up ) {
			  /* move zoom box down by half window height */
			  mve = (x1endb - x1begb)/mvefac ;
			  x1begb = x1begb - mve ;
			  x1endb = x1endb - mve ;
			  msg="move "; 
			  /* check for bounds of full window */
			  if (x1begb < x1beg) {
			    if ( lock ) { x1begb = x1begb + mve ;
			                  x1endb = x1endb + mve ;
					  msg="limit ";
					  mve=0;
			    } else { x1begb = x1beg ;}
			  }
			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);
				
			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);
			
			/* note that image is out of date */
			imageOutOfDate = 1;
								
			} else if (keysym==XK_o || keysym==XK_KP_Subtract ) {
			  /* zoom out .... vertical*/
			  mve = (x1endb - x1begb)/mvefac ;
			  x1begb = x1begb - mve ;
			  x1endb = x1endb + mve ;
			  /* check for bounds of full window */
			  if (x1begb < x1beg) x1begb = x1beg ;
			  if (x1endb > x1end) x1endb = x1end ;
			  /*   .... and horizontal */
			  mve = (x2endb - x2begb)/mvefac ;
			  x2begb = x2begb - mve ;
			  x2endb = x2endb + mve ;
			  /* check bounds of original image */
			  if (x2begb < x2beg) x2begb = x2beg ;
			  if (x2endb > x2end) x2endb = x2end ;
			 	
			  /* clear area and force an expose event */
		   	  XClearArea(dpy,win,0,0,0,0,True);
			 
			  /* note that image is out of date */
			  imageOutOfDate = 1;
								
			} else if (keysym==XK_i || keysym==XK_KP_Add ) {
			  /* zoom in .... vertical*/
			  mve = (x1endb - x1begb)/(2.*mvefac) ;
			  x1begb = x1begb + mve ;
			  x1endb = x1endb - mve ;
			  /*   .... and horizontal */
			  mve = (x2endb - x2begb)/(2.*mvefac) ;
			  x2begb = x2begb + mve ;
			  x2endb = x2endb - mve ;

			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);
			 
			  /* note that image is out of date */
			  imageOutOfDate = 1;
								
			} else if (keysym==XK_c || keysym==XK_Page_Down) {
		  		
				/* Change clip for image */
 		       		clip += clip/10. ;
				if (verbose) fprintf(stderr,"clip=%g\n",clip);
 				/* note that image is out of date */
				 imageOutOfDate = 1;				
				 
			} else if (keysym==XK_a || keysym==XK_Page_Up) {

				/* Change clip for image */
			        clip -= clip/10. ;
				if (verbose) fprintf(stderr,"clip=%g\n",clip);
				/* note that image is out of date */
				imageOutOfDate = 1;				

			} else {
				continue;
			}

		/* else if button down (1 == zoom, 2 == mouse tracking, 3==window pick */
		} else if (event.type==ButtonPress) {
			/* if 1st button: zoom */
			if (event.xbutton.button==Button1) {
				int xb,yb,wb,hb;

				/* track pointer and get new box */
				xRubberBox(dpy,win,event,&xb,&yb,&wb,&hb);
			
				/* if new box has tiny width or height */
				if (wb<4 || hb<4) {
				
					/* reset box to initial values */
					x1begb = x1beg;
					x1endb = x1end;
					x2begb = x2beg;
					x2endb = x2end;
			
				/* else, if new box has non-zero width */
				/* if new box has zero width or height */
				} else {
			
					/* calculate new box parameters */
					zoomBox(x,y,width,height,
						xb,yb,wb,hb,
						x2begb,x2endb,
						x1begb,x1endb,
						&x2begb,&x2endb,
						&x1begb,&x1endb,
                                                style);
				}

				/* clear area and force an expose event */
				XClearArea(dpy,win,0,0,0,0,True);
			
				/* note that image is out of date */
				imageOutOfDate = 1;
			
			/* else if 2nd button down: display mouse coords */
			} else if (event.xbutton.button==Button2) {

				
				SeismicPick ppick=MouseLocate(dpy,win,
					  event,style,
					  x,y,width,height,x1begb,x1endb,
					  x2begb,x2endb,p2beg,p2end);
				ppick.trace_number=
					this->get_trace_number(ppick.x2);
				cout << ppick.trace_number<<endl;
				int isample=this->member[ppick
					.trace_number]
					.sample_number((double)ppick.x1);

				ppick.amplitude=this->member[ppick.trace_number].s[isample];
				return(ppick);
		

			// To pick a window use button 3 
			} else if (event.xbutton.button==Button3) {
				int ixb,ixe;
				int wnow;
				pick_x1_interval(dpy,win,event,style,
					x,y,width,height,
					&ixb,&ixe,NULL, NULL);
				// Not correct, quick and dirty for testing
				TimeWindow tw;
				if(style==SEISMIC)
				{
					tw.start=screen_to_time(ixb,
						y,height,x1begb,x1endb);
					tw.end=screen_to_time(ixe,
						y,height,x1begb,x1endb);
				}
				else
				{
					tw.start=screen_to_time(ixb,
						x,width,x1begb,x1endb);
					tw.end=screen_to_time(ixe,
						x,width,x1begb,x1endb);
				}
				return(SeismicPick(tw));

			} else {
				continue;
			}

		/* else if pointer has moved */
		} else if (event.type==MotionNotify) {
			
			/* if button2 down, show mouse location */
			if (showloc)
				xMouseLoc(dpy,win,event,style,True,
					x,y,width,height,x1begb,x1endb,
					x2begb,x2endb,p2beg,p2end);

		/* else if button2 released, stop tracking */
		} else if (event.type==ButtonRelease &&
			   event.xbutton.button==Button2) {
			showloc = 0;
		}

	} /* end of event loop */
}

}

