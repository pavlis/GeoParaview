#include <string>

#include "elog.h"
#include "xplot.h"
#include "gclgrid.h"
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include "metadata.h"
#include "seispp.h"
#include "SeismicPlot.h"
#define MAXPLOTSAMPLES 100000
using namespace std;
using namespace SEISPP;

namespace SEISPP
{

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
	i1beg = SEISPP::nint((x1beg-f1)/d1);
	i1beg = MAX(0,MIN(n1-1,i1beg));
	x1beg = f1+i1beg*d1;
	i1end = SEISPP::nint((x1end-f1)/d1);
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

void xMousePrint(XEvent event, int style, FILE *mpicksfp,
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
	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = p2beg+x2begb+(p2end+x2endb-x2begb-p2beg)*
			(event.xmotion.x-x)/width;
	}

	/* write string indicating mouse location */
	fprintf(mpicksfp, "%0.6g  %0.6g\n", x1, x2);
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
		bval=md.get_bool("WingleTrace");
		if(bval)
			wt=1;
		else
			wt=0;
		bval=md.get_bool("VariableArea");
		if(bval)
			va=1;
		else
			va=0;
		
		clip_data=md.get_bool("clip_data");
		if(clip_data)
			clip=static_cast<float>(md.get_double("clip_level"));
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
		if(s=="seismic") 
			style=SEISMIC;
		else
			style=NORMAL;
			
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
		s=md.get_string("styles");
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
}
SeismicPlot::SeismicPlot(TimeSeriesEnsemble& tse, Metadata md)  : TimeSeriesEnsemble(tse)
{
	int i;
	int nmember=member.size();
	this->SetParameters(md);
	for(i=0;i<nmember;++i)
		curvecolor.push_back(default_curve_color);
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
	/* connect to X server */
	if ((dpy=XOpenDisplay(NULL))==NULL)
		throw SeisppError("XOpenDisplay failed");
	scr = DefaultScreen(dpy);
	unsigned long black, white;
	black = BlackPixel(dpy,scr);
	white = WhitePixel(dpy,scr);
	
	/* set endian for display.  Original a bit more complex and involved
	an apparent default passed as a parameter.  This allowed bypass of 
	BitmapBitOrder call below.*/
	if(BitmapBitOrder(dpy)==LSBFirst)
		endian=0;
	else if(BitmapBitOrder(dpy)==MSBFirst)
		endian=1;

	/* create window */
	win = xNewWindow(dpy,xbox,ybox,wbox*member.size(),
		hbox,(int) black,(int) white,windowtitle);
		
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

SeismicPlot::SeismicPlot(ThreeComponentEnsemble& tce, int comp, Metadata md) 
{
	int i;
	int nmember=tce.member.size();
	this->SetParameters(md);
	member.reserve(nmember);
	// this loads component comp as the member to be plotted.
	for(int i=0;i<nmember;++i)
	{
		TimeSeries *c;
		c=ExtractComponent(tce.member[i],comp);
		member.push_back(*c);
		delete c;
		curvecolor.push_back(default_curve_color);
	}
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
	/* connect to X server */
	if ((dpy=XOpenDisplay(NULL))==NULL)
		throw SeisppError("XOpenDisplay failed");
	scr = DefaultScreen(dpy);
	unsigned long black, white;
	black = BlackPixel(dpy,scr);
	white = WhitePixel(dpy,scr);
	
	/* set endian for display.  Original a bit more complex and involved
	an apparent default passed as a parameter.  This allowed bypass of 
	BitmapBitOrder call below.*/
	if(BitmapBitOrder(dpy)==LSBFirst)
		endian=0;
	else if (BitmapBitOrder(dpy)==MSBFirst)
		endian=1;

	/* create window */
	win = xNewWindow(dpy,xbox,ybox,wbox*member.size(),
		hbox,(int) black,(int) white,windowtitle);
		
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


void SeismicPlot::LoadNewData(ThreeComponentEnsemble& tse,int comp) 
{
	member.resize(tse.member.size());
	for(int i=0;i<tse.member.size();++i)
	{
		TimeSeries *c;
		c=ExtractComponent(tse.member[i],comp);
		member[i]=(*c);
		delete c;
	}
}

void SeismicPlot::LoadNewData(TimeSeriesEnsemble& tse)
{
	member.resize(tse.member.size());
	for(int i=0;i<tse.member.size();++i)
		member[i]=tse.member[i];
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
		if(data[j].endtime() < t0)
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
					if( (id>=0) && (id<data[j].ns) )
						z[iz]=static_cast<float>(data[j].s[i]);
				}
			}
		}
	}
}
void SeismicPlot::draw() {
	int i,j;
	int i1,i2;
	int nmembers=member.size();
	if(nmembers<=0) throw SeisppError("Seisplot::draw not data to plot\n");
	// First compute range of x1 and x2
	double x1min,x1max,x2min,x2max;
	float mve,mvefac;
	int lock;
	char *msg;
	float d2,f1,f2;

	x1min=member[0].t0;
	x1max=member[0].endtime();
	for(i=0;i<nmembers;++i)
	{
		x1min=min(member[i].t0,x1min);
		x1max=max(member[i].endtime(),x1max);
	}
	x2min=member[0].s[0];
	x2max=member[0].s[0];
	// minor inefficiency here, but negligible 
	for(j=0;j<nmembers;++j)
	{
		for(i=0;i<member[j].ns;++i)
		{
			x2min=min(member[j].s[i],x2min);
			x2max=max(member[j].s[i],x2max);
		}
	}
	
	if(style==SEISMIC) 
	    
	    d2=(x2max-x2min)/wbox*nmembers;
	else
	    d2=(x2max-x2min)/hbox*nmembers;
  
  
  
	//
	// Set x1beg and x1end for auto scaling.  Not else currently
	// because in that case assume x1beg and x1end are set in construction
	//
	if(time_scaling=="auto")
	{
		x1beg=x1min;
		x1end=x1max;
	}
	f1=x1beg;
	//
	// parent code set a bunch of parameters now stored with the
	// SeismicPlot object.  Look at xwigb main for original 
	//
	/* initialize zoom box parameters */
	float x1begb,x1endb,x2begb,x2endb;
	x1begb = static_cast<float>(x1beg);	 
	x1endb = static_cast<float>(x1end);
	x2begb = static_cast<float>(x2beg);	 
	x2endb = static_cast<float>(x2end);

					
	/* determine good size for axes box */
	int x,y,width,height;
	xSizeAxesBox(dpy,win,
		labelfont,titlefont,style,
		&x,&y,&width,&height);
	/* note that image is out of date */
	int imageOutOfDate = 1;
	int width1=width;
	int height1=height;
	// Determine size.  Depends on setting of time_scaling.  If auto get it fromt he
	// data, otherwise use the internal parameters.
	//
	// These were parameters in the original code.  This is a translation 
	// from seispp symbols to symbols used in original SU code
	int n1,n2;
	float d1=static_cast<float>(member[0].dt);
	if(time_scaling=="auto")
	{
		n1=SEISPP::nint( (x1max-x1min)/d1 );
		n2=member.size();
	}
	else
	{
		n1=SEISPP::nint((x1beg-x1end)/d1);
		n2=member.size();
	}
	//
	// Get the number of samples from x1beg and x1end.  
	// With above logic this should work if time_scaling is auto or manual
	//
	n1=SEISPP::nint((x1beg-x1end)/d1);
        n2=member.size();

	// Sanity checks on n1
	if(n1<2)
	{
		char message[256];
		sprintf(message,"SeismicPlot::RasterDraw:  Plot range error.\n%d samples in requested time range of %lf to %lf\n",
				n1,x1beg,x1end);
		throw SeisppError(string(message));
	}
	else if(n1>MAXPLOTSAMPLES)
	{
		char message[128];
		sprintf(message,"SeismicPlot::RasterPlot: Plot range request too large\nTime range of %lf to %lf requires %n samples.  Sanity check overrides request.  Check input parameters\n",
			x1beg,x1end,n1);
		throw SeisppError(string(message));
	}
	// d2 and f2 are defined in this object by less obscure names
	d2=static_cast<float>(this->trace_spacing);
	f2=static_cast<float>(this->first_trace_offset);
	// This buffer is needed in plot loop below.  We construct a fortran-type matrix
	// of floats with the data aligned by member[i].t0 values.  Traces are in columns
	// of the matrix.  Depend on new throwing an exception if alloc fails.
	//
	float *z=new float[n1*n2];
	// internal function used here
	load_z_matrix(this->member,z,n1,n2,x1beg,d1);
	// Original code had these as a pair of 2d arrays.  Seems here they are redundant
	// so I'm going to make them simple float vectors and create them from z when
	// needed
	float *x1curve=new float[n1];
	float *x2curve=new float[n1];
	for(i=0;i<n1;++i) x1curve[i]=x1beg+d1*static_cast<float>(i);
	/* main event loop */
	float p2beg,p2end;
	XEvent event;
	XImage *image=NULL;
	// This seens to force an initialization on first pass of these variables
	int winwidth=-1;  
	int winheight=-1;
	int showloc=0;
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

			/* draw curve on top of image.  Original code had x2curve
			independent of z.  Here I (glp) force x2curve to be 
			copied from z.  curvecolor will be a useful way to flag
			different data for different programs. */
			for (i=0; i<n2; i++)
			{
				scopy(n1,z+i*n1,1,x2curve,1);
				xDrawCurve(dpy,win,
					   x,y,width,height,
					   x1begb,x1endb,0.0,0.0,
					   x2begb,x2endb,p2beg,p2end,
					   x1curve,x2curve,n1,
					   const_cast<char *>(curvecolor[i].c_str()),
						style);
			}

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
				      p2beg, p2end);
			  
				
			} else if (keysym==XK_l ) {
				/* set lock */		  
			  lock = 1 ;
			  if (verbose) fprintf(stdout,"zoom lock set  %d\n",lock);

 			} else if (keysym==XK_u ) {
				/* unset lock */		  
			  lock = 0 ;
			  if (verbose) fprintf(stdout,"zoom lock released %d\n",lock);

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

			} else if (keysym==XK_q || keysym==XK_Q) {
			/* This is the exit from the event loop */
				break;
			} else if (keysym==XK_p || keysym==XK_P) {
			// This method requires SU to be installed.  Builds a system call
			// to pswigb
			/* invoke pswigb with appropriate data */
				char cmdtemp[256];
				
				char *cmdline;
				float cmdfloat;
				int nbpi;
				float num;
				short debug=1;

				FILE *plotfp;	/*fp for plot data*/
							
				cmdline = (char *) malloc(n2+BUFSIZ);
				if(cmdline==NULL) 
					throw SeisppError("RasterDraw:  malloc failure in postscript plot code segment");
				strcpy(cmdline,"pswigb ");
				//
				// SU code passed argv to get these.  We have to do 
				// it by translation
				sprintf(cmdtemp," d2num=%f ",d2num);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x1beg=%f ",x1beg);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x1end=%f ",x1end);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x2beg=%f ",x2beg);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x2end=%f ",x2end);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," title=%s ",title);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," label1=%s ",label1);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," label2=%s ",label2);
				strcat(cmdline,cmdtemp);
				if(style==SEISMIC)
					sprintf(cmdtemp," style=seismic");
				else
					sprintf(cmdtemp," style=normal");
				strcat(cmdline,cmdtemp);
				
				/* override incompatible args */
				sprintf(cmdtemp," axescolor=%s",labelcolor);
				strcat(cmdline,cmdtemp);
				// Variable in SU.  Here we fix resolution at 300 dpi
				nbpi = 300; 
				sprintf(cmdtemp," nbpi=%d",nbpi);
				strcat(cmdline,cmdtemp);
				cmdfloat = DisplayWidthMM(dpy,scr)/25.4;
				cmdfloat /= DisplayWidth(dpy,scr);
				sprintf(cmdtemp," wbox=%g", cmdfloat*width);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," xbox=%g", 0.5+cmdfloat*xbox);
				strcat(cmdline,cmdtemp);
				cmdfloat = DisplayHeightMM(dpy,scr)/25.4;
				cmdfloat /= DisplayHeight(dpy,scr);
				sprintf(cmdtemp," hbox=%g", cmdfloat*height);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," ybox=%g", 0.5+cmdfloat*ybox);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x1beg=%g", x1begb);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x1end=%g", x1endb);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x2beg=%g", x2begb);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x2end=%g", x2endb);
				strcat(cmdline,cmdtemp);
				num=(x2endb-x2begb)/4;
				sprintf(cmdtemp," d2num=%g", num);
				strcat(cmdline,cmdtemp);
				strcat(cmdline," style=normal ");
				strcat(cmdline," title=\"");
				strcat(cmdline,title); strcat(cmdline,"\"");
				strcat(cmdline," label1=\"");
				strcat(cmdline,label1); strcat(cmdline,"\"");
				strcat(cmdline," label2=\"");
				strcat(cmdline,label2); strcat(cmdline,"\"");
				sprintf(cmdtemp," > %s ", plotfile.c_str());
				strcat(cmdline,cmdtemp);
				// Slight modification from SU version.  SU version
				// used an epopen.  Here I use C++ exceptions to trap
				// problems.  popen forks the process, runs command cmdline,
				// with stdin connected to a pipe plotfp
				plotfp=popen(cmdline,"w");
				if(plotfp==NULL)
					throw SeisppError("SeismicPlot::RasterDraw:  failure of popen in attempting to write postscript file");

				free(cmdline);
				int nwritten;
				int nz=n1*n2;
				nwritten=fwrite(z,sizeof(float),nz,plotfp);
				pclose(plotfp);				
				if(nwritten!=nz)
					cerr << "SeismicPlot::RasterDraw (Warning):  postscript plot output problem.  Output appears to have been truncated"<<endl;
		
	
			} else {
				continue;
			}

		/* else if button down (1 == zoom, 2 == mouse tracking */
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

				showloc = 1;
				
				xMouseLoc(dpy,win,event,style,showloc,
					  x,y,width,height,x1begb,x1endb,
					  x2begb,x2endb,p2beg,p2end);
		

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
	delete [] z;
	delete [] x1curve;
	delete [] x2curve;
}
//
// Under construction note.  This code can be derived from SU xpicker program.  It has
// many similarities to xwigb used for the RasterDraw method above.  There are some 
// difference, however, so my idea is to not convert it until I get the simpler wigb
// conversion above done.  Then cross-comparision of the two main programs should make
// it easy to clone appropriate parts of the RasterDraw method above with changes from
// the xpicker added.
SeismicPick SeismicPlot::pick()
{
	cerr << "pick method not yet implemented"<<endl;
	return(SeismicPick());
}
}

