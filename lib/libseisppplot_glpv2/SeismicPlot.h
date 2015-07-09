#ifndef SEISMICPLOT_H
#define SEISMICPLOT_H


#include "xplot.h"
#include <X11/Xatom.h>
#include <X11/keysym.h>

// Motif includes, used in adjust_scrollbar call 
#include <Xm/Xm.h>
#include <Xm/ScrollBar.h>

#include "seispp.h"
#include <string>
#include <vector>

#include <sstream>

namespace SEISPP
{
using namespace std;
using namespace SEISPP;

typedef struct PointPick {
	double time;
	double amplitude;
} PointPick;

enum PickType {POINT, WINDOW, PICKEND, UNDEFINED};

enum SplotEventType {B1Down, B2Down, B3Down, B1Motion, B2Motion, B3Motion, B1Up, B2Up, B3Up, KDown, KUp};
enum SplotDirectionType {ZLeft, ZRight, ZUp, ZDown, ZIn, ZOut};

enum SplotCommandType {NoCmd, SinglePick, MultiplePicks};

class SeismicPick {
public:
	PickType type;
	double time;
	double amplitude;
	int trace_number;
	TimeWindow twin;

	SeismicPick();
	SeismicPick(float x1in,float x2in);
	SeismicPick(TimeWindow tw);
	SeismicPick(const SeismicPick& p);
	SeismicPick& operator=(const SeismicPick& p);
	TimeWindow get_window();
	PointPick get_point();
	int get_trace_number();
	void set_point(double t, double a);
	friend class SeismicPlot;
private:
	float x1,x2;
	bool point_set;
	bool window_set;
};
	

class SeismicPlot : public TimeSeriesEnsemble
{
public:
	//Basic constructor builds a window with no data
	SeismicPlot();
	//A semi-copy constructor.  Data loaded from ensemble, plot parameters
	//extracted from Metadata object
	//Design note:  should add parameters from Metadata to object metatdata are
	//for flexibility -- methods can extract some parameters metadata component.
	SeismicPlot(TimeSeriesEnsemble &,Metadata);
	// Similar to above, but extract one component to build display
	SeismicPlot(ThreeComponentEnsemble & ,int,Metadata);
	//destructor is nontrivial and cannot be defaulted
	~SeismicPlot();
	// set plot parameters 
	void SetParameters(Metadata);

	int get_width() {return width;}
	int get_height() {return height;}

	// These two methods modify only the ensemble data.  Design detail is
	// undecided -- may or may not redraw window.
	void LoadNewData(TimeSeriesEnsemble &);
	void LoadNewData(ThreeComponentEnsemble &,int);
	void draw();
	SeismicPick pick();
	TimeWindow pick_time_window();
	PointPick pick_point();
	int select_member();
	vector<int> select_members();
	void remove_traces(vector<int>);
	void remove_trace(int);

private:

	// standardized start function called by constructors.
	void X_init();
	// Peng had this public.  Not desirable in public interface.
	// Useful only internally.
	void adjust_nmember() {nmember=member.size();}

	// Internal function to compute member number.  simple
	// for equal spaced traces, tough with variable offset.
	int get_trace_number(float x2val);
	// These are graphics parameters.  
	// Most are direct interface to SU routines, some were altered
	// by glp for clarity 
	int xbox,ybox,wbox,hbox;
	// assorted drawing parameters
	char *labelfont,*titlefont;
	char *labelcolor,*titlecolor;
	char *label1,*label2,*title,*windowtitle;
	float labelsize,titlesize;
	char *gridcolor;
	int style;
	int grid1,grid2;
	bool use_variable_trace_spacing;
	string trace_axis_attribute;  // used only when use_variable_trace_spacing true
	double trace_spacing;  // d2 in SU code.  Used only when not variable
	double first_trace_offset;  // f2 in SU code
	int wt,va;
	int n1tic;
	int n2tic;
	// GLP addition.  
	string time_scaling;  // If "auto" use data range, otherwise use x1beg and x1end for time range
	string trace_axis_scaling;  // If "auto" use data range, otherwise use x2beg and x2end for time range
	// These are only used when scaling is set as fixed
	float x1beg,x1end,x2beg,x2end;
	float d1num,f1num,d2num,f2num;
	int interp;
	bool clip_data;
	float perc;
	int wigclip;
	float xcur;
	// This was originally a char **
	vector<string>(curvecolor);
	string plotfile;
	bool verbose;
	string default_curve_color;
	int endian;  // Established by query, but conveniently stored.
	//
	// These are X11 related quantities used by SU xplot
	// Others probably need to be added.  
	//
	Display *dpy;
	Window win;
	Widget wdgt;  // WARNING:  Not currently handled correctly
	int scr;
	KeySym keysym;
	XComposeStatus keystat;
	XImage *image;
	GC gci; 
	// x2 is position of zero line for each trace.  Length nmember
	float *x2;   
	// From here on are draw variables needed to pass relevant
	// data between calls.  In a procedural version they might
	// be thought of as extern variables or at least variables
	// in file scope.
	int nmember;
	float d1,d2,f1,f2;
        float x1begb,x1endb,x2begb,x2endb;
        int x,y,width,height;
	int n1,n2;
	float *z;  // float buffer, waste of space. 
	float clip;
	int imageOutOfDate;
	// Moved to private by GLP.  We can't have X specific
	// variables in the public interface as X is too implementation
	// specific.  We must hide these in the private interface if
	// they must be in this object.
	
	/*************************************************
	// Old prototype from Peng 
	SeismicPick process_callbacks(Display *, Window, SplotEventType, XEvent *, 
		stringstream *, Widget *);
	*******************************************/
	/**************** Changed 12/30/2005 by GLP to this */
	SeismicPick process_callbacks(SplotEventType, XEvent *,stringstream *);
	void process_commands(SplotCommandType, stringstream *);
	SplotCommandType get_command() {return scmd;}
	// this should be private also more because render_image is 
	// not general enough.  One might choose vector graphics to
	// display the data.
 	void render_image();

	//added by Peng Wang, 9/20/2005
	void handle_zoom_keys(float, stringstream *);
	void handle_zoom_move(SplotDirectionType, stringstream *);
	SplotCommandType scmd;
	int showloc;
	float p2beg, p2end;
	int lock;
	float mvefac;
};

} // End Namespace SEISPP
#endif
