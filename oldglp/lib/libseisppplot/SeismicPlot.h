
#ifndef SEISMICPLOT_H
#define SEISMICPLOT_H


#include "xplot.h"
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include "seispp.h"
#include <string>
#include "SeismicPick.h"

using namespace SEISPP;

class SeismicPlot : TimeSeriesEnsemble
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
	// These two methods modify only the ensemble data.  Design detail is
	// undecided -- may or may not redraw window.
	void LoadNewData(TimeSeriesEnsemble &);
	void LoadNewData(ThreeComponentEnsemble &,int);
	// draws data.  Constructors build the window frame.
	void draw();
	// Generalized picker interface.  Return a SeismicPick object.
	SeismicPick pick();
private:
  // These are X11 related quantities used by SU xplot
  Display *dpy;
  Window win;
  KeySym keysym;
  XComposeStatus keystat;
  XImage *image;
  GC gci;
  int scr;  // screen number 
  // These are graphics parameters.  Defaults will be loaded from Metadata object
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
  float *x2;   
};

#endif
