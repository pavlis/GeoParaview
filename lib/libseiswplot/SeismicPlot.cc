#include <X11/Xthreads.h>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include "SeismicPlot.h"
using namespace std;
using namespace SEISPP;
#define APP_CLASS "seismicplot"
void continue_callback(Widget w, XtPointer client_data, XtPointer call_data)
{
    SeismicPlot *plot_handle=reinterpret_cast<SeismicPlot *>(client_data);

    plot_handle->ExitDisplay();
}
/* Default constructor gets its parameter information from the
standard antelope parameter file location under a frozen name 
(SeismicPlot.pf).  */
SeismicPlot::SeismicPlot()
{
    Pf *pf;
    const string pfglobal_name("SeismicPlot");
    if(pfread(const_cast<char *>(pfglobal_name.c_str()),&pf))
        throw SeisppError("SeismicPlot constructor:  pfread failed for "
                + pfglobal_name);
    try {
	Metadata md(pf);
	pffree(pf);
	*this=Metadata::operator=(md);
	/* from here down this code is identical to the 
	constructor with a Metadata argument.  Beware in
	maintenance */
        for(int k=0;k<3;++k)
        {
            rshell[k]=NULL;
            seisw[k]=NULL;
        }
        argc=1;
        argv[0]=strdup("SeismicPlotWidget");
        XtSetLanguageProc(NULL, NULL, NULL);
        XtToolkitThreadInitialize();
        toplevel = XtVaAppInitialize(&AppContext, (char *)"seismicplot",NULL,0,
                                &argc,argv, NULL,NULL);
        main_w=XmCreateForm(toplevel,(char *) "seismicplot",NULL,0);
        EventLoopIsActive=false;
        Arg args[5];
        int n=0;
        XmString str;
        const string button_string("Continue");
        str=XmStringCreateLocalized(const_cast<char *>(button_string.c_str()));
        n=0;
        XtSetArg (args[n], XmNlabelString, str); n++;
        XtSetArg (args[n], XmNmultiClick, XmMULTICLICK_DISCARD); n++;
        continue_button=XmCreatePushButton(main_w,
            const_cast<char *>(button_string.c_str()),args,n);
        XmStringFree(str);
        XtAddCallback(continue_button,XmNactivateCallback,continue_callback,
                (XtPointer)this);
        XtManageChild(continue_button);
        XtManageChild(main_w);
        XtRealizeWidget(toplevel);
    } catch (...) {throw;}
}
SeismicPlot::SeismicPlot(Metadata& md) : Metadata(md)
{
    try {
        /* This is a complete duplicate of above.  There is probably 
           a tricky way to handle this some other way, but for now we
           do it this way.  Be careful if changes are made to one to 
           change the other too. */
        /* These have to be initialized */
        for(int k=0;k<3;++k)
        {
            rshell[k]=NULL;
            seisw[k]=NULL;
        }
        argc=1;
        argv[0]=strdup("SeismicPlotWidget");
        XtSetLanguageProc(NULL, NULL, NULL);
        XtToolkitThreadInitialize();
        toplevel = XtVaAppInitialize(&AppContext, (char *)"seismicplot",NULL,0,
                                &argc,argv, NULL,NULL);
        main_w=XmCreateForm(toplevel,(char *) "seismicplot",NULL,0);
        EventLoopIsActive=false;
        Arg args[5];
        int n=0;
        XmString str;
        const string button_string("Continue");
        str=XmStringCreateLocalized(const_cast<char *>(button_string.c_str()));
        n=0;
        XtSetArg (args[n], XmNlabelString, str); n++;
        XtSetArg (args[n], XmNmultiClick, XmMULTICLICK_DISCARD); n++;
        continue_button=XmCreatePushButton(main_w,
            const_cast<char *>(button_string.c_str()),args,n);
        XmStringFree(str);
        XtAddCallback(continue_button,XmNactivateCallback,continue_callback,
                (XtPointer)this);
        XtManageChild(continue_button);
        XtManageChild(main_w);
        XtRealizeWidget(toplevel);
    } catch(...) {throw;}
}
void SeismicPlot::plot(TimeSeriesEnsemble& d)
{
    try{
        if(rshell[0]==NULL)
        {
            int n;
            Arg args[3];
            rshell[0]=XtVaCreatePopupShell("Component 0",
                    topLevelShellWidgetClass,toplevel,
                    XmNtitle,"Component 0 Plot",
                    XmNallowShellResize,True,
                    XmNwidth,1000,
                    XmNheight,800,
                    XmNdeleteResponse,XmUNMAP,NULL);
            n=0;
            XtSetArg(args[n],(char *) ExmNdisplayOnly,1); n++;
            seisw[0]=ExmCreateSeisw(rshell[0],(char *)"Seisw0",args,n);
            XtManageChild(seisw[0]);
            XtPopup(rshell[0],XtGrabNone);
        }
	XtAppLock(AppContext);
        XtVaSetValues(seisw[0],ExmNseiswEnsemble, (XtPointer) (&d),
            ExmNseiswMetadata,(XtPointer)(dynamic_cast<Metadata*>(this)), NULL);
	XtAppUnlock(AppContext);
        if(!EventLoopIsActive) this->launch_Xevent_thread_handler();
    }
    catch(...){throw;}
}
void SeismicPlot::refresh()
{
    try {
        int i;
	XtAppLock(AppContext);
        for(i=0;i<3;++i)
        {
            if(seisw[i]!=NULL) XtVaSetValues(seisw[i],ExmNseiswMetadata, 
                    (XtPointer)(dynamic_cast<Metadata*>(this)), NULL);
        }
        XtAppUnlock(AppContext);
    }catch(...) {throw;}
}
/* TimeSeries and ThreeComponentSeismograms are plotted in window 1
   using the same method as above.  Do this by simply creating a temporary
   ensemble object and calling that method */
void SeismicPlot::plot(TimeSeries& d)
{
    try {
        TimeSeriesEnsemble e;
        e.member.push_back(d);
        this->plot(e);
    } catch(...){throw;}
}
void SeismicPlot::plot(ThreeComponentSeismogram& d)
{
    try {
        TimeSeriesEnsemble e;
        TimeSeries *comp;
        int i;
        for(i=0;i<3;++i)
        {
            comp=ExtractComponent(d,i);
            e.member.push_back(*comp);
            delete comp;
        }
        this->plot(e);
    }catch(...){throw;}
}
void SeismicPlot::plot(ThreeComponentEnsemble& d)
{
    try {
        int k;
        for(k=0;k<3;++k)
        {
          if(rshell[k]==NULL)
          {
            int n;
            Arg args[3];
            char title[30];
            sprintf(title,"Component %d Plot",k);
            rshell[k]=XtVaCreatePopupShell(title,
                    topLevelShellWidgetClass,toplevel,
                    XmNtitle,title,
                    XmNallowShellResize,True,
                    XmNwidth,1000,
                    XmNheight,800,
                    XmNdeleteResponse,XmUNMAP,NULL);
            n=0;
            XtSetArg(args[n],(char *) ExmNdisplayOnly,1); n++;
            seisw[k]=ExmCreateSeisw(rshell[k],(char *)"Seisw0",args,n);
            XtManageChild(seisw[k]);
            XtPopup(rshell[k],XtGrabNone);
          }
        }
	XtAppLock(AppContext);
        auto_ptr<TimeSeriesEnsemble> comp[3];
        for(k=0;k<3;++k)
        {
            comp[k]=ExtractComponent(d,k);
            XtVaSetValues(seisw[k],ExmNseiswEnsemble, (XtPointer) (comp[k].get()),
                    ExmNseiswMetadata,(XtPointer)(dynamic_cast<Metadata*>(this)), NULL);
        }

	XtAppUnlock(AppContext);
        if(!EventLoopIsActive) this->launch_Xevent_thread_handler();
    }catch(...){throw;}
}
void EventHandler(SeismicPlot *plot_handle)
{
    XtAppLock(plot_handle->AppContext);
        do {
            XEvent event;
            XtAppNextEvent(plot_handle->AppContext,&event);
            XtDispatchEvent(&event);
        }while (plot_handle->WindowIsActive());
    XtAppUnlock(plot_handle->AppContext);
}
void SeismicPlot::launch_Xevent_thread_handler()
{
    EventLoopIsActive=true;  // This may need a mutex
    boost::thread Xevthrd(boost::bind(&EventHandler,this));
    Xevthrd.join();
    EventLoopIsActive=false;
}
