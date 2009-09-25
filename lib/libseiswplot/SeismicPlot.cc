#include <X11/Xthreads.h>
#include "SeismicPlot.h"
using namespace std;
using namespace SEISPP;
#define APP_CLASS "seismicplot"
/* The unparameterized constructor uses a rather odd logic in that it is 
   actually just an interface to the pf constructor with a default pf. 
   This is odd because the default constructor usually has a rather different
   meaning, but to me, at least, makes sense in this context. */
SeismicPlot::SeismicPlot()
{
    Pf *pf;
    const string pfglobal_name("SeismicPlot");
    if(pfread(const_cast<char *>(pfglobal_name.c_str()),&pf))
        throw SeisppError("SeismicPlot constructor:  pfread failed for "
                + pfglobal_name);
    try {
        /* Normally main initializes the plot system in Motif.  I'm not
           sure this is essential, but we fake this hear with an interal
           argc and argv.  We keep them in the object as well as private 
           members as I am not sure of the memory status of these*/
        argc=1;
        argv[0]=strdup("SeismicPlotWidget");
        /* Old- does not work 
        XtSetLanguageProc(NULL, NULL, NULL);
        shell = XtVaAppInitialize(&AppContext, APP_CLASS,NULL,0,
                &argc,argv, NULL,NULL);
        form=XmCreateForm(shell,(char *) APP_CLASS,NULL,0);
        Arg args[3];
        int nargs=0;
        XtSetArg(args[nargs],(char *) ExmNzoomFactor,100); nargs++;
        XtSetArg(args[nargs],(char *) XmNpaneMaximum,5000); nargs++;
        XtSetArg(args[nargs],(char *) XmNpaneMinimum,500); nargs++;
        seisw[0]=ExmCreateSeisw(form,(char *)"Seisw",args,nargs);
        XtManageChild(seisw[0]);
        TimeSeriesEnsemble tsestartup(1,100);
        XtVaSetValues(seisw[0],ExmNseiswMetadata, 
                ExmNseiswEnsemble,(XtPointer)(&tsestartup),
                (XtPointer)dynamic_cast<Metadata*>(this),NULL);
        nwindows=1;
        XtRealizeWidget(shell);
        */
    } catch (...) {throw;}
}
/* This crazy function is required as a workaround for a bug in the present
   version of the seisw widget. It appears to fail on startup if it does not 
   have at least an empty ensemble to to plot.  This plots 100 straight line
   traces */
auto_ptr<TimeSeriesEnsemble> BuildDummyData()
{
    auto_ptr<TimeSeriesEnsemble> result(new TimeSeriesEnsemble(10,100));
    int i,j;
    for(i=0;i<10;++i)
    {
        /* Set required data to be valid */
        result->member[i].dt=1.0;
        result->member[i].t0=0.0;
        result->member[i].live=true;
        result->member[i].ns=100;
        result->member[i].tref=relative;
        result->member[i].put("samprate",1.0);
        result->member[i].put("time",0.0);
        result->member[i].put("nsamp",100);
        result->member[i].s.clear();
        for(j=0;j<100;++j) result->member[i].s.push_back(sin(2.0*M_PI*(double)(j)/20));
    }
    return(result);
}
SeismicPlot::SeismicPlot(Metadata& md) : Metadata(md)
{
    try {
        /* This is a complete duplicate of above.  There is probably 
           a tricky way to handle this some other way, but for now we
           do it this way.  Be careful if changes are made to one to 
           change the other too. */
        argc=1;
        argv[0]=strdup("SeismicPlotWidget");
        /* Old - does not work 
        XtSetLanguageProc(NULL, NULL, NULL);
        if(!XtToolkitThreadInitialize())
        {
            cerr << "Warning:  XtToolkitThreadInitialize failed on this platform."
                <<endl
                <<"Do not run an application with multiple, concurrent SeismicPlot objects."
                <<endl;
        }
        shell = XtVaAppInitialize(&AppContext, APP_CLASS,NULL,0,
                &argc,argv, NULL,NULL);
        form=XmCreateForm(shell,(char *) APP_CLASS,NULL,0);
        Arg args[3];
        int nargs=0;
        XtSetArg(args[nargs],(char *) ExmNzoomFactor,100); nargs++;
        XtSetArg(args[nargs],(char *) XmNpaneMaximum,5000); nargs++;
        XtSetArg(args[nargs],(char *) XmNpaneMinimum,500); nargs++;
        seisw[0]=ExmCreateSeisw(form,(char *)"Seisw",args,nargs);
        TimeSeriesEnsemble tsestartup(1,100);
        XtVaSetValues(seisw[0],ExmNseiswMetadata, 
                ExmNseiswEnsemble,(XtPointer)(&tsestartup),
                (XtPointer)dynamic_cast<Metadata*>(this),NULL);
        nwindows=1;
        XtRealizeWidget(shell);
        XtManageChild(seisw[0]);
        this->launch_Xevent_thread_handler();
        */
        XtSetLanguageProc(NULL, NULL, NULL);
        XtToolkitThreadInitialize();
        toplevel = XtVaAppInitialize(&AppContext, (char *)"seismicplot",NULL,0,
                                &argc,argv, NULL,NULL);
        main_w = XtVaCreateManagedWidget ((char *)"main window",
                xmMainWindowWidgetClass, toplevel,
                XmNscrollBarDisplayPolicy, XmAS_NEEDED,
                XmNscrollingPolicy,        XmAUTOMATIC,
                XmNwidth,1000,
                XmNheight,800,
                NULL);

        
        Arg args[4];
        int nargs=0;
        XtSetArg(args[nargs],(char *) ExmNzoomFactor,100); nargs++;
        XtSetArg(args[nargs],(char *) ExmNdisplayOnly,1); nargs++;
        seisw[0]=ExmCreateSeisw(main_w,(char *)"Seisw",args,nargs);
        auto_ptr<TimeSeriesEnsemble>x1=BuildDummyData();
        TimeSeriesEnsemble *tse=x1.get();
        XtVaSetValues(seisw[0],
                ExmNseiswMetadata,(XtPointer)(dynamic_cast<Metadata*>(this)),
                ExmNseiswEnsemble,(XtPointer)(tse),NULL);
        XtRealizeWidget(toplevel);
        //XtAppMainLoop(AppContext);
        this->launch_Xevent_thread_handler();
    } catch(...) {throw;}
}
void SeismicPlot::plot(TimeSeriesEnsemble& d)
{
    try{
        XtVaSetValues(seisw[0],ExmNseiswEnsemble, (XtPointer) (&d),
            ExmNseiswMetadata,(XtPointer)(dynamic_cast<Metadata*>(this)), NULL);
    
    }
    catch(...){throw;}
}
void SeismicPlot::refresh()
{
    try {
        int i;
        for(i=0;i<nwindows;++i)
        {
            XtVaSetValues(seisw[i],ExmNseiswMetadata, 
                    (XtPointer)(dynamic_cast<Metadata*>(this)), NULL);
        }
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
/* This is the hard one.  It does a kind of weird thing that may not survive
   a test of time.  It creates two popup shells and sets nwindow to 3 when
   first called.  Otherwise it splits the ensemble into three components
   and plots each in a different window */
void SeismicPlot::plot(ThreeComponentEnsemble& d)
{
    try {
        stringstream serr;
        serr << "ThreeComponentEnsemble plot method:  ";
        int k;
        if(nwindows==1)
        {
            rshell[0]=XtVaCreatePopupShell("Component 2",
                    topLevelShellWidgetClass, toplevel,
                    XmNtitle,"ThreeComponent Data - Component 2",
                    XmNallowShellResize,True,XmNwidth,800,XmNheight,800,
                    XmNdeleteResponse,XmUNMAP,NULL);
        seisw[1] = XtVaCreateManagedWidget ((char *)"mainw_c2",
                xmMainWindowWidgetClass, rshell[0],
                XmNscrollBarDisplayPolicy, XmAS_NEEDED,
                XmNscrollingPolicy,        XmAUTOMATIC,
                XmNwidth,1000,
                XmNheight,800,
                NULL);
            XtManageChild(seisw[1]);
            ++nwindows;
            XtPopup(rshell[0],XtGrabNone);
            rshell[1]=XtVaCreatePopupShell("Component 3",
                    topLevelShellWidgetClass, toplevel,
                    XmNtitle,"ThreeComponent Data - Component 3",
                    XmNallowShellResize,True,XmNwidth,800,XmNheight,800,
                    XmNdeleteResponse,XmUNMAP,NULL);
        seisw[2] = XtVaCreateManagedWidget ((char *)"mainw_c3",
                xmMainWindowWidgetClass, rshell[0],
                XmNscrollBarDisplayPolicy, XmAS_NEEDED,
                XmNscrollingPolicy,        XmAUTOMATIC,
                XmNwidth,1000,
                XmNheight,800,
                NULL);
            ++nwindows;
            XtManageChild(seisw[2]);
            XtPopup(rshell[1],XtGrabNone);
        }
        if(nwindows==3)
        {
            auto_ptr<TimeSeriesEnsemble> comp[3];
            for(k=0;k<3;++k)
            {
                comp[k]=ExtractComponent(d,k);
                /* this requires a fix eventually.  As I'm writing this
                   I know this won't work but without access to the web
                   I don't know how to convert an auto_ptr into a raw
                   pointer required by this routine.  This will also 
                    leak memory badly
                 */
                TimeSeriesEnsemble *rawptr=new TimeSeriesEnsemble(*comp[k]);
                XtVaSetValues(seisw[k],ExmNseiswEnsemble, 
                        (XtPointer) rawptr,NULL);
            }
        }
        else
        {
            serr << "nwindows private variable ="<<nwindows
                << " which is an illegal value"<<endl
                <<"Coding error.  Program has probably overwritten itself."
                <<"  Contact author"<<endl;
            throw SeisppError(serr.str());
        }
    }catch(...){throw;}
}
void SeismicPlot::launch_Xevent_thread_handler()
{
    cout << "Entering launch_Xevent_thread_handler"<<endl;
    if(pthread_create((pthread_t *)(&event_handler_thread_id), NULL,
                Xevent_loop,(void *)(&AppContext)))
    {
        throw SeisppError(string("pthread_create failed trying to create Xevent handler thread."));
    }
}

double SeismicPlot::PickTime()
{
    cerr << "PickTime method has not been implemented"<<endl;
}
TimeWindow SeismicPlot::PickWindow()
{
    cerr << "PickWindow method has not been implemented"<<endl;
}

