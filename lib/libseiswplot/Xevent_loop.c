#include <X11/Intrinsic.h>
void *Xevent_Loop(void *arg);
void *Xevent_loop(void *arg)
{
        XtAppContext *appc;
        appc=(XtAppContext *)(arg);
        XtAppMainLoop(*appc);
        return(NULL);
}
