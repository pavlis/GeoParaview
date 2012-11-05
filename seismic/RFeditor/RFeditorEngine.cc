#include <set>
#include <list>
#include <algorithm>
#include "stack.h"
#include "RFeditorEngine.h"
using namespace std;
using namespace SEISPP;
const int StackThreshold(5);   // do not stack unless count above this
/* Helpers for edit method. */
vector<int> stack_and_sort(TimeSeriesEnsemble& d)
{
    try{
        vector<int> result;
        /* Silently do nothing if the ensemble is empty */
        if(d.member.size()<=0) return result;
        /* This works for RF data from EARS that all have a common
           length.  DO NOT transport this code without making this
           more general */
        TimeWindow twin(d.member[0].t0,d.member[0].endtime());
        /* This is a convenient way to post original order to the 
           ensemble.  Used before returning to produce result vector */
        vector<TimeSeries>::iterator dptr;
        int i;
        for(i=0,dptr=d.member.begin();dptr!=d.member.end();++dptr,++i)
            dptr->put("member_number",i);
        /* Always use the robust method with the stack and robust
           window the same */
        Stack s(d,twin,twin,RobustSNR);
        /* Frozen name a problem, but label this stack so clear it
           is such */
        s.stack.put(string("sta"),string("stack"));
        /* Sort always by stack weight */
        sort(d.member.begin(),d.member.end(),less_stackwt<TimeSeries>());
        result.reserve(d.member.size());
        for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        {
            int im;
            im=dptr->get_int("member_number");
            result.push_back(im);
        }
        /* push two copies of the stack to the front of the vector
           container.  Mark one dead to serve as a spacer.  This will
           make stack appear at the bottom of the plot with a gap.  This 
           is a bit dangerous and could scramble some things later so
           caution if the algorithm is later modified. */
        d.member.insert(d.member.begin(),s.stack);
        d.member.insert(d.member.begin(),s.stack);
        d.member[1].live=false;
        return result;
    } catch(...){throw;};
}
/* Returns a new ensemble reordered by ordering.   Does this 
   blindly without checking bounds for efficiency as it assumes
   this is intimately linked to the procedure immediately above.
   Do not transport blindly.  Also not very efficient because of
   double copy operation*/
TimeSeriesEnsemble reorder_transverse(TimeSeriesEnsemble& t,
        vector<int> ordering)
{
    try{
        TimeSeriesEnsemble result(t);
        int i;
        /* This is not the most obvious way to do this, but
           I think this will be more efficient than indexing
           the other way */
        result.member.clear();
        int nm=t.member.size();
        result.member.reserve(nm);
        for(i=0;i<nm;++i)
        {
            int im=ordering[i];
            result.member.push_back(t.member[im]);
        }
        return(result);
    } catch(...){throw;};
}


RFeditorEngine::RFeditorEngine(Metadata& params) 
    : Rwindow(params), Twindow(params)
{
    try {
        string edit_mode;
        edit_mode=params.get_string("editing_mode");
        if(edit_mode=="cutoff")
        {
            cutoff_editing=true;
            Rwindow.toggle_edit_mode();
            Twindow.toggle_edit_mode();
        }
        else
            cutoff_editing=false;
    }catch(...){throw;};
}
/* Helpers for edit method.  */ 
/* The first (apply_edits) takes index numbers returned by the 
   editing interface and works through the ensemble tse using
   those indices.  It forces a kill of those member, although 
   that may be redundant as the widget does this too I think.
   More importantly it extracts the event id (keyed by evidkey)
   of each seismogram killed. A set container of these evid
   values is returned. */
set<long>  apply_kills(TimeSeriesEnsemble& tse, set<int> killset)
{
    try {
        set<int>::iterator iptr;
        set<long> evids_to_kill;
        for(iptr=killset.begin();iptr!=killset.end();++iptr)
        {
            /* Skip stack data */
            if(tse.member[*iptr].get_string("sta")=="stack") continue;
            tse.member[*iptr].live=false;
            long evid=tse.member[*iptr].get_long(evidkey);
            evids_to_kill.insert(evid);
            cout << "Killing evid "<<evid<<endl;
        }
        return(evids_to_kill);
    }catch(...){throw;};
}
/* Followup companion helper to apply_kills.  This procedure takes the
set of event ids returned by apply_edits and uses them to kill any
members of tse passed to this function with matching evid.  
In this code this is called to the complement of radial or transverse
depending on which was being edited. */
void apply_kills_to_other(TimeSeriesEnsemble& tse, set<long> evids_to_kill)
{
    try {
        vector<TimeSeries>::iterator tptr;
        for(tptr=tse.member.begin();tptr!=tse.member.end();++tptr)
        {
            // Skip stack trace when present 
            string sta=tptr->get_string("sta");
            if(sta=="stack") continue;
            int evid=tptr->get_long(evidkey);
            if(evids_to_kill.find(evid)!=evids_to_kill.end())
            {
                cout << "Transverse killing evid "<<evid<<endl;
                tptr->live=false;
            }
        }
    }catch(...){throw;};
}
int RFeditorEngine::edit(TimeSeriesEnsemble& rd, TimeSeriesEnsemble& td)
{
    try{
        /* This procedure runs the robust stacker and sorts the radial
           data by stack weight.   Returns a vector of ints linking 
           original member numbers to sorted member numbers.  We use
           this is to reorder transverse data before plotting. */
        if(rd.member.size()>StackThreshold)
        {
            vector<int> reordering=stack_and_sort(rd);
            td=reorder_transverse(td,reordering);
        }
        /* This changes the plot titles for each window.  Assumes sta
           is defined in ensemble metadata */
        string sta=rd.get_string("sta");
        string title=string("Transverse RF data for station ")+sta;
        Twindow.put("title",title);
        title=string("Radial RF data for station ")+sta;
        Rwindow.put("title",title);
        //Rwindow.refresh();
        /* This displays the transverse but immediately exits the event loop */
        Twindow.plot(td,false);
        Twindow.refresh();
        /* Plot the radial data in edit mode */
        Rwindow.plot(rd,true);
        Twindow.ExitDisplay(); // kill the event thread or it will multiply
        set<int> rkills;
        rkills=Rwindow.report_kills();
        set<long> events_to_kill;
        if(rkills.size()>0)
        {
            events_to_kill=apply_kills(rd,rkills);
            apply_kills_to_other(td,events_to_kill);
        }
        set<int> tkills;
        Twindow.plot(td,true);
        tkills=Twindow.report_kills();
        if(tkills.size()>0)
        {
            events_to_kill=apply_kills(td,tkills);
            apply_kills_to_other(rd,events_to_kill);
        }
    }catch(...){throw;};
}
