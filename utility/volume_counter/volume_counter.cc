#include <string>
#include <fstream>
#include <list>
#include <vector>
#include "gclgrid.h"
#include "seispp.h"
#include "dbpp.h"
#include "PfStyleMetadata.h"
#include "EventCatalog.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "volume_counter db out [-ss subset_condition -pf pffile -v] < infile"<<endl
        << "   db is Antelope db  (needs event->origin and will write gcl tables"<<endl
        << "   out is name for output gclfield data (file name and in gclfield table)"<<endl
        << "   (note out is vectorfield with 0 count and 1 area/volume ration estimate)" <<endl
        << "   -ss defines subset condition applied before processing (can be important for efficiency)"<<endl
        << "   -pf for alternate pf (check pf before running of course)"<<endl
        << "   -v to turn on verbose error reporting"<<endl;;
    exit(-1);
}
class DistanceTest : public unary_function<Hypocenter,bool>
{
    public:
        DistanceTest(GCLgrid3d *gclg,int i, int j, int k,double rc)
        {
            g=gclg;
            x1p=gclg->x1[i][j][k];
            x2p=gclg->x2[i][j][k];
            x3p=gclg->x3[i][j][k];
            rcutoff=rc;
        };
        bool operator() (const Hypocenter& h) const
        {
            double ssq,delta;
            Cartesian_point cp=g->gtoc(h.lat,h.lon,r0_ellipse(h.lat)-h.z);
            delta=x1p-cp.x1;
            ssq=delta*delta;
            delta=x2p-cp.x2;
            ssq+=delta*delta;
            delta=x3p-cp.x3;
            ssq+=delta*delta;
            if(sqrt(ssq)<=rcutoff) return true;
            return false;
        };
    private:
        double x1p,x2p,x3p;
        /* This is the cutoff distance */
        double rcutoff;
        /* A pointer for efficiency.   Would be a horrible
           overhead otherwise. */
        GCLgrid3d *g;
};

double sum_areas(auto_ptr<EventCatalog> events,string magname,double nullvalue)
{
    try {
        double sumresult(0.0);
        events->rewind();
        int nevents=events->size();
        int i;
        /* EventCatalog has an internal counter - odd increment operator 
           is a syntatic oddity of usign an auto_ptr */
        for(i=0;i<nevents;++i,++(*events))
        {
            Metadata aux=events->current_aux();
            double mag=aux.get_double(magname);
            /* assume anything negative is null */
            if(mag<0.0) mag=nullvalue;
            /* glp see momentscalingtheory.pdf for derivation of this formula*/
            double area=pow(10.0,(mag-3.296));
            sumresult += area;
        }
        return sumresult;
    }catch(...){throw;};
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<3) usage();
    string dbname(argv[1]);
    string outfield(argv[2]);
    string pffile("volume_counter");
    bool apply_subset(false);
    string subset_string;
    int i,j,k;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pffile=string(argv[i]);
        }
        else if(sarg=="-ss")
        {
            ++i;
            if(i>=argc) usage();
            apply_subset=true;
            subset_string=sarg;
        }
        else if(sarg=="-v")
            SEISPP_verbose=true;
        else
            usage();
    }
    pffile=pffile+".pf";
    ios::sync_with_stdio();
    try{
        DatascopeHandle dbh(dbname,false);
        if(apply_subset) 
        {
            if(SEISPP_verbose)
                cout << "Original db size before subset application"
                    << dbh.number_tuples()<<endl
                    << "Applying this subset condition:  "
                    << apply_subset<<endl;
            dbh.subset(subset_string);
            if(SEISPP_verbose)
                cout << "Number of tuples after subset = "
                    << dbh.number_tuples()<<endl;
        }
        PfStyleMetadata control=SEISPP::pfread(pffile);
        if(SEISPP_verbose)
        {
            cout << "Parameter file values read"<<endl;
            control.pfwrite(cout);
        }
        double search_radius=control.get_double("search_radius");   // this is in km
        string dir=control.get_string("output_directory");
        /* This is essential the same as makegclgrid but using the new PfStyleMetadata object */
        double lat0,lon0,r0,azimuth_x,dx1,dx2,dx3;
        int n1,n2,n3,i0,j0;
        lat0=control.get_double("origin_latitude");
        lon0=control.get_double("origin_longitude");
        lat0=rad(lat0);
        lon0=rad(lon0);
        r0=r0_ellipse(lat0);
        azimuth_x=control.get_double("x_axis_azimuth");
        double azimuth_y=azimuth_x-90.0;
        azimuth_y=rad(azimuth_y);
        dx1=control.get_double("delta_x1");
        dx2=control.get_double("delta_x2");
        dx3=control.get_double("delta_x3");
        n1=control.get_int("n1");
        n2=control.get_int("n2");
        n3=control.get_int("n3");
        i0=control.get_int("origin_offset_x1");
        j0=control.get_int("origin_offset_x2");
        string gridname=control.get_string("gridname");
        bool use_netmag=control.get_bool("use_netmag");  
        string magnitude_attribute=control.get_string("magnitude_attribute");
        double null_magnitude_default=control.get_double("null_magnitude_default");
        if(use_netmag)
        {
            if(SEISPP_verbose)
                cout << "Attempting to join netmag table"<<endl
                    << "Current view size before join="<<dbh.number_tuples()
                    <<endl;;
            dbh.natural_join(string("netmag"));
            if(SEISPP_verbose)
                cout << "Number of rows after joining netmag="
                    << dbh.number_tuples()<<endl;
        }
        /* We always load ms, mb, and ml as they are part of origin*/
        MetadataList maglist;
        Metadata_typedef md;
        md.tag=string("origin.ms");
        md.mdt=MDreal;
        maglist.push_back(md);
        md.tag=string("origin.mb");
        maglist.push_back(md);
        md.tag=string("origin.ml");
        maglist.push_back(md);
        if(use_netmag)
        {
            md.tag=string("netmag.magnitude");
            maglist.push_back(md);
        }
        AttributeMap am("css3.0");
        EventCatalog eqs(dbh,maglist,am);
        /* Use new here so we can save some memory once
           we create the field object */
        GCLgrid3d *g=new GCLgrid3d(n1,n2,n3,gridname,lat0,lon0,r0,
                azimuth_y,dx1,dx2,dx3,i0,j0);
        GCLvectorfield3d f(*g,2);
        /* this is not essential, but best to be explicit */
        f.zero();
        delete g;

        /* Now loop through all the points and count */
        const double volume(4.0*M_PI*search_radius*search_radius*search_radius/3.0 );
        list<Cartesian_point>::iterator points;
        for(i=0;i<f.n1;++i)
            for(j=0;j<f.n2;++j)
            {
                if(SEISPP_verbose)
                        cout<<"Starting on index i,j="<<i<<" "<<j<<endl
                            << "Time is "<<strtime(now())<<endl;
                for(k=0;k<f.n3;++k)
                {
                    DistanceTest closeones(&f,i,j,k,search_radius);
                    auto_ptr<EventCatalog> locals=eqs.subset<DistanceTest>(closeones);
                    int count=locals->size();
                    double asum=sum_areas(locals,magnitude_attribute,null_magnitude_default);
                    /* Note we convert to count per unit volume */
                    f.val[i][j][k][0]=static_cast<double>(count)/volume;
                    f.val[i][j][k][1]=asum/volume;
                }
            }
        /* This always writes the field output data to a file set by the name outfield */
        f.save(dbh,dir,dir,outfield,outfield);
    }catch(std::exception& err)
    {
        cerr << "Something threw the following error:"<<endl
            << err.what()<<endl;
    }
    catch(SeisppError& serr)
    {
        serr.log_error();
    }
}






