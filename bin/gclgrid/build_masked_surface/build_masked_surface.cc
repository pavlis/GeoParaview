#include "seispp.h"
#include "gclgrid.h"
#include "GCLMasked.h"
#include "MultiMask.h"
using namespace SEISPP;
enum infile_type{GRID, SCALARFIELD, VECTORFIELD};
void usage()
{
    cerr << "build_masked_surface (-g indata | -sf indata | -vf indata ) -o outdata < list_of_polygons"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc<5) usage();
    string infile,outfile;
    infile_type intype;
    int i,j;
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-g")
        {
            ++i;
            if(i>=argc)usage();
            infile=string(argv[i]);
            intype=GRID;
        }
        else if(sarg=="-sf")
        {
            ++i;
            if(i>=argc)usage();
            infile=string(argv[i]);
            intype=SCALARFIELD;
        }
        else if(sarg=="-vf")
        {
            ++i;
            if(i>=argc)usage();
            infile=string(argv[i]);
            intype=VECTORFIELD;
        }
        else if(sarg=="-o")
        {
            ++i;
            if(i>=argc)usage();
            outfile=string(argv[i]);
        }
        else
            usage();
    }
    try{
        MultiMask polymask(cin);
        /* We cast all types to this base to exploit polymorphism */
        GCLgrid *g;
        GCLscalarfield *sf;
        GCLvectorfield *vf;
        switch(intype)
        {
            case GRID:
                g=new GCLgrid(infile);
                break;
            case SCALARFIELD:
                sf=new GCLscalarfield(infile);
                g=dynamic_cast<GCLgrid*>(sf);
                break;
            case VECTORFIELD:
                vf=new GCLvectorfield(infile);
                g=dynamic_cast<GCLgrid*>(vf);
                break;
            default:
                cerr << "Coding error:  input type has not been set correctly"<<endl;
                usage();
        }
        /* This will hold the mask we will apply.  We clone it's geometry from g to guarantee
           consistency */
        GCLMask mask(*g);
        /* now compute the mask */
        for(i=0;i<g->n1;++i)
            for(j=0;j<g->n2;++j)
            {
                double lat,lon;
                lat=g->lat(i,j);
                lon=g->lon(i,j);
                if(polymask.is_inside(lat,lon))
                    mask.enable_point(i,j);
                else
                    mask.mask_point(i,j);
            }
        GCLMaskedGrid *mg;
        GCLMaskedScalarField *msf;
        GCLMaskedVectorField *mvf;
        switch(intype)
        {
            case SCALARFIELD:
                msf=new GCLMaskedScalarField(*sf,mask);
                break;
            case VECTORFIELD:
                mvf=new GCLMaskedVectorField(*vf,mask);
                break;
            case GRID:
            default:
                mg=new GCLMaskedGrid(*g,mask);
        }
        /* Finally save the data to a special file format for 
           masked GCLgrid objects */
        switch(intype)
        {
            case SCALARFIELD:
                msf->save(outfile);
                break;
            case VECTORFIELD:
                mvf->save(outfile);
                break;
            case GRID:
            default:
                mg->save(outfile);
        }
    }catch(std::exception& stex)
    {
        cerr << stex.what();
        exit(-1);
    }
}

