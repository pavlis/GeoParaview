#include "seispp.h"
#include "gclgrid.h"
#include "GCLMasked.h"
#include "MultiMask.h"
using namespace SEISPP;
enum infile_type{GRID, SCALARFIELD, VECTORFIELD};
void usage()
{
    cerr << "build_masked_surface (-g indata | -sf indata | -vf indata ) -o outdata [-maxval xxx -minval xxx]< list_of_polygons"<<endl
        << "Use -g for GCLgrid only, -sf for GCLscalarfield, and "
        << "-vf for GCLvectorfield"<<endl
        << "Use -o to specify output masked 2d field object "
        <<"(default MaskedGrid)"<<endl
        <<"-maxval and -minval are a ceiling and floor for field values"<<endl
        <<"(cells with values outside minval to maxval range will be masked"
        <<" - this feature is off by default)"
        <<endl;

    exit(-1);
}
/* This series of procedures implement a masked defined by points exceeding
   a min and max range.  More or less like a clip by value function.
   All return the number of points masked by that method. */
int mask_range(GCLMaskedScalarField *f,double vmin, double vmax)
{
    int count(0);
    int i,j;
    for(i=0;i<f->n1;++i)
        for(j=0;j<f->n2;++j)
        {
            double vtest=f->val[i][j];
            if( (vtest<vmin) || (vtest>vmax) )
            {
                f->mask_point(i,j);
                ++count;
            }
        }
    return(count);
}
/* For a vector field we use the L2 norm of components for size */
int mask_range(GCLMaskedVectorField *f,double vmin, double vmax)
{
    int count(0);
    int i,j,k;
    for(i=0;i<f->n1;++i)
        for(j=0;j<f->n2;++j)
        {
            double vtest,vtmp;
            for(k=0,vtest=0.0;k<f->nv;++k) 
            {
                vtmp=f->val[i][j][k];
                vtest += vtmp*vtmp;
            }
            vtest=sqrt(vtest);
            if( (vtest<vmin) || (vtest>vmax) )
            {
                f->mask_point(i,j);
                ++count;
            }
        }
    return(count);
}
/* For a plain grid we use the depth method */
int mask_range(GCLMaskedGrid *g,double vmin, double vmax)
{
    int count(0);
    int i,j;
    for(i=0;i<g->n1;++i)
        for(j=0;j<g->n2;++j)
        {
            double vtest=g->depth(i,j);
            if( (vtest<vmin) || (vtest>vmax) )
            {
                g->mask_point(i,j);
                ++count;
            }
        }
    return(count);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc<5) usage();
    string infile,outfile;
    infile_type intype;
    bool range_masking(false);
    /* Arbitrary but avoids need to test for both min and max in
       arg list */
    double minval(-0.99e99),maxval(0.99e99);
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
        else if(sarg=="-minval")
        {
            ++i;
            if(i>=argc)usage();
            range_masking=true;
            minval=atof(argv[i]);
        }
        else if(sarg=="-maxval")
        {
            ++i;
            if(i>=argc)usage();
            range_masking=true;
            maxval=atof(argv[i]);
        }
        else
            usage();
    }
    try{
        if(range_masking)
        {
            cout << "Range masking has been enabled"<<endl
                << "Field values smaller than "<<minval
                <<" or greater than "<<maxval<<" will be masked"
                <<endl
                << "Note:  for GCLgrid objects computed depths"
                << " will be tested this way for masking"<<endl;
        }
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
        /* there might be a way to use a polymophic construct here, but
           I couldn't figure it out. */
        if(range_masking)
        {
          int nmasked;
          switch(intype)
          {
            case SCALARFIELD:
                nmasked=mask_range(msf,minval,maxval);
                break;
            case VECTORFIELD:
                nmasked=mask_range(mvf,minval,maxval);
                break;
            case GRID:
            default:
                nmasked=mask_range(mg,minval,maxval);
          }
          cout << "Range masking set the mask on "<<nmasked<<" grid points"
              <<endl;
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

