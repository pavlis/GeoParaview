#include <float.h>
#include "dmatrix.h"
#include "GCLMVFSmoother.h"
GCLMVFSmoother::GCLMVFSmoother() : firfilter()
{
    nrow=0;
    ncol=0;
    i0=0;
    j0=0;
}
GCLMVFSmoother::GCLMVFSmoother(dmatrix& fircoef, int i0in, int j0in)
         : firfilter(fircoef)
{
    nrow=firfilter.rows();
    ncol=firfilter.columns();
    i0=i0in;
    j0=j0in;
}
GCLMaskedScalarField GCLMVFSmoother::apply(GCLMaskedScalarField& parent)
{
    try{
        GCLMaskedScalarField field(parent);
        int i,j;
        int ii,jj;
        int is,js;  // variable start for averaging
        int ie,je;  // variable end
        int ifilt,jfilt;  // index in firfilter
        /* to handle variable hits we have to normalize each point
           independently.  This variable holds that sum for each point */
        double sumwt;  
        double smoothval;
        for(i=0;i<field.n1;++i)
            for(j=0;j<field.n2;++j)
            {
                if(field.point_is_valid(i,j))
                {
                    is=i-i0;
                    if(is<0) is=0;
                    js=j-j0;
                    if(js<0) js=0;
                    ie=i+nrow;
                    if(ie>field.n2) ie=field.n1;
                    je=j+ncol;
                    if(je>field.n2) je=field.n2;
                    sumwt=0.0;
                    smoothval=0.0;
                    /*loops intentionally focus on ii and jj as 
                      this is safer*/
                    for(ii=is,ifilt=0;ii<ie;++i,++ifilt)
                        for(jj=js,jfilt=0;jj<je;++jj,++jfilt)
                        {
                            if(field.point_is_valid(ii,jj))
                            {
                                sumwt+=firfilter(ifilt,jfilt);
                                smoothval+=firfilter(ifilt,jfilt)
                                            *parent.val[ii][jj];
                            }
                        }
                    /*Careful of null values.  Possible if user gives 0
                    coefficients in firfilter when smoother interacts 
                    with irregular shapes */
                    if(fabs(sumwt)<FLT_EPSILON)
                        field.val[i][j]=GCLFieldNullValue;
                    else
                        field.val[i][j]=smoothval/sumwt;
                }
            }
        return(field);
    }catch(...){throw;};
}
