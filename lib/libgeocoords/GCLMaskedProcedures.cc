double CartesianDistance(Cartesian_point& cp1, Cartesian_point& cp2)
{
    double d,dx;
    dx=cp2.x1-cp1.x1;
    x+=dx*dx;
    dx=cp2.x2-cp1.x2;
    x+=dx*dx;
    dx=cp2.x3-cp1.x3;
    x+=dx*dx;
    return(sqrt(x));
}
// These procedures should probably set values outside mask to some
// stock not define value 
GCLMaskedVectorfield ComputeNormals(GCLmaskedgrid& raw)
{
    try {
        GCLMaskedVectorField result(raw,3);
        int i,j,k;
        double dx1[3],dx2[3],x1crossx2[3];
        for(i=0;i<(n1-1);++i)
            for(j=0;j<(n2-1);++i)
            {
                dx1[0]=raw.x1[i+1][j]-raw.x1[i][j];
                dx1[1]=raw.x2[i+2][j]-raw.x2[i][j];
                dx1[2]=raw.x3[i+3][j]-raw.x3[i][j];
                dx2[0]=raw.x1[i][j+1]-raw.x1[i][j];
                dx2[1]=raw.x2[i][j+1]-raw.x2[i][j];
                dx2[2]=raw.x3[i][j+1]-raw.x3[i][j];
                d3cros(dx1,dx3,x1crossx3);
                double nmag;
                nmag=dr3mag(x1crossx3); // to make unit normals 
                for(k=0;k<3;++k) result.val[i][j][k]=x1crossx2[k]/nmag;
            }
        return(result);
    } catch(...){throw;};
}
GCLMaskedVectorField ComputeLocalVerticals(GCLMaskedGrid& raw)
{
    try {
        GCLMaskedVectorField result(raw,3);
        int i,j,k;
        const double Roffset(20.0);  // radius offset to compute normal 
        Cartesian_point cp;
        Geographic_point gp;
        double rawvector[3];
        for(i=0;i<(n1-1);++i)
            for(j=0;j<(n2-1);++i)
            {
                gp=raw.ctog(raw.x1[i][j],raw.x2[i][j],raw.x3[i][j]);
                gp.r+=Roffset;
                cp=raw.gtoc(gp);
                rawvector[0]=cp.x1-raw.x1[i][j];
                rawvector[1]=cp.x2-raw.x2[i][j];
                rawvector[2]=cp.x3-raw.x3[i][j];
                double nrmraw;
                nrmraw=dr3mag(rawvector);
                for(k=0;k<3;++k) result.val[i][j][k]=rawvector[k]/nrmraw;
            }
        // fill outer border using last valid value
        j=n2-1;
        for(i=0;i<n1;++i)
            for(k=0;k<3;++k) result.val[i][j][k]=result.val[i][j-1][k];
        i=n1-1;
        for(j=0;j<n2;++j)
            for(k=0;k<3;++k) result.val[i][j][k]=result.val[i-1][j][k];
    } catch(...){throw;};
}
/* The procedurs below use a derivative based on a horizontal (constant r)
   delta from a point i,j in the grid.   This standardizes that 
   calculation in one place. 

Arguments:
g - grid to compute dx1 from
i,j - index into grid

A positive value is indicates a normal return. Test for a negative
value to detect errors in calculation.   Two ways this happens:
(1) trying to access i value outside grid range, and (2) when 
i or i+1 is in a masked (invalid) region.
   */
double dx1_offset(GCLMaskedGrid& g,int i, int j)
{
    double dx1;
    const double BadResult(-99999.9);  // arbitrary negative number
    if(i>=(g.n1-1)) return(BadResult);
    Cartesian_point cp1,cp2r1;
    Geographic_point gp1,gp2,gp2r1;

    if(g.point_is_valid(i,j) && g.point_is_valid(i+1,j))
    {
        cp1.x1=g.x1[i][j];
        cp1.x2=g.x2[i][j];
        cp1.x3=g.x3[i][j];
        gp1=g.ctog(g.x1[i][j],g.x2[i][j],
                g.x3[i][j]);
        gp2=g.ctog(g.x1[i+1][j],g.x2[i+1][j],
                g.x3[i+1][j]);
        gp2r1=gp2;
        gp2r1.r=gp1.r;
        cp2r1=g.gtoc(gp2r1);
        dx1=CartesianDistance(cp1,cp2r1);
    }
    else
        dx1=BadResult;
    return(dx1);
}
/* Companion to above. Returns the geographic points of the point
   i+1,j from i,j at constant r.   i.e. point is same lat,lon as
   i+1,j but r is copied form point at i,j.   Same as above
   procedure making this somewhat repetitious, but another common
   calculation here simplified this way.  This function 
   does NOT do any error checking assuming it is ONLY called 
 AFTER calling the previous procedure. */
Geographic_point constant_r_project(GCLMaskedGrid& g,int i, int j)
{
    Cartesian_point cp1,cp2r1;
    Geographic_point gp1,gp2,gp2r1;
    gp1=g.geo_coordinates(i,j);
    gp2=g.geo_coordinates(i+1,j);
    gp2.r=gp1.r;
    return(gp2);
}


/* This procedure computes directional derivative dr/dx1 at each 
   grid point by a forward difference formula.   j=n2-1 is a copy 
   of n2-2.  Implementation using simple forward difference scheme
   was chosen because original purpose  of this procedure is locked
   to slabmodel program where this is sensible because x2=0 is 
   normally a trench position which is a hard constraint.   High 
   accuracy is also not a huge issue.   

   Also this is based on an approximation that is good only when the 
   dx1 intervals are not huge.   that is dx1 increments are measured 
   with cartesian coordinates not a great circle path geometry.  
   Unlikely to be a huge issue. 
 
 raw is the grid used to generate the field.  mindip is the 
 minimum dip (in radians) allowed.  That is, if the computed
 dip along the x1 direction is less than this value it will be
 truncated.  To turn this off just set mindip to -M_PI_2*/
   
GCLMaskedScalarField compute_drdx1(GCLmaskedgrid& raw,double mindip)
{
    try{
        GCLMaskedScalarField result(raw);
        int i,j;
        Cartesian_point cp1,cp2r1;
        Geographic_point gp1,gp2,gp2r1;
        double dx1;
        //dx1/dr is a tangent.  For efficiency we run mindip compare
        // against the tangent
        double mindiptan;
        // This is needed to allow setting mindip to M_PI_2 to turn 
        // of the test.   Needed because tan(M_PI_2) is -infinity.
        if(mindip==M_PI_2)
        {
            mindiptan=-99999999.9;
        else
            mindiptan=tan(mindip);
        for(j=0;j<n2;++j)
            for(i=0;i<(n1-1);++i)
            {
                /* ASsumes this procedure returns a negative
                   value when one of the points is masked*/
                dx1=dx1_offset(result,i,j);
                if(dx1>0.0)
                {
                    gp1=result.geo_coordinates(i,j);
                    gp2=result.geo_coordinates(i+1,j);
                    /* notice this is normally negative */
                    result.val[i][j]=(gp2.r-gp1.r)/dx1;
                    if(-result.val[i][j]<mindiptan)
                        result.val[i][j]=mindiptan;
                }
                else
                {
                    // Not essential but better to be sure initialized
                    result.val[i][j]=GCLFieldNullValue;
                }
            }
        /* Fill i=n1-1 with copies of n1-2 */
        i=n1-1;
        for(j=0;j<n2;++j)
        {
            if(result.point_is_valid(i,j))
                result.val[i][j]=result.val[i-1][j];
            else
                result.val[i][j]=GCLFieldNullValue;
        }

    } catch(...){throw;};
}

GCLMaskedScalarField x1_integrator(GCLMaskedGrid& parent,
                GCLMaskedScalarField& drdx1field)
{
    try {
        double dx1,dr;
        Geographic_point gp0,newgp;
        Cartesian_point cp;
        GCLMaskedScalarField ifld(parent);
        int i,j;
        for(j=0;j<ifld.n2;++j)
            for(i=1;i<ifld.n1;++i)
            {
                /* if this first several points are masked this
                   will skip until both points are valid */
                if(ifld.point_is_valid(i-1,j))
                {
                    dx1=x1_offset(parent,i-1,j);
                    /* do nothing when i,j is masked */
                    if(dx1<0.0) continue;
                    gp0=ifld.geo_coordinates(i-1,j);
                    newgp=ifld.geo_coordinates(i,j);
                    dr=drdx1field.val[i-1][j]*dx1;
                    newgp.r=gp0+dr;
                    cp=ifld.gtoc(newgp);
                    ifld.x1[i][j]=cp.x1;
                    ifld.x2[i][j]=cp.x2;
                    ifld.x3[i][j]=cp.x3;
                }
            }
        return(ifld);
    } catch(...){throw;};
}
/* small helper to avoid many repetitions of this */
Cartesian_point get_point(GCLgrid& g, int i, int j)
{
    /* Intentionally do not test i,j assuming this is only used
       inside this file  where it should not happen */
    Cartesian_point cp;
    cp.x1=g.x1[i][j];
    cp.x2=g.x2[i][j];
    cp.x3=g.x3[i][j];
    return(cp);
}
GCLMaskedGrid LinearDepthCorrection(GCLMaskedGrid& raw, 
                vector<TiePoints>& tpv)
{
    try{
        GCLMaskedGrid result;
        if(tpv.size()!=raw.n2) 
        {
            ostringstream sserr;
            sserr << "LinearDepthCorrection procedure:  size mismatch"<<endl
                <<"Vector of time point size="<<tpv.size()<<endl
                <<"Grid dimension in x2 direction is "<<raw.nw<<endl
                <<"They must match (FATAL ERROR)"<<endl;
            throw GCLgridError(sserr.str());
        }
        Geographic_point gp,gptie;
        Cartesian_point cp1,cp2;
        /* This stores total arc lengths used as basis for linear 
           drift like correction to points to match each tie point */
        vector<double> arclengths;
        arclengths.reserve(raw.n1);
        int i,j;
        for(j=0;j<raw.n2;++j)
        {
            arclengths.clear();
            double arc,dr,drds;
            for(i=0,arc=-1.0;i<raw.n1;++i)
            {
                if(raw.point_is_valid(i,j))
                {
                    if(arc<0.0)
                    {
                        arclengths.push_back(0.0);
                        cp1=get_point(raw,i,j);
                    }
                    else
                    {
                        cp2=get_point(raw,i,j);
                        arc += CartesianDistance(cp1,cp2);
                        arclengths.push_back(arc);
                        cp1=cp2;
                    }
                }
                else
                {
                /*Push a negative number for any masked point */
                    arclengths.push_back(-1.0);
                }
            }
            if(raw.point_is_valid(tpv[j].i_tie,tpv[j].j_tie))
            {
                gptie=raw.geo_coordinates(tpv[j].i_tie,tpv[j].j_tie);
                dr=gptie.r-tpv[j].radius;
                drds=dr/arc;
                for(i=0;i<raw.n1;++i)
                {
                    if(raw.point_is_valid(i,j))
                    {
                        gp=raw.geo_coordinates(i,j);
                        gp.r += arclengths[i]*drds;
                        cp1=raw.gtoc(gp);
                        result.x1[i][j]=cp1.x1;
                        result.x2[i][j]=cp1.x2;
                        result.x3[i][j]=cp1.x3;
                    }
                }
            }
            else
            {
                cerr << "LinearDepthCorrection procedure (WARNING):  "
                    <<"tie point at grid position ("<<i_tie<<","
                       <<j_tie<<") is in a masked region"<<endl
                       <<"x2=const curve on grid line with j="
                       <<j_tie<< " may contain an offset"<<endl;
            }
        }
        return(result);
    } catch(...){throw;};
}




