/* This object is the core of the implementation of this synthetic generator */
SphericalRayPathArray::SphericalRayPathArray(VelocityModel1d& vmod, vector<double> rayp,
        double zminin, double zmax,  double dzin)
{
    const string base_error("SphericalRayPathArray constructor:  ");
    vector<double>::iterator rayptr;
    int nrays=rayp.size();
    if(nrays<=1) throw SeisppError(base_error
                +"input ray parameter vector is empty");
    for(rayptr=rayp.begin()+1;rayptr!=rayp.end();++rayptr)
        if((*rayptr)-(*(rayptr-1))<0.0) throw SeisppError(base_error
                + "ray parameter vector must be an increasing sequence");
    /* rayp vector must be an increasing sequence from 0 */
    if(fabs(ray[0])>FLT_EPSILON) throw SeisppError(base_error
                + "First point in ray parameter vector must be 0");
    zmin=zminin;
    dz=dzin;
    rays.reserve(nrays);
    try {
        const string mode("z");
        const vfortmax(2.0);  // needs to be smaller than anything realistic 
        tmax=zmax/vfortmax;
        slow=rayp;
        for(rayptr=rayp.begin();rayptr!=rayp.end();++rayptr)
        {
            rays.push_back(RayPathSphere(vmod,*rayptr,zmax,tmax,dz,mode));
        }
        /* Next need to calculate ray theory amplitudes -- need to consult references.
        first step is to determine nz as the minimum depth found for all rays */
        vector<RayPathSphere>::iterator raypaths,nextpath;
        raypaths=rays.begin();
        nz=raypaths->npts;
        ++raypaths;
        for(;raypaths!=rays.end();++raypaths)
        {
            if( (raypaths->npts)<nz) nz=raypaths->npts;
        }
        zfloor=dz*static_cast<double>(nz-1);
        if(SEISPP_verbose) cout << "SphericalRayPathArray:  using grid of "
                            << nz << " points spaced at "<<dz
                                <<" for total depth of "<<zfloor<<endl;
        /* Now set the amplitudes using ray geometric spreading */
        amps=dmatrix(nz,nrays);
        amps.zero();
        int i,j;
        double z,tani0,dx,dpdx;
        const double nullamp(-99999.);
        /* nrays -1 is the loop range because we use a forward 
        difference formula for amplitude calculation.  Last two
        columns will be identical.  Not ideal, but adequate
        for planned use of this */
        for(j=0,raypaths.begin();j<(nrays-1);++j,++raypaths)
        {
            for(i=0;i<nz;++i)
            {
                z=raypaths->depth(i);
                /* flag amp values above zmin with a negative value */
                if(z<zmin)
                    amp(i,j)=nullamp;
                else
                {
                    nextpath=raypath;
                    ++nextpath;
                    dpdx=(nextpath->p) - (raypath->p);
                    // here dx is delta increment
                    dx=(nextpath->delta[i])-(raypaths->delta[i]);
                    dx*=raypaths->r[0];
                    if(dx<=0.0) throw SeisppError(base_error
                            + "Ray grid error.  incremental delta .le. 0.0");
                    dpdx/=dx;
                    // now dx is used for local tangent computation 
                    if(i==(nz-1))
                        dx=raypaths->delta[nz-1]-raypaths->delta[nz-2];
                    else
                        dx=raypaths->delta[i+1]-raypaths->delta[i];
                    dx*=raypaths->r[0];
                    if(j==0)
                        amp(i,j)=dpdx;
                    else
                    {
                        tani0=dz/dx;
                        amp(i,j)=tani0*dpdx/(raypaths->delta[i]*raypaths->r[0]);
                    }
                }
            }
        }
        /* copy last column amplitude */
        for(i=0;i<nz;++i) amps(i,nrays-1)=amps(i,nrays-2);
        /* Amplitude right now is energy.  Normalize to p=0 value at first
           valid value and convert to amplitude by square root factor */
        double ampnormalizer;
        for(int i=0;i<nz;++i)
            if(amps(i,0)>0.0)
            {
                ampnormalizer=amps(i,0);
                break;
            }
        for(j=0;j<nrays;++j)
            for(i=0;i<nz;++i)
                if(amps(i,j)>0.0)
                    amps(i,j)=sqrt(amps(i,j)/ampnormalizer);

    } catch(...){throw;};
}
SphericalRayPathArray::SphericalRayPathArray(VelocityModel1d& vmod, 
        double p0, double dp, int np, 
        double zminin, double zmax,  double dzin)
{
    vector<double> pvector;
    pvector.reserve(np);
    double p;
    int i;
    for(i=0,p=0.0;i<np;++i,p+=dp)pvector.push_back(p);
    *this=SphericalRayPathArray(vmod,pvector,zminin,zmax,dzin);
}

pair<double,double> SphericalRayPathArray::time_and_amp(double delta, double z);
{
    const string base_error("SphericalRayPathArray::time_and_amp method:  ");
    if(z<zmin) throw SeisppError(base_error
            + "Requested point source depth is smaller than zmin.\n"
          + "Edit pf file.");
    if(z>zfloor) throw SeisppError(base_error
            + "Requested depth is below depth floor for this ray grid ");
    double time,amp;
    // scan the constant z position immediately above z to bracket the delta requested
    vector<RayPathSphere>::iterator raypaths,rayupper,raylower;
    int iupper,ilower;
    int iz=static_cast<int>(z/dz);
    for(iupper=0,raypaths=rays.begin();raypaths!=rays.end();++raypaths,++iupper)
        if(raypaths->delta[iz]>delta) break;
    rayupper=raypaths-1;
    for(ilower=0,raypaths=rays.begin();raypaths!=rays.end();++raypaths,++ilower)
        if(raypaths->delta[iz+1]>delta) break;
    raylower=raypaths-1;
    double wt[2][2];
    double a,dx,adx;
    a=1.0/(((rayupper+1)->delta[iz]-rayupper->delta[iz]));
    dx=delta-(rayupper->delta[iz]);
    adx=a*dx;
    wt[0][0]=1.0-adx;
    wt[0][1]=adx;
    a=1.0/(((raylower+1)->delta[iz+1]-raylower->delta[iz+1]));
    dx=delta-(raylower->delta[iz+1]);
    adx=a*dx;
    wt[1][0]=1.0-adx;
    wt[1][1]=adx;
    a=1.0/dz;
    dx=z-static_cast<double>(iz-1)*dz;
    adx=a*dx;
    wt[0][0]*=(1.0-adx);
    wt[0][1]*=(1.0-adx);
    wt[1][0]*=adx;
    wt[1][1]*=adx;
    double sumwt(0.0);
    for(int=0;i<2;++i)
       for(int j=0;j<2;++j) sumwt+=wt[i][j];
    time=0.0;
    time+=wt[0][0]*rayupper->t[iz];
    time+=wt[0][1]*(rayupper+1)->t[iz];
    time+=wt[1][0]*raylower->t[iz+1];
    time+=wt[1][1]*(raylower+1)->t[iz+1];
    amp=0.0;
    amp+=wt[0][0]*amp(iupper,iz);
    amp+=wt[0][1]*amp(iupper,iz);
    amp+=wt[1][0]*amp(ilower,iz+1);
    amp+=wt[1][1]*amp(ilower,iz+1);
    return(pair<double,double>(time,amp));
}
SphericalRayPathArray::SphericalRayPathArray(const SphericalRayPathArray& parent)
{
    nrays=parent.nrays;
    rays=parent.rays;
    slow=parent.slow;
    zmin=parent.zmin;
    nz=parent.nz;
    dz=parent.dz;
    amps=parent.amps;
    zfloor=parent.zfloor;
}
SphericalRayPathArray& SphericalRayPathArray::operator=(
        const SphericalRayPathArray& parent)
{
    if(this!=parent)
    {
	    nrays=parent.nrays;
	    rays=parent.rays;
	    slow=parent.slow;
            zmin=parent.zmin;
	    nz=parent.nz;
	    dz=parent.dz;
	    amps=parent.amps;
	    zfloor=parent.zfloor;
    }
    return(*this);
}

PointSourcePSSynthetic::PointSourcePSSynthetic(VelocityModel_1d& vmod,
        Pf *pf)
{
    const base_error("PointSourcePSSynthetic constructor:  ");
    try{
        Tbl *t;
        t=pfget_tbl(pf,"ray_parameter_grid");
        if(t==NULL) throw SeisppError(base_error
                + "Missing required parameter ray_parameter_grid");
        int i;
        vector<double> raypvector;
        for(i=0;i<maxtbl(t);++i)
        {
            double rayp;
            char *line=(char *)gettbl(t,i);
            istringstream ss(line);
            ss >> rayp;
            raypvector.push_back(rayp);
        }
        freetbl(t,0);
        Metadata control(pf);
        double zmin=control.get_double("zmin");
        double zmax=control.get_double("zmax");
        double dz=control.get_double("dz");
        rays=SphericalRayPathArray(vmod,raypvector,zmax,dz);
        distance_cutoff=control.get_double("distance_cutoff");
        t=pfget_tbl(pf,"point_sources");
        if(t==NULL) throw SeisppError(base_error
                + "Missing required parameter point_sources");
        points.reserve(maxtbl(t));
        amp.reserve(maxtbl(t));
        for(i=0;i<maxtbl(t);++i)
        {
            char *line=(char *)gettbl(t,i);
            double lat,lon,depth,r,ampin;
            Geographic_point gp;
            istringstream ss(line);
            ss >> lat;
            ss >> lon;
            ss >> depth;
            ss >> ampin;
            gp.lat=rad(lat);
            gp.lon=rad(lon);
            r0=r0_ellipse(lat);
            gp.r=r0-depth;
            points.push_back(gp);
            amp.push_back(ampin);
        }
    } catch(...){throw;};
}

TimeSeries ComputeScalar(Hypocenter& hypo,
         double rlat, double rlon, double rz,string type)
{
    throw SeisppError(string("ComputeScalar method not implemented"));
}
TimeSeries ComputeScalar(const TimeSeries& parent,
         Hypocenter& hypo,double rlat, double rlon, double rz,string type)
{
    throw SeisppError(string("ComputeScalar method not implemented"));
}
ThreeComponentSeismogram Compute3C(Hypocenter& hypo,
         double rlat, double rlon, double rz,string units)
{
    throw SeisppError(string("Compute3C plain method not implemented"));
}
ThreeComponentSeismogram Compute3C(const ThreeComponentSeismogram& parent,
         Hypocenter& hypo,double rlat, double rlon, double rz,string units)
{
    try{
    } catch(...){throw;};
}

