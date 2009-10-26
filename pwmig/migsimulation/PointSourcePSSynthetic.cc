/* This object is the core of the implementation of this synthetic generator */
SphericalRayPathArray::SphericalRayPathArray(VelocityModel1d& vmod, vector<double> rayp,
        double zmax,  double dzin)
{
    const string base_error("SphericalRayPathArray constructor:  ");
    vector<double>::iterator rayptr;
    int nrays=rayp.size();
    if(nrays<=1) throw SeisppError(base_error
                +"input ray parameter vector is empty");
    /* rayp vector must be an increasing sequence from 0 */
    if(rayp[0]!=0.0) throw SeisppError(base_error
                + "First point in ray parameter vector must be 0");
    for(rayptr=rayp.begin()+1;rayptr!=rayp.end();++rayptr)
        if((*rayptr)-(*(rayptr-1))<0.0) throw SeisppError(base_error
                + "ray parameter vector must be an increasing sequence");
    dz=dzin;
    rays.reserve(nrays);
    try {
        const string mode("z");
        const vfortmax(2.0);  // needs to be smaller than anything realistic 
        tmax=zmax/vfortmax;
        slow=rayp;
        for(rayptr=rayp.begin();rayptr!=rayp.end();++rayptr)
        {
            rays.push_back(vmod,*rayptr,zmax,tmax,dz,mode);
        }
        /* Next need to calculate ray theory amplitudes -- need to consult references.
        first step is to determine nz as the minimum depth found for all rays */
        vector<RayPathSphere>::iterator raypaths;
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
    } catch(...){throw;};
}

pair<double,double> SphericalRayPathArray::time_and_amp(double delta, double z);
{
    const string base_error("SphericalRayPathArray::time_and_amp method:  ");
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

