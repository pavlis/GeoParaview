class RayPathSphere
{
public:
	int npts;
	double *r,*delta,*t;
	RayPathSphere(){r=NULL,delta=NULL,t=NULL;};
	RayPathSphere(int n)
	{npts=n; r=new double[n], delta=new double[n], t=new double[n];};
	// This fully parametrized version constructs a full path
	RayPathSphere(Velocity_Model&,double p, double zmax, double tmax,
			        double del, string mode);
	~RayPathSphere(){if(r!=NULL)delete[]r;
		if(delta!=NULL)delete[]delta;
		if(t!=NULL)delete[]t;};
};
