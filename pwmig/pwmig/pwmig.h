class Ray_Transformation_Operator
{
public:
	int npoints;
	Ray_Transformation_Operator(int np);
	// constructor for constant azimuth 
	Ray_Transformation_Operator(GCLgrid&g, dmatrix& path, double azimuth);
	// constructor for depth dependent operator
	Ray_Transformation_Operator(GCLgrid&g, dmatrix& path);
	Ray_Transformation_Operator(const Ray_Transformation_Operator& pat);
	dmatrix *apply(dmatrix& in);
private:
	vector<dmatrix> U;
};
// function prototypes 
Slowness_vector get_stack_slowness(Three_Component_Ensemble& ensemble);
Hypocenter get_event_hypocenter(Three_Component_Ensemble& ensemble);
GCLgrid3d *Build_GCLraygrid(bool fixed_u,
	GCLgrid& parent,Slowness_vector u,Hypocenter h,
	Velocity_Model_1d& vmod,double zmax, double tmax, double dt);
