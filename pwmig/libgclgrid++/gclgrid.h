#include "db.h"

#ifdef	__cplusplus

class GCLgrid
{
	public:
		char name[10];  // name assigned to this coordinate grid
		double lat0, lon0, r0;  // geographical location of origin 
		double azimuth_y;  // Azimuth of positive y axis 
		double dx1_nom, dx2_nom;  // nominal grid spacing 
		int n1,n2;  // grid size in each component direction
		int i0, j0;  // origin location in grid 
		double xlow, xhigh, ylow, yhigh, zlow, zhigh;// bounding box 
		int cartesian_defined, geographic_defined;
		double **x1, **x2, **x3; //cartesian coordinates of nodes
		double **lat, **lon, **r;  //geographical coordinates of nodes

		GCLgrid(){
			n1=0;n2=0;x1=NULL;x2=NULL;x3=NULL;lat=NULL;lon=NULL;r=NULL;
		};
		GCLgrid(int n1size, int n2size);
		GCLgrid(int n1size, int n2size, char *n, double la0, double lo0,
			double radius0, double az, double dx1n, double dx2n, 
			int iorigin, int jorigin);
		GCLgrid(Dbptr db, char *nm);  // acquire from Antelope database 
		GCLgrid(const GCLgrid&);  //standard copy constructor
		GCLgrid& operator=(const GCLgrid& );
		void dbsave(GCLgrid&, Dbptr, char *);
		~GCLgrid();
};
//3d version is identical except it requires 3 indexes instead of 2 for
//coordinates.  We use inheritance to simply this description.
class GCLgrid3d: GCLgrid
{
	public:
		double dx3_nom;
		int n3;
		int k0;
		double ***x1, ***x2, ***x3;
		double ***lat, ***lon, ***r;

		GCLgrid3d(){
			n1=0;n2=0;n3=0;
			x1=NULL;x2=NULL;x3=NULL;lat=NULL;lon=NULL;r=NULL;
		};
		GCLgrid3d(int n1size, int n2size, int n3size);
		GCLgrid3d(int n1size, int n2size, int n3size, 
			char *n, double la0, double lo0,
			double radius0, double az, 
			double dx1n, double dx2n, double dx3n,
			int iorigin, int jorigin);
		GCLgrid3d(Dbptr db, char *nm); 
		GCLgrid3d(const GCLgrid3d&); 
		GCLgrid3d& operator=(const GCLgrid3d& );
		void dbsave(GCLgrid3d&, Dbptr, char *);
		~GCLgrid3d();
};	  		
#endif

#ifdef	__cplusplus
extern "C" {
#endif

	

//
// helper function prototypes 
//
double r0_ellipse(double);
double ***create_3dgrid_contiguous(int, int, int);
double **create_2dgrid_contiguous(int, int);
#ifdef  __cplusplus
}
#endif

