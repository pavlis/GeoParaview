
/* projects a ray in RayPathSphere object, which stored in
spherical coordinate form as pairs of delta and depth, to the 
Cartesian reference frame for the input GCLgrid (grid) at the
grid position i,j,k.  k should normally be the surface, but ther
algorithm does not require this.  It just does this blindly.

Returns a dmatrix of Cartesian coordinates for the path with
the coordinatres for each node stored in the columns of the
output matrix.
Author:  Gary Pavlis
Written:  May 2003
*/
dmatrix GCLproject_ray(GCLgrid3d *grid,RayPathSphere *ray,
	int i, int j, int k)
{

	dmatrix coords();
	dmatrix U();
	dmatrix newpath();

	coords = GCLgrid_Ray_gtoc(grid,ray);
	U=compute_unit_spherical_transformation(grid,
			grid->lat[i][j][k],grid->lon[i][j][k]);
	newpath=U*coords;
	for(int i=0;i<newpath.nc;++i)
	{
		newpath(0,i)-=grid->x1[i][j][k];
		newpath(1,i)-=grid->x2[i][j][k];
		newpath(2,i)-=grid->x3[i][j][k];
	}
	return(newpath);
}
		

