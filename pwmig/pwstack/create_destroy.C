#include "pf.h"
#include "pwstack.h"
Rectrangular_Slowness_Grid::Rectangular_Slowness_Grid(Pf *pf)
{
	name = pfget_string(pf,"slowness_grid_name");
	uxlow = pfget_double(pf,"uxlow");
	uylow = pfget_double(pf,"uylow");
	nux = pfget_double(pf,"number_ux");
	nuy = pfget_double(pf,"number_uy");
}
Top_Mute::Top_Mute(Pf *pf)
{
	char *s;
	t0e = pfget_double(pf,"top_mute_endtime");
	t1 = pfget_double(pf,"top_mute_endtaper");
	s=pfget_string(pf,"top_mute_time_reference");
	if(!strcmp(s,(char *)"absolute")
		reftype = absolute;
	else
		reftype = relative;
}
