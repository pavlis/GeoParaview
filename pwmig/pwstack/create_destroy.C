#include "stock.h"
#include "pf.h"
#include "pwstack.h"
Rectangular_Slowness_Grid::Rectangular_Slowness_Grid(Pf *pf)
{
	name = pfget_string(pf,(char *)"slowness_grid_name");
	uxlow = pfget_double(pf,(char *)"uxlow");
	uylow = pfget_double(pf,(char *)"uylow");
	nux = pfget_int(pf,(char *)"number_ux");
	nuy = pfget_int(pf,(char *)"number_uy");
}
Top_Mute::Top_Mute(Pf *pf)
{
	char *s;
	t0e = pfget_double(pf,(char *)"top_mute_endtime");
	t1 = pfget_double(pf,(char *)"top_mute_endtaper");
	s=pfget_string(pf,(char *)"top_mute_time_reference");
	if(!strcmp(s,(char *)"absolute"))
		reftype = absolute;
	else
		reftype = relative;
}
