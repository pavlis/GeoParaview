#include "stock.h"
#include "pf.h"
#include "metadata.h"

Slowness_vector get_stack_slowness(Three_Component_Ensemble& ensemble)
		throw(Metadata_error)
{
	Slowness_vector slow;

	// We extract these from the global metadata area for the ensemle.
	// conceptually correct, but requires more careful implementation and assumes
	// the data are there.
	try {
		slow.ux = ensemble.md.get_double("plane_wave_ux");
		slow.uy = ensemble.md.get_double("plane_wave_uy");
	} catch (Metadata_error mde)
	{
		throw(mde);
	}
	return(slow);
}

Hypocenter get_event_hypocenter(Three_Component_Ensemble& ensemble)
{
	Hypocenter h;

	try {
		h.lat = ensemble.md.get_double("origin.lat");
		h.lon = ensemble.md.get_double("origin.lon");
		h.z = ensemble.md.get_double("origin.depth");
		h.time = ensemble.md.get_double("origin.time");
	} catch (Metadata_error mde)
	{
		throw mde;
	}
	return(h);
}
