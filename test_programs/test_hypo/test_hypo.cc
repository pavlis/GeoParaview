#include "seispp.h"
using namespace SEISPP;

int main(int argc, char **argv)
{
	Hypocenter hypo;
	hypo.lat=0.0;
	hypo.lon=1.0;
	hypo.z=10.0;
	hypo.time=100000.0;
	int i;
	double lat,lon;

	try {
		cout << "Testing distance, azimuth, and ptime functions" <<endl;
		cout << "distance esaz seaz ptime" << endl;
		for(i=0;i<20;++i)
		{
			lon=1.0;
			lat=0.0+0.25*static_cast<double>(i);
			cout << hypo.distance(lat,lon) 
				<< " " << hypo.esaz(lat,lon)
				<< " " << hypo.seaz(lat,lon)
				<< " " << hypo.ptime(lat,lon,0.0)
				<< endl;
		}
		cout << "testing slowness vector functions" << endl;
		cout << "distance esaz seaz pslow pazimuth Sslow Sazimuth "<<endl;
		Slowness_vector Pslow,Sslow;
		for(i=0;i<20;++i)
		{
			lon=1.0;
			lat=0.0+0.25*static_cast<double>(i);
			Pslow = hypo.pslow(lat,lon,0.0);
			Sslow = hypo.phaseslow(lat,lon,0.0,string("S"));
			cout << hypo.distance(lat,lon) 
				<< " " << hypo.esaz(lat,lon)
				<< " " << hypo.seaz(lat,lon)
				<< " " << Pslow.mag()
				<< " " << Pslow.azimuth()
				<< " " << Sslow.mag()
				<< " " << Sslow.azimuth()
				<< endl;
		}
	}
	catch (seispp_error err)
	{
		err.log_error();
		elog_die(0,"seispp_error");
	}
}
				

