typedef struct spectrum_
{
	char sta[8],chan[10];
	char sname[6];
	int arid;
	double time;
	int chanid;
	int jdate;
	double endtime;
	char phase[10];
	double freq0;
	int nfreq;
	double df;
	double rayleigh;
	double tbp;
	double sscale;
	char auth[15];
	char pstype[10];
	char datatype[4];
	char dir[66];
	char dfile[40];
	int index;
	float *spec;
}  Spectrum;
Spectrum read_spectrum();
Spectrum find_spectrum();
Spectrum resample_spectrum();
