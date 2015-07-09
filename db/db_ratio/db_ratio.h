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

typedef struct spec_arrays_
{
	int nfreq;
	int nstation[4];
	int n1gather;  /* This is the fixed Fortran style leading dimension of char arrays*/
	float *chan[4];
}  Spectrum_gather;

#define SP1C "sp1c"
#define SP3C "sp3c"
#define RATIO1 "rps1"
#define RATIO3 "rps3"
/* Need these function prototypes */
Spectrum Spectrum_copy(Spectrum);
Spectrum_gather data_gather(Spectrum *spec,
        int nspec,
        char *sdenom,
        char *chanids[6],
        Spectrum threeCspec[],
        int *nthreeCspec);
int init_plot(char *);
Arr *MTSload_station_table(Dbptr);
float median_float(float *, int);
Spectrum read_spectrum(Dbptr);
long write_spectrum(float *,int,char *,char *);
int save_spectrum(Dbptr,Spectrum *);
void save_spectrum_error(Spectrum *,int);
Spectrum resample_spectrum(Spectrum, Spectrum);
