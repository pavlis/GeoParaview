namespace SEISPP
{
enum Filter_Type {highpass,lowpass,bandpass,WAA,WAV,WAD,DIF,DIF2,INT,INT2,DEMEAN};

class Time_Invariant_Filter
public:
	Filter_Type type;
	double fmin();
	double fmax();
	int fmin_poles();
	int fmax_poles();
	string type_description();
	Time_Invariant_Filter(string);
	void apply(int ns, double *s);
	void apply(int ns, float *s);
	void apply(Time_Series& ts);
	void apply(Three_Component_Seismogram& tce);
	void apply(Dbptr tr);

private:
	string filter_spec;
	// keep these here to avoid having to parse this 
	double f1,f2;  // low and high corner respectively
	int npole1,npole2;
};

}  // End namespace SEISPP


