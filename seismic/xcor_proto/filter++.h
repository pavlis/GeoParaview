namespace SEISPP
{
enum Filter_Type {highpass,lowpass,bandpass,WAA,WAV,WAD,DIF,DIF2,INT,INT2,DEMEAN};

class Time_Invariant_Filter
{
public:
	Filter_Type type;
	double fmin();
	double fmax();
	int fmin_poles();
	int fmax_poles();
	string type_description();
	Time_Invariant_Filter(){f1=0.0;f2=0.0;npole1=0;npole2=0;};
	Time_Invariant_Filter(string);
	Time_Invariant_Filter(double flow, int npl, double fhigh, int fph);
	Time_Invariant_Filter(const Time_Invariant_Filter&);
	Time_Invariant_Filter& operator=(const Time_Invariant_Filter&);
	void apply(int ns, double *s,double dt);
	void apply(int ns, float *s,double dt);
	void apply(Time_Series& ts);
	void apply(Three_Component_Seismogram& tce);
	void apply(Dbptr tr);

private:
	string filter_spec;
	// keep these here to avoid having to parse this 
	double f1,f2;  // low and high corner respectively
	int npole1,npole2;
};
// helpers
void Filter_Ensemble(Time_Series_Ensemble& ensemble,
		Time_Invariant_Filter& filter);
void Filter_Ensemble(Three_Component_Ensemble& ensemble,
		Time_Invariant_Filter& filter);
}  // End namespace SEISPP


