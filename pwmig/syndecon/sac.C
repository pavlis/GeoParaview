#include <iostream>
#include <fstream>
#include "seispp.h"
using namespace std;
const double sac_magic_real=-12345.0;
const int sac_magic_int=-12345;
const string sac_magic_string="-12345";

Time_Series* Read_SAC_ascii(string fname)
{
	double dhdr;
	int ihdr;
	string shdr; 
	ifstream infile(fname.c_str(), ios::in);
	if(!infile)throw SAC_data_error("Cannot open SAC ascii file");
	

	Time_Series *s = new Time_Series();
	Time_Series& ssac = *s;

	// The parameters in the header are in a fixed order.
	// first are reals, then ints, then strings 
	// The constants are used to skip nulls.

	infile >> dhdr;
	if(dhdr!=sac_magic_real) 
	{
		ssac.dt = dhdr;
		ssac.md.put_metadata("DELTA",dhdr);
	}
	else
	{
		delete s;
		throw SAC_data_error("Missing Required parameter dt");
	}
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("DEPMIN",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("DEPMAX",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("SCALE",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("ODELTA",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) 
	{
		ssac.t0 = dhdr;
		ssac.md.put_metadata("B",dhdr);
		ssac.tref=absolute;
	}
	else
	{
		ssac.t0=0.0;
		ssac.tref=relative;
	}
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("E",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("O",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("A",dhdr);
	infile >> dhdr;
	// skipped as format declares this one "internal"
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T0",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T1",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T2",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T3",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T4",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T5",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T6",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T7",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T8",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("T9",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("F",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP0",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP1",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP2",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP3",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP4",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP5",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP6",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP7",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP8",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("RESP9",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("STLA",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("STLO",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("STEL",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("STDP",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("EVLA",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("EVLO",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("EVEL",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("EVDP",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("MAG",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER0",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER1",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER2",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER3",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER4",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER5",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER6",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER7",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER8",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("USER9",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("DIST",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("AZ",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("BAZ",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("GCARC",dhdr);
	// format has two more internals here
	infile >> dhdr;
	infile >> dhdr;
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("DEPMEN",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("CMPAZ",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("CMPINC",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("XMINIMUM",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("XMAXIMUM",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("YMINIMUM",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("YMAXIMUM",dhdr);
	infile >> dhdr;
	if(dhdr!=sac_magic_real) ssac.md.put_metadata("ADJTM",dhdr);
	// format has 6 straight unused entries
	infile >> dhdr;
	infile >> dhdr;
	infile >> dhdr;
	infile >> dhdr;
	infile >> dhdr;
	infile >> dhdr;
	// now we do the ints 
	// should crack this date and turn it into an absolute time
	// but until this thing is really used with data I'll not bother.
	// my first task is only with synthetics.
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NZYEAR",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NZJDAY",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NZHOUR",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NZMIN",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NZSEC",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NZMSEC",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NVHDR",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NORID",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NEVID",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) 
	{
		ssac.ns =ihdr;
		ssac.md.put_metadata("NPTS",ihdr);
	}
	else
	{
		delete s;
		throw SAC_data_error("Missing Required parameter ns");
	}
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NSPTS",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NWFID",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NXSIZE",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("NYSIZE",ihdr);
	infile >> ihdr;
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IFTYPE",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IDEP",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IZTYPE",ihdr);
	infile >> ihdr;
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IINST",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("ISTREG",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IEVREG",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IEVTYP",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IQUAL",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("ISYNTH",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IMAGTYP",ihdr);
	infile >> ihdr;
	if(ihdr!=sac_magic_int) ssac.md.put_metadata("IMAGSRC",ihdr);
	infile >> ihdr;
	infile >> ihdr;
	infile >> ihdr;
	infile >> ihdr;
	infile >> ihdr;
	infile >> ihdr;
	infile >> ihdr;
	infile >> ihdr;

	infile >> ihdr;
        if(ihdr!=sac_magic_int) ssac.md.put_metadata("LEVEN",ihdr);
	infile >> ihdr;
        if(ihdr!=sac_magic_int) ssac.md.put_metadata("LPSPOL",ihdr);
	infile >> ihdr;
        if(ihdr!=sac_magic_int) ssac.md.put_metadata("LOVROK",ihdr);
	infile >> ihdr;
        if(ihdr!=sac_magic_int) ssac.md.put_metadata("LCALDA",ihdr);
	infile >> ihdr;
	// now the strings 
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KSTNM",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KEVNM",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KHOLE",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KO",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KA",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT0",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT1",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT2",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT3",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT4",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT5",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT6",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT7",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT8",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KT9",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KF",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KUSER0",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KUSER1",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KUSER2",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KCMPNM",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KNETWK",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KDATRD",shdr);
	infile >> shdr;
	if(shdr!=sac_magic_string) ssac.md.put_metadata("KINST",shdr);
//cout << fname << endl;
//ssac.md.print_all_metadata();
	ssac.live=true;

	ssac.s = new double[ssac.ns];
	for(int i=0;i<ssac.ns;++i)
	{
		infile >> ssac.s[i];
		if(infile.eof())
			throw SAC_data_error("EOF encountered prematurely");
	}
	return(s);
}
