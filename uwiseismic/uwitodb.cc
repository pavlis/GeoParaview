/* in redhat, at least, this defines int16_t16_t and int16_t32_t used below */
#include <stdint.h>
#include <stdio.h>
#include <ios>
#include <iostream>
#include <string>
#include <list>
#include <map>
#include "stock.h"
#include "elog.h"
#include "db.h"
#include "pf.h"
using namespace std;
string  day_to_date(int32_t day)
{
	int32_t	d;
	int32_t	y;
	int32_t	m;

	y = (4 * day - 1) / 146097;
	day = (4 * day - 1) % 146097;
	d = day / 4;
	day = (4 * d + 3) / 1461;
	d = (4 * d + 3) % 1461;
	d = (4 + d) / 4;
	m = (5 * d - 3) / 153;
	d = (5 * d - 3) % 153;
	d = (5 + d) / 5;
	y = 100 * y + day;
	if (m < 10)
		m += 3;
	else	{
		m -= 9;
		y++;
	}
	char out[128];
	//original
	//strfmt(out, "%04d/%02d/%02d", (int16_t) y, (int16_t) m, (int16_t) d);
	//
	sprintf(out, "%02d/%02d/%04d", m,d,y);
	string result(out);
	return(result);
}
string time_to_date(int32_t csecs)
{
	int	h = (csecs / 360000);
	int	m;
	int	s;
	int	ss;

	csecs = csecs % 360000;
	m = (csecs / 6000);
	csecs = csecs % 6000;
	s =  (csecs / 100);
	ss = (csecs % 100);
	char out[128];
	/*
	//strfmt(out, "%02d:%02d:%02d.%02d", h, m, s, ss);
	*/
	sprintf(out,"%02d:%02d:%02d.%02d", h, m, s, ss);
	string result(out);
	return(result);
}
int32_t	julian_day(int16_t y,int16_t m,int16_t d)
{
	int16_t	c;
	int16_t	ya;
	int32_t	j;

	if (m > 2)
		m -= 3;
	else	{
		m += 9;
		y--;
	}
	c = y / 100;
	ya = y % 100;
	j = (146097 * c) / 4 + (1461 * ya) / 4 + (153 * m + 2) / 5 + d;
	return(j);
}

int32_t day_number(int8_t *time)
{
	int16_t	y = 100 * static_cast<int16_t>(*time) + static_cast<int16_t>(*(time + 1));
	int16_t	m = static_cast<int16_t>(*(time + 2));
	int16_t	d = static_cast<int16_t>(*(time + 3));

   return(julian_day(y, m, d));
}
int32_t	centiseconds(int8_t *time)
{
	int16_t	h = static_cast<int16_t>(*(time + 4));
	int16_t	m = static_cast<int16_t>(*(time + 5));
	int16_t	s = static_cast<int16_t>(*(time + 6));
	int16_t	ss = static_cast<int16_t>(*(time + 7));

	return(100 * (60 * (60 * h + m) + s) + ss);
}
void usage(char *prog)
{
	cerr << prog << " [-pf pffile] file1 file2 ... filen db"<<endl;
	exit(-1);
}
typedef struct UWIHEADER {
	int32_t next;
	int32_t prev;
	int16_t junk1;
	int16_t junk2;
	int16_t segment_no;
	int8_t dtype;
	char name[5];
	int16_t count;
	int16_t flag;
	int32_t time;
	int32_t window_no;
	int8_t debris[16];
} UWIHEADER;
void print_header(UWIHEADER& head)
{
  cout << "segment_no="<<head.segment_no<<", "
     << "name="<<head.name<<", "
     << "nsamp="<<(head.count)*256<<", "
     << "time="<<day_to_date(head.window_no)<<" "
	<< time_to_date(head.time)<<endl;
}
/* csec can be negative or larger than 24 hours so this
system requires this code to repair it.  Adapted from newpb
source code. Note window_no here holds the day field*/
void repair_times(UWIHEADER& head) 
{
	while(head.time<0)
	{
		head.window_no--;
		head.time += 8640000;
	}
	while(head.time>= 8640000)
	{
		head.window_no++;
		head.time -= 8640000;
	}
}
int main(int argc, char **argv)
{
	ios::sync_with_stdio();
	if(argc<3) usage(argv[0]);
	string dbname(argv[argc-1]);
	cout << "Writing data to database = "<<dbname<<endl;
	Dbptr db;
	if(dbopen(const_cast<char *>(dbname.c_str()),"r+",&db)==dbINVALID) 
	{
		cerr << "dbopen failed on "<<dbname<<endl;
		usage(argv[0]);
	}
	db=dblookup(db,0,"wfdisc",0,0);
	list<string> fnlist;
	string pffile("uwitodb");
	int i;
	for(i=1;i<argc-1;++i)
	{
		string arg(argv[i]);
		if(arg=="-pf")
		{
			++i;
			pffile=string(argv[i]);
		}
		else
		{
			fnlist.push_back(string(argv[i]));
		}
	}
	if(fnlist.empty()) usage(argv[0]);
	// parse parameter file
	Pf *pf;
	if(pfread(const_cast<char *>(pffile.c_str()),&pf) <0)
	{
		cerr << "pfread error on pf file "<< pffile<<endl;
		exit(-1);
	}
	// We build a map keyed by uwistation name for stations and channels
	map<string,string> stations;
	map<string,string> channels;
	Tbl *t;
	char *line;
	t=pfget_tbl(pf,"station_channel_map");
	if(t==NULL)
	{
		cerr << "required station_channel_map Tbl& is not in parameter file"<<endl;
		exit(-1);
	}
	for(i=0;i<maxtbl(t);++i)
	{
		char sta[10],chan[10],uwista[6];
		line=static_cast<char *>(gettbl(t,i));
		sscanf(line,"%s%s%s",uwista,sta,chan);
		stations[uwista]=string(sta);
		channels[uwista]=chan;
	}
		
	// constants for this format
	double samprate=100.0;
	char *datatype="i2";
	list<string>::iterator fname;
	FILE *fp;
	for(fname=fnlist.begin();fname!=fnlist.end();++fname)
	{
		int nread,nskip;
		UWIHEADER head;
		string sta;
		string chan;
		int nsamp;
		double time;
		double etime;
		long int foff;
		char dir[128],dfile[64];
		// Used by find functions applied to map
		map<string,string>::iterator ptrsta,ptrchan;
		cout << "Working on file "<<*fname<<endl;
		fp=fopen((*fname).c_str(),"r");
		if(fp==NULL)
		{
			cerr << "Open failed with file "<< *fname << " skipped"<<endl;
		}
		else
		{
			while((nread=fread(static_cast<void*>(&head),1,48,fp))>0)
			{
				repair_times(head);
				print_header(head);
				string tstring(day_to_date(head.window_no)+" "+time_to_date(head.time));
				time=str2epoch(const_cast<char *>(tstring.c_str()));
				etime=time+static_cast<double>(head.count*256)/samprate;
				foff=ftell(fp);
				// skip now as logic below can jump out
				if(fseek(fp,static_cast<long int>(head.count)*512,SEEK_CUR) )
				{
					cerr << "Seek error on file "<< *fname << endl
						<< "Fatal error.  File may be corrupted"<<endl;
					exit(-1);
				}
				string teststa(head.name);
				ptrsta=stations.find(teststa);
				ptrchan=channels.find(teststa);
				if( (ptrsta==stations.end()) || (ptrchan==channels.end()) )
				{
					sta=string(head.name);
					chan="EHZ";
				}
				else
				{
					sta=ptrsta->second;
					chan=ptrchan->second;
				}
				nsamp=head.count*256;
				dirbase(const_cast<char *>((*fname).c_str()),dir,dfile);
				cout << head.name << " converted to "
					<< sta << ":"<<chan<<endl;
				//
				// calib currently not computed.  Need clarification
				//
/*
				if(dbaddv(db,0,"sta",sta.c_str(),
						"chan",chan.c_str(),
						"time",time,
						"endtime",etime,
						"nsamp",nsamp,
						"samprate",samprate,
						"datatype",datatype,
						"dir",dir,
						"dfile",dfile,
						"foff",foff,0) < 0)
				{
					elog_complain(0,"problem saving %s:%s at time %s",
						sta.c_str(),chan.c_str(),
						tstring.c_str());
					cout << "dbaddv for "  
						<< sta<<":"<<chan<<endl
						<< "Will probably be missing from database" << endl;
				}
*/
				db.record=dbaddnull(db);
				dbputv(db,0,"sta",sta.c_str(),
                                                "chan",chan.c_str(),
                                                "time",time,
                                                "endtime",etime,
                                                "nsamp",nsamp,
                                                "samprate",samprate,
                                                "datatype",datatype,
                                                "dir",dir,
                                                "dfile",dfile,
                                                "foff",foff,0);

			}
		}
	}
}
