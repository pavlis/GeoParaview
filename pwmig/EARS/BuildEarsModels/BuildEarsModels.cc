#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include "stock.h"
#include "coords.h"
//#include "dbpp.h"
using namespace std;
//using namespace SEISPP;
/* Assume awk is used used to convert Ears csv table to only the 
   following attributes:

   net  sta   lat   lon elev  crust_thickness crust_vp  crust_vs 

   Program reads from file with that structure and builds a 
   directory chain (defined by arg) with station files - names
   equal to station names.

   Complexity that made me do this in C++ is that we have to have a way to
   check for duplicate station names and reassign them names with net code
   when duplicates are found */

class EarsStaGeom
{
    public:
        EarsStaGeom(){};
        EarsStaGeom(char *n,char *s,double latin,double lonin,
                double elevin)
        {
            net=string(n);
            sta=string(s);
            lat=latin;  lon=lonin; elev=elevin;
        };
        EarsStaGeom(const EarsStaGeom& parent)
        {
            net=parent.net;
            sta=parent.sta;
            lat=parent.lat;
            lon=parent.lon;
            elev=parent.elev;
            net=parent.net;
        };
        string net;
        string sta;
        double lat,lon,elev;
        bool operator==(const EarsStaGeom& other);
};
bool EarsStaGeom::operator==(const EarsStaGeom& other)
{
    if(this->sta!=other.sta) return(false);
    /* cutoff distance in km */
    const double match_cutoff(0.5);  
    double delta,azimuth;
    /* I violate my standard here and do not convert
       lat and lon values to radians on input so here
       we have to it in this test. */
    dist(rad(lat),rad(lon),rad(other.lat),rad(other.lon),
            &delta,&azimuth);
    double deltakm=deg2km(deg(delta));
    if(deltakm>match_cutoff)
        return false;
    else
        return true;
}

string seedname(string net, string sta)
{
    string longname;
    if(sta.size()<4)
        longname=net+"_"+sta;
    else
        longname=net+sta;
    return(longname);
}
void usage()
{
    cerr <<"Usage error"<<endl
        << "BuildEarsModels infile directory "<<endl
        << " where infile is reformatted EARS csv file of crustal data" <<endl
        << " directory is the directory to write model files"<<endl;
    exit(-1);
}
int main(int argc, char **argv)
{
  try {
    if(argc!=3) usage();
    FILE *infile=fopen(argv[1],"r");
    if(infile==NULL)
    {
        cerr << "Open failed on file "<<argv[1]<<endl;
        usage();
    }
    /* Make sure we have a working directory to place files in */
    cout << "Creating directory "<<argv[2]<<endl;
    if(makedir(argv[2]))
    {
        cerr << "Cannot create that directory"<<endl;
        usage();
    }
    string dir(argv[2]);
    /*
    DatascopeHandle dbh(argv[3],true);
    AttributeMap am("css3.0");
    list<string> matchkeys;
    matchkeys.push_back(string("sta"));
    DatascopeMatchHandle dbm(dbh,string("site"),matchkeys,am);
    */
    /* We use this to index EarsStaGeom by sta - used to 
       check for duplicates. */
    map<string,EarsStaGeom> allstations;
    map<string,EarsStaGeom>::iterator it;
    char sta[20],net[5];
    double lat,lon,elev;
    double vpcrust,vscrust,hcrust;
    /* These need to be check against ears paper */
    const double vpmantle(8.2);
    const double vsmantle(4.5);
    const double rhocrust(2.7);
    const double rhomantle(3.3);
    while(fscanf(infile,"%s%s%lf%lf%lf%lf%lf%lf",
                net,sta,&lat,&lon,&elev,&hcrust,&vpcrust,&vscrust)
              != EOF)
    {
        string key(sta);
        string system_command;
        EarsStaGeom esgread(net,sta,lat,lon,elev);
        it=allstations.find(key);
        if(it==allstations.end())
            allstations[key]=esgread;
        else
        {
            cout << "Found duplicate for station name "<<key<<endl;
            cout << "Current line has net="<< net <<" sta="
                << sta<<" at coordinates:  "
                <<lat<<", "<<lon<<endl;
            cout << "Existing entry has net="<<(*it).second.net
                <<" sta="<<it->second.sta
                << " at coordinates:  "
                << it->second.lat <<", "<<it->second.lon<<endl;
            key=seedname(net,sta);
            cout <<"Saving as station name="<<key<<endl;
            esgread.sta=key;
            //allstations[key]=esgread;
            /* We need to also change the name of the one we 
               created first in this situation.  This may 
               fail but I think it will continue.*/
            system_command="mv "+dir+"/"+sta+" "
                +dir+"/"+seedname(it->second.net,it->second.sta);
            cout <<"STRIPTHISLINE: "<< system_command<<endl;

        }
        /* Now we write the velocity model using the station 
           name as they key (actually we use key) */
        string fname=dir+"/"+key;
        FILE *modfile=fopen(fname.c_str(),"w");
        if(modfile==NULL)
            cerr << "fopen failed on file = "<<fname<<endl
                <<"No model file written for "<<key<<endl;
        else
        {
            fprintf(modfile,"%lf %lf %lf %lf\n",
                    vpcrust,vscrust,rhocrust,hcrust);
            fprintf(modfile,"%lf %lf %lf 0.0\n",
                    vpmantle,vsmantle,rhomantle);
            fclose(modfile);
        }
    }
  } catch(std::exception& e)
  {
      cerr << "Fatal Error:  message returned follows."<<endl
          << e.what()<<endl;
  }
}


