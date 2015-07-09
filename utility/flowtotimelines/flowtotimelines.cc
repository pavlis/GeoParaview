#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "coords.h"
/* This program takes the output of the slabmodel program and 
   builds a complementary mesh of lines that are time lines.  
   It works by assuming all lines originate at time zero and 
   are defined at a regular sample interval in time.   

   The output can optionally be decimated.
   */
using namespace std;
void usage()
{
    cerr << "flowtotimelines [-d decfactor] < in > out "<<endl;
    exit(-1);
}
int main(int argc, char **argv)
{
    int i;
    int decfac(1);
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-d")
        {
            ++i;
            std::istringstream sstr(argv[i]);
            sstr>>decfac;
            cerr << "Time lines will be decimated by a factor of "
                << decfac<<endl;
        }
        else
            usage();
    }
    typedef vector<double> Dvec;
    vector<Dvec> lat,lon,depth;
    vector<int> flowlinesize;
    int segsize,nseg;
    int maxsize(0);
    char line[256];
    /* Initialize these variables at top of the main read loop */
    Dvec fllat,fllon,flz;
    i=0;  nseg=0;
    while(!cin.eof())
    {
        double latin,lonin,zin;
        cin.getline(line,256);
        if(line[0]=='#') continue;
        if(cin.eof()) break;
        if(line[0]=='>')
        {
            segsize=fllat.size();
            if(segsize==0)
                cerr << "Warning:  zero length line segment number "<<nseg
                    <<endl << "Ignored."<<endl;
            else
            {
                maxsize=max(maxsize,segsize);
                flowlinesize.push_back(segsize);
                lat.push_back(fllat);
                lon.push_back(fllon);
                depth.push_back(flz);
                fllat.clear();
                fllon.clear();
                flz.clear();
                ++nseg;
            }
        }
        else
        {
            sscanf(line,"%lf%lf%lf",&lonin,&latin,&zin);
            fllat.push_back(latin);
            fllon.push_back(lonin);
            flz.push_back(zin);
        }
    }
    int j;
    for(i=0;i<maxsize;i=i+decfac)
    {
        segsize=0;
        for(j=0;j<nseg;++j)
        {
            if(i<flowlinesize[j])
            {
                cout << lon[j][i]<<" "
                    <<lat[j][i]<<" "
                    <<depth[j][i]<<endl;
                ++segsize;
            }
            else if(segsize>1)
            {
                // This will silently skip segments with only one point
                cout << ">"<<endl;
                segsize=0;
            }
        }
        cout << ">"<<endl;
    }
}

