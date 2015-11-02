#include "MultiMask.h"
MultiMask::MultiMask(istream& in)
{
    try {
        char linebuf[512];
        string fname,polarity;
        while(in.getline(linebuf,512))
        {
            stringstream ss(linebuf);
            ss >> fname;
            ss >> polarity;
            GeoPolygonRegion poly(fname);
            polygon.push_back(poly);
            if(polarity=="outside")
                inside.push_back(false);
            else
                inside.push_back(true);
        }
    }catch(...){throw;};
}
MultiMask::MultiMask(const MultiMask& parent) : polygon(parent.polygon), inside(parent.inside)
{
}
MultiMask& MultiMask::operator=(const MultiMask& parent)
{
    if(this!=&parent)
    {
        polygon=parent.polygon;
        inside=parent.inside;
    }
    return *this;
}

bool MultiMask::is_inside(double lat0, double lon0)
{
    vector<GeoPolygonRegion>::iterator pptr;
    int i;
    for(pptr=polygon.begin(),i=0;pptr!=polygon.end();++pptr,++i)
    {
        bool testone=pptr->is_inside(lat0,lon0);
        cout << "Polygon "<<i<<" point "<<deg(lat0)<<" "<<deg(lon0)
            << " is_inside returned "<<testone<<endl;
        if(inside[i])
        {
            if(testone)
                continue;
            else
                return false;
        }
        else
        {
            /* block for polygons tagged to keep points outside */
            if(testone)
                return false;
            else
                continue;
        }
    }
    return true;
}

