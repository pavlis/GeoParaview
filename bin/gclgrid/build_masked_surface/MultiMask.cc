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
    bool result(false);
    vector<GeoPolygonRegion>::iterator pptr;
    int i;
    for(pptr=polygon.begin(),i=0;pptr!=polygon.end();++pptr,++i)
    {
        bool testone=pptr->is_inside(lat0,lon0);
        if(inside[i])
            result = result & testone;
        else
            result = result & (!testone);
    }
    return result;
}

