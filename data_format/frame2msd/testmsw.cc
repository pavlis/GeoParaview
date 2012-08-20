#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "MSWriter.h"
using namespace std;
int main(int argc, char **argv)
{
    string ifile("testi.dat");
    string dfile("testd.dat");
    string ffile("testf.dat");
    string sfile("tests.dat");
    
    vector<int> di;
    vector<short> ds;
    vector<float> df;
    vector<double> dd;

    int nsamp(10000);
    double period,dt,t0;
    dt=0.05;
    t0=1000000.0;  // arbitrary
    int i;
    for(i=0;i<nsamp;++i)
    {
        short int val;
        val=static_cast<short int>(5000.0*sin(2.0*3.14*dt
                    *static_cast<double>(i)/2.0));
        ds.push_back(val);
        di.push_back(static_cast<int>(val));
        df.push_back(static_cast<float>(val));
        dd.push_back(static_cast<double>(val));
    }
    /*
    for(i=0;i<nsamp;++i)
    {
        cout << i <<" " <<ds[i] << " " 
            << di[i] << " "
            << df[i] << " "
            << dd[i] << endl;
    }
    */
    // these are fixed.  sta changes for each writer
    string net("TT"),chan("BHZ"),loc("00");
    double samprate=1/dt;
    MSWriter *msw;

    try{
    cout << "Trying to write float data file"<<endl;
    msw=new MSWriter(ffile);
    msw->write<float>(&(df[0]),samprate,nsamp,t0,net,string("FSTA"),chan,loc);
    // write second record with fixed offset 
    msw->write<float>(&(df[0]),samprate,nsamp,t0+10000.0);
    delete msw;

    cout << "Trying to write double data file"<<endl;
    msw=new MSWriter(dfile);
    msw->write<double>(&(dd[0]),samprate,nsamp,t0,net,string("DSTA"),chan,loc);
    // write second record with fixed offset 
    msw->write<double>(&(dd[0]),samprate,nsamp,t0+10000.0);
    delete msw;
    cout << "Trying to write int data file"<<endl;
    msw=new MSWriter(ifile);
    msw->write<int>(&(di[0]),samprate,nsamp,t0,net,string("ISTA"),chan,loc);
    msw->write<int>(&(di[0]),samprate,nsamp,t0+10000.0);
    delete msw;

    cout << "Trying to write short int data file"<<endl;
    msw=new MSWriter(sfile);
    msw->write<short>(&(ds[0]),samprate,nsamp,t0,net,string("SSTA"),chan,loc);
    // write second record with fixed offset 
    msw->write<short>(&(ds[0]),samprate,nsamp,t0+10000.0);
    delete msw;

    cout << "Trying to reopen int data file and add data for another station"
        <<endl;
    msw=new MSWriter(ifile);
    msw->write<int>(&(di[0]),samprate,nsamp,t0+20000.0,
            net,string("IST2"),chan,loc);
    msw->write<int>(&(di[0]),samprate,nsamp,t0+40000.0);
    delete msw;

    } catch (exception& err)
    {
        cout << "Error:  message throw follows"<<endl
            << err.what()<<endl;
    }
}
