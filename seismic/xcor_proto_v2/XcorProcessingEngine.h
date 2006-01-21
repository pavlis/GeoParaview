//  Peng Wang, HPC/RAC/UITS
//  Indiana University
//
//  Copyright 2005, The Trustees of Indiana University.
//  Last Modified: 12/1/2005

#ifndef __XCOR_ENGINE_H
#define __XCOR_ENGINE_H

#include <sstream>
#include <vector>

#include "stock.h"
#include "pf.h"
#include "Metadata.h"
#include "resample.h"
#include "ensemble.h"
#include "seismicarray.h"
#include "filter++.h"
#include "MultichannelCorrelator.h"
#include "SeismicPlot.h"
#include "AnalysisSetting.h"
// Maintenance problem here
// This keyword is used in array_get_data.cc. 
// There was no clear include file that it could be placed
// into.
//
const string predicted_time_key("predarr.time");


using namespace std;
using namespace SEISPP;


class XcorProcessingEngine {
public:
	XcorProcessingEngine(Pf * global_pf, AnalysisSetting asinitial,
		string waveform_db_name, string result_db_name);
	~XcorProcessingEngine();  // necessary unless we can get rid of mcc raw pointer
        void change_analysis_setting(AnalysisSetting a) {analysis_setting=a; if(!analysis_setting.rw_set) analysis_setting.robust_tw=analysis_setting.beam_tw;}

	//Display results using SeismicPlot
	// These will likely change what they do.  
	// This programming interface may change some, but should not be dramatic.
        void display_data();
	void display_beam();
	void display_correlations();
	PointPick pick_beam();
	void do_picks();
	void edit_data();
	void shift_arrivals(double tshift);
	// End interface routines likely to change
	
	//save the resulting beam to the output database
        void save_results(int evid, int orid);

	void load_data(Hypocenter& hypo);
	// We should use a shared_ptr here so we keep an internal
	// reference as well as one visible externally.  This
	// avoids copying and uses resource management class to
	// simplify memory management.  I elected to not do this
	// for now as Sun's stock compiler does not seem to have std::tr1
	// For now beware this returns a pointer to the result stored in 
	// the private area.
	MultichannelCorrelator *analyze();

	// this method applies another filter to current data
	// Use change_analysis setting to alter base filter
	void filter_data(TimeInvariantFilter f);
	
	//validate analysis setting
	bool validate_setting(AnalysisSetting & a) { return (a.aw_set && a.bw_set && a.rt_set); } 

	TimeSeries get_beam();
	TimeSeriesEnsemble get_waveforms();
	TimeSeriesEnsemble get_raw_data();
	MultichannelCorrelator *get_correlation_results();

private:

	//This is the meta data object that stores the global static parameters.
	Metadata global_md;
	DatascopeHandle waveform_db_handle;
	DatascopeHandle result_db_handle;
	// Save the table references for each output table used by this program.
	// This avoids constant lookups.  These are saved as raw datascope Dbptr
	// structures instead of a DatascopeHandle because they are hidden behind
	// the interface anyway.
	Dbptr dbassoc;
	Dbptr dbarrival;
	Dbptr dbxcorarrival;
	Dbptr dbxcorbeam;
	Dbptr dbevlink;
	Dbptr dbwfprocess;
	bool save_arrival;  // arrival/assoc are not saved unless this is true
	//regular_gather is raw data saved here to avoid needing to reconstruct it
	// repeatedly when filtering changes.  regular_gather has uniform sampling.
	// Maintained as as an auto_ptr for convenience and efficiency.
	auto_ptr<TimeSeriesEnsemble> regular_gather;
	// This holds the working data.  
	TimeSeriesEnsemble waveform_ensemble;

	AnalysisSetting analysis_setting;
	MultichannelCorrelator * mcc;
	
	//This is the meta data object to hold SeismicPlot display parameters.
	Metadata data_display_md;
	// For now I use raw pointers to hold plot objects.  This will
	// need to change eventually.  Allows me to get past current Xwindow
	// problem in Peng's fancier interface.
	SeismicPlot *dataplot;
	Metadata xcor_display_md;
	SeismicPlot *xcorplot;
	Metadata beam_display_md;
	SeismicPlot *beamplot;

	// Updates geometry when needed.  Multiple methods need this so 
	// and require access to the private area so we put it here.
        void UpdateGeometry(TimeWindow twin);
	SeismicArray stations;

        // These are easily wrapped into constructor.  Method exist to build them
        // from a pf.   Their purpose is to define what attributes are extracted from
        // or written to the database.
        MetadataList ensemble_mdl;
        MetadataList trace_mdl;
        MetadataList beam_mdl;  // Used for array trace output
        AttributeMap am;
	string beam_directory;
	string beam_dfile;

	ResamplingDefinitions rdef;
	double target_dt;
        bool trim_resampled;  // set true of edge transients automatically cut in resampling
	// Added to private data 1/2/2006 by GLP.  Part of analysis setting now.
	TimeWindow current_data_window;
	TimeWindow raw_data_twin,regular_gather_twin;
	string netname;
	double xcorpeak_cutoff,coherence_cutoff,stack_weight_cutoff;
	double time_lag_cutoff;
};

#endif
