/* This is an experiment to think about defining a GCLgrid-based array of 
three-component seismograms to hold data for plane wave migration algorithm */

class pwmig_input_data : public GCLgrid
{
	public:
// needs a special constructor to build this from an MPI input stream
// constructor should read data and simultaneously build the index.
		pwmig_input_data(){index=NULL;};
		pwmig_input_data(Pf *pf,Dbptr db);
		pwmig_input_data(Pf *pf,GCLgrid& g);
		~pwmig_input_data();
		get_seismogram(int i, int j); // key function, returns copy of seismogram at position i,j
	private:
		Three_Component_Ensemble seismicdata;
		int **index;  // element k = index(i,j) is seismicdata.tcse[k];
}
/*  Destructor only neds to handle the index array.  the seismicdata object is
built seperately and the GCLgrid components are handled by the GCLgrid destructor.
*/
pwmig_input_data::~pmmig_input_data()
{
	if(index!=NULL)
	{
		for(int i=0;i<n2;++i) delete [] index[i];
		delete index;
	}
}
/* In this constructor the grid is assumed to be defined in the database, db, 
and a named found in Pf is used to match the right grid to this object. 
I'm building this up front even though I may not use it because I'm not
sure if I can make the one below it (which I intend to use) work right.
*/
pwmig_input_data::pwmig_input_data(Pf *pf, Dbptr db)
{
	int i;
	char *gridname;
	try{
		gridname = pfget_string(pf,"gridname");
		if(gridname == NULL)
		{
			elog_die(0,(char *)"pwmig_input_data constructor failure:  Cannot find gridname parameter in parameter space\n");
		}
		GCLgrid gtmp(db,gridname);
		*this = gtmp;
	} catch (int ierr)
	{
		elog_die(0,(char *)"pwmig_input_data constructor failure:  database driven constructor threw error %d\n",ierr);
	}
	// Now the stuff dependent upon seismic data 
	// first create space for the index array
	new *int *index[n1];
	for(i=0;i<n2;++i)
	{
		new int index[i][n2];
	}
