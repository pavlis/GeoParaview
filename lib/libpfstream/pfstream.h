typedef struct Pf_ensemble {
	Tbl *ensemble_keys;
	Tbl *group_keys;
	int ngroups,nmembers;
	int *group_start,*group_end;
	Pf **pf;
} Pf_ensemble;
typedef struct Attribute_map {
	char *dbname;
	int type;
	char *pfname;
} Attribute_map;
extern int DB2PFS_verbose;
enum Dbtype {DBINT, DBREAL, DBTIME, DBSTR};

#define ENDPF_SENTINEL "__EOF__"
#define END_OF_DATA_SENTINEL "__EOI__"

Pf_ensemble *create_Pf_ensemble(int, int);
void free_Pf_ensemble(Pf_ensemble *);
Pf *pfstream_read(char *);
Pf_ensemble *pfget_Pf_ensemble(Pf *,char *);
int pfwrite_stream(Pf *,char *,FILE *,int);
Pf *build_ensemble(Pf *,Pf_ensemble *,char *);
void pfensemble_put_double(Pf_ensemble *,char *,double);
void pfensemble_put_time(Pf_ensemble *,char *,double);
void pfensemble_put_int(Pf_ensemble *,char *,int);
void pfensemble_put_string(Pf_ensemble *,char *,char *);
Tbl *pfget_Attribute_map(Pf *,char *);

/* These are necessary because of deficiencies in pfput_double which
doesn't store enough digits with an epoch time */
double pfget_time(Pf *,char *);
void pfput_time(Pf *,char *,double);

/* this prototype declaration is necessary because pf.h does not 
define it for some reason */
Pf *pfdup(Pf *);

