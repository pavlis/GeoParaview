typedef struct Attribute_map {
	char *dbname;
	int type;
	char *oname;
} Attribute_map;
extern int DB2PFS_verbose;
enum Dbtype {DBINT, DBREAL, DBSTR};

Tbl *pfget_Attribute_map(Pf *,char *);
