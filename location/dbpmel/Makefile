BIN=dbpmel
DATADIR=	# directory under $(ANTELOPE)/data in which to put $(DATA) files
DATA=		# data files to be placed in data/$(DATADIR)
PF=dbpmel.pf

MAN1=	dbpmel.1
ldflags= 	# flags to linker/loader
ldlibs=-lgenloc -lsunperf $(DBLIBS)

include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 

OBJS=dbpmel.o dbform_working_view.o dbpmel_load_stations.o \
	dbsubs.o ftest.o matrix_subs.o pmel.o pmel_process.o scmatrixsubs.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
