BIN=dbpmel
DATADIR=	# directory under $(ANTELOPE)/data in which to put $(DATA) files
DATA=		# data files to be placed in data/$(DATADIR)
PF=dbpmel.pf

MAN1=	dbpmel.1
ldflags= 	# flags to linker/loader
#  This list is longer than necessary because the get_required_pf function is loaded
# from libmultiwavelet.  I need a glp utility library for this type function
ldlibs=-lgenloc -lmultiwavelet -xlic_lib=sunperf $(DBLIBS) $(F77LIBS) $(TRLIBS) -ltt -ldl -lm

include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 

OBJS=dbpmel.o dbform_working_view.o dbpmel_load_stations.o \
	dbsubs.o fixlist_utilities.o ftest.o shape.o matrix_subs.o \
	pmel.o pmel_process.o scmatrixsubs.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
