BIN=dbpmel
DATADIR=schemas2
DATA=genloc1.1
PF=dbpmel.pf

#MAN1=	dbpmel.1
ldflags=-L/opt/SUNWspro/lib
ldlibs= -lgenloc -lglputil $(DBLIBS) $(F77LIBS) \
-ltt -ldl -lm -xlic_lib=sunperf,fui,fai,fsu -dalign
include $(ANTELOPEMAKE)  	# This line must appear following the definitions above, 

OBJS=dbpmel.o dbpmel_load_stations.o \
	dbsubs.o fixlist_utilities.o ftest.o shape.o \
	pmel.o pmel_process.o scmatrixsubs.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
