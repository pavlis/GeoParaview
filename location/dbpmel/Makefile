BIN=dbpmel
PF=dbpmel.pf
MAN1=		# man pages for commands
cflags=		# -g or -O is the most common option here and in fflags

include $(ANTELOPEMAKE)  	
OBJS=	dbpmel.o \
	dbpmel_load_stations.o \
	dbsubs.o \
	matrix_subs.o \
	pmel.o \
	pmel_process.o \
	scmatrixsubs.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
