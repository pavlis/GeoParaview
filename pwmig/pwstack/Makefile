BIN=makegclgrid
PF=makegclgrid.pf

MAN1=makegclgrid.1
cflags=-g
ldflags=$(DBLIBS)

include $(ANTELOPEMAKE) 
OBJS=pwstack.o \
	apsubs.o \
	dbsubs.o \
	pwstack_process.o \
	trace_subs.o 
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
