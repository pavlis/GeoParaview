BIN=makegclgrid
INCLUDE=
PF=makegclgrid.pf
#MAN1=		# man pages for commands

cflags=-g
ldlibs=-xlic_lib=sunperf $(DBLIBS) -lgclgrid

include $(ANTELOPEMAKE)
DIRS=

OBJS=makegclgrid.o dr3.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
