BIN=makegclgrid
INCLUDE=
PF=makegclgrid.pf
#MAN1=		# man pages for commands

cflags=-g -DSUNPERF
ldlibs= $(DBLIBS) -lgclgrid -xlic_lib=sunperf

include $(ANTELOPEMAKE)
DIRS=

OBJS=makegclgrid.o dr3.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
