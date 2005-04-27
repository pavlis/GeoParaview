BIN=dbremapgrid
INCLUDE=
PF=
#MAN1=		# man pages for commands
CC=CC
LD=CC

#cflags=-g -DSUNPERF
cflags=-g 
ldlibs= -L/opt/SUNWspro/lib $(DBLIBS) -lgclgrid -lperf

include $(ANTELOPEMAKE)
DIRS=

OBJS=dbremapgrid.o
$(BIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
