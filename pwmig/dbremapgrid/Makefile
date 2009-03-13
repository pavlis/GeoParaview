BIN=dbremapgrid
INCLUDE=
PF=
#MAN1=		# man pages for commands
CC=CC
LD=CC

#cflags=-g -DSUNPERF
cflags=-g 
ldlibs= -lseispp -lgclgrid -ltrvltm $(DBLIBS) -lperf -lm -lgclgrid -lseispp

include $(ANTELOPEMAKE)
DIRS=

OBJS=dbremapgrid.o
$(BIN) : $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
