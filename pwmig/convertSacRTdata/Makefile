BIN=convertRTdata
INCLUDE=
PF=
#MAN1=		# man pages for commands

#cflags=-g -DSUNPERF
cflags=-g 
ldlibs= -lseispp -lgclgrid -ltrvltm $(DBLIBS) $(TRLIBS) -lperf -lm

include $(ANTELOPEMAKE)
DIRS=

OBJS=convertRTdata.o
$(BIN) : $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
