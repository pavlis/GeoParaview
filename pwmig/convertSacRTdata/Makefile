BIN=convertRTdata
INCLUDE=
PF=
#MAN1=		# man pages for commands

#cflags=-g -DSUNPERF
# use this for Sun using old matlab version.  
#cflags=-g -DOLDMATLAB
cflags=-g
ldlibs= -lseispp -lgclgrid -ltrvltm $(DBLIBS) $(TRLIBS) -lperf -lm -lmatlabprocessor  $(MATLABLIB) -leng -lmx -lut -lmat

include $(ANTELOPEMAKE)
include antelopemake.local
CXXFLAGS += $(MATLABINCLUDE)
DIRS=

OBJS=convertRTdata.o
$(BIN) : $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
