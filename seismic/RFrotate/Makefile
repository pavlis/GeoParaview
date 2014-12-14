BIN=RFrotate
MAN1=RFrotate.1
DATADIR=schemas/css3.0.ext
DATA=tmatrix
CXXFLAGS=-g
PF=RFrotate.pf

ldlibs=-lseispp -lgclgrid -ltrvltm $(DBLIBS) $(TRLIBS) -lperf -lm
include $(ANTELOPEMAKE)
OBJS=RFrotate.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
