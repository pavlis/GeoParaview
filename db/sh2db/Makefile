BIN=sh2db
DATADIR=schemas/css3.0.ext
DATA=orient
PF=sh2db.pf
MAN1=sh2db.1

ldlibs=-lseispp -lgclgrid -ltrvltm $(TRLIBS) $(DBLIBS) -lperf -lm
include $(ANTELOPEMAKE)  	
OBJS=sh2db.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
