BIN=xcor_proto
PF=xcor_proto.pf
#MAN1=
#ldlibs=-lseispp -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid
ldlibs=-lseispp -ltrvltm $(TRLIBS) $(DBLIBS) -lperf -lgclgrid -ltrvltm

CXX=g++ -g
include $(ANTELOPEMAKE)  	
OBJS=xcor_proto.o \
	filter.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
