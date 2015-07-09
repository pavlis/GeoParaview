BIN=sac2su
MAN1=sac2su.1
PF=sac2su.pf
cflags=-g
ldflags=-lseispp -lgclgrid -ltrvltm $(TRLIBS) $(DBLIBS) -lperf -lm

include $(ANTELOPEMAKE)
OBJS=	sac2su.o SacFileHandle.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
