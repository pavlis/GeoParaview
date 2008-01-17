BIN=project1dmod

cxxflags=-g 
ldflags=-lseispp -ltrvltm -lgclgrid -lseispp $(TRLIBS) $(DBLIBS) -lperf -lm

include $(ANTELOPEMAKE)

OBJS=project1dmod.o
$(BIN) : $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
