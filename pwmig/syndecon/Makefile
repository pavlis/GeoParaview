BIN=syndecon

cflags=-g 
CPPFLAGS=-g
ldlibs= -L/opt/SUNWspro/lib -lseispp -ltrvltm -lpfstream -lbrttutil $(DBLIBS) -lperf -lgclgrid
include $(ANTELOPEMAKE)
DIRS=
OBJS=syndecon.o sac.o
$(BIN) : $(OBJS)
	$(CXX) $(CFLAGS)  -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
