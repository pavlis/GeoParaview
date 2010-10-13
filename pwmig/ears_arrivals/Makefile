BIN=ears_arrivals

cflags=-g
ldflags=$(DBLIBS) -lseispp

include $(ANTELOPEMAKE)  	
OBJS=ears_arrivals.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
