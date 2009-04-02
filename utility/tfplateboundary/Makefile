BIN=tfplateboundary
PF=tfplateboundary.pf
#MAN1=tfplateboundary.1
cflags=-g
ldlibs=$(DBLIBS) -lseispp

include $(ANTELOPEMAKE)  	
OBJS=tfplateboundary.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
