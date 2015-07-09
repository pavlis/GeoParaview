BIN=compare_solutions
cflags=	
ldflags= $(DBLIBS)

include $(ANTELOPEMAKE)  	
OBJS=compare_solutions.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
