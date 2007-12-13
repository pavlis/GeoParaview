BIN=test

ldlibs=	$(DBLIBS)

include $(ANTELOPEMAKE) 

OBJS=test.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
