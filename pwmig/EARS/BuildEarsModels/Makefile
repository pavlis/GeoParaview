BIN=BuildEarsModels

#MAN1=

cflags= -g
ldlibs= $(DBLIBS)

include $(ANTELOPEMAKE)  	

OBJS=BuildEarsModels.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
