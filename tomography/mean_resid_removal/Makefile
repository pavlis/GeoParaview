BIN=mean_resid_removal
#MAN1=mean_resid_removal.o
ldlibs=$(DBLIBS) -lseispp

include $(ANTELOPEMAKE)  	
OBJS=mean_resid_removal.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -g -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
