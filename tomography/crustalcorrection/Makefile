BIN=crustalcorrection
PF=crustalcorrection.pf
MAN1=crustalcorrection.1
ldlibs=$(DBLIBS) -lseispp

include $(ANTELOPEMAKE)  	
OBJS=crustalcorrection.o \
   StationVariableVelocityModel.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
