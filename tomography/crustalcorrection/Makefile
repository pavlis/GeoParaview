BIN=crustalcorrection
PF=crustalcorrection.o
MAN1=crustalcorrection.o
ldlibs=$(DBLIBS) -lseispp

include $(ANTELOPEMAKE)  	
OBJS=crustalcorrection.o \
   StationVariableVelocityModel.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
