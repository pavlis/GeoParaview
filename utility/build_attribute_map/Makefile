BIN=build_attribute_map
include $(ANTELOPEMAKE)
OBJS=build_attribute_map.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
