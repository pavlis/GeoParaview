BIN=testmsw

    CXXFLAGS += -I/Users/pavlis/foreignsrc/libmseed
    LDFLAGS += -L/Users/pavlis/foreignsrc/libmseed -lmseed
OBJS=testmsw.o MSWriter.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
