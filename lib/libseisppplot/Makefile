BIN=testseisplot
ldlibs=$(X11LIBS) -lseispp $(TRLIBS) -lm

include $(ANTELOPEMAKE)  	
OBJS=testseisplot.o SeismicPick.o \
        SeismicPlot.o \
        colormap.o \
        image.o \
        window.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
