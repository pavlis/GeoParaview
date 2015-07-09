LIB=libseisppplot.a
DLIB=$(LIB:.a=$(DSUFFIX)) 
INCLUDE=SeismicPlot.h cwpcmaps.h xplot.h

cflags=-G
ldflags=
ldlibs=

include $(ANTELOPEMAKE)

OBJS=SeismicPlot.o SeismicPick.o colormap.o image.o window.o
$(LIB) : $(OBJS)
	$(RM) $@
	$(CXXAR) $(CXXARFLAGS) $@ $(OBJS)  
	$(RANLIB) $@

$(DLIB) : $(OBJS)
	$(RM) $@
	$(CXXDLD) $(CXXDLDFLAGS) -o $(DLIB) $(OBJS) $(DLDLIBS) 
