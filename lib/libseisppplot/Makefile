LIB=libseisppplot.a
DLIB=$(LIB:.a=$(DSUFFIX)) 
# For now only install xplot.h There are other include file here but think they are only used here.
MAN3=


ldlibs=$(X11LIBS)
include $(ANTELOPEMAKE)

OBJS=SeismicPick.o \
	SeismicPlot.o \
	axesbox.o \
	colormap.o \
	drawcurve.o \
	garnish.o \
	image.o \
	legendbox.o \
	rfwtva.o \
	rfwtvaint.o \
	rubberbox.o \
	window.o \
	xContour.o
$(LIB) : $(OBJS)
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LORDER) $(OBJS) $(TSORT)
	$(RANLIB) $@
