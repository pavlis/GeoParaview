LIB=libgeocoords.a
INCLUDE=GeoCoordError.h \
  GeoPath.h \
  GeoSplineSurface.h \
  GeoSurface.h \
  GeoTriMeshSurface.h \
  LatLong-UTMconversion.h \
  PLGeoPath.h \
  PlateBoundaryPath.h \
  RegionalCoordinates.h
cflags=-g`
ldlibs=-lgclgrid -lcoords -lcgeom $(DBLIBS)

include $(ANTELOPEMAKE)

OBJS=GeoSplineSurface.o \
  GeoTriMeshSurface.o \
  LatLong-UTMconversion.o \
  PLGeoPath.o \
  PlateBoundaryPath.o \
  RegionalCoordinates.o

$(LIB) : $(OBJS)
	$(RM) $@
	$(CXXAR) $(CXXARFLAGS) $@ $(OBJS)  # the order "$(CXXARFLAGS) $@" is required
	$(RANLIB) $@
