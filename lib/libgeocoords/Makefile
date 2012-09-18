LIB=libgeocoords.a
INCLUDE=GeoCoordError.h \
  GeoPath.h \
  GeoPolygonRegion.h \
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
GeoSplineSurface.cc : GeoSplineSurface.h
GeoTriMeshSurface.cc : GeoTriMeshSurface.h
PLGeoPath.cc : GeoPath.h PLGeoPath.h
PlateBoundaryPath.cc : GeoPath.h PlateBoundaryPath.h
RegionalCoordinates.cc : RegionalCoordinates.h

OBJS=GeoSplineSurface.o \
  GeoTriMeshSurface.o \
  GeoPolygonRegion.o \
  LatLong-UTMconversion.o \
  PLGeoPath.o \
  PlateBoundaryPath.o \
  RegionalCoordinates.o

$(LIB) : $(OBJS)
	$(RM) $@
	$(CXXAR) $(CXXARFLAGS) $@ $(OBJS)  # the order "$(CXXARFLAGS) $@" is required
	$(RANLIB) $@
