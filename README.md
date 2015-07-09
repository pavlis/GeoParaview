# GeoParaview
This is a colletion of programs and libraries to build georeferenced files that can be visualized in 3D with the open source 
package called paraview.  This is necessary because paraview has limited supported for georeferencing data.  When I started
developing this package it had none.   

The package uses a concept encapsulated in an object found in libgeocoords with the name RegionalCoordinates.  A 
RegionalCoordinates object can be use to convert between geographic coordinates and a local cartesian coordinate
system defined by three parameters:  origin latitude, origin longitude, origin radius, and a local azimuth to (if 
desired) rotate the coordinate system.   The later can be useful, for example, when local geography for a study area has
an orientation that is not cardinal.   An example is a coastline. 

The package has library with what I hope are some useful C++ objects that simplify certain geographic manipulations.
There is a collection of programs with very different themes.  There are paraview vtk converters that are the primary 
content of this collection.  There is also a modeling directory that has some (poorly document) research code that
uses libgeocoords to build some useful plate tectonics oriented 3d models for use with paraview.
