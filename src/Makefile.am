
# Add any debug flags

SUBDIRS = bft fvm . fortran

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

AM_CPPFLAGS = $(DEBUG_CPPFLAGS) -I$(top_srcdir)/src/bft  -I$(top_srcdir)/src/fortran -I$(top_srcdir)/src/fvm @MPI_CPPFLAGS@
AM_CFLAGS =
AM_LDFLAGS = @MPI_LDFLAGS@

# Public header files (to be installed)

cwipiincludedir = $(includedir)
cwipiinclude_HEADERS = \
cwipi.h 

noinst_HEADERS = \
quickSort.h \
cwipi_timer.h \
solve_ax_b_4.h

# Library source files

lib_LTLIBRARIES = libcwipi.la

noinst_LTLIBRARIES =

libcwipi_la_LIBADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
@MPI_LIBS@

libcwipi_la_SOURCES = applicationProperties.cxx \
applicationPropertiesDataBase.cxx \
locationToLocalMesh.cxx \
locationToDistantMesh.cxx \
cwipi.cxx \
quickSort.c \
solve_ax_b_4.c \
oldMesh.cxx \
oldCoupling.cxx \
couplingDataBase.cxx \
geomUtilities.cxx \
cwipi_timer.c \
creeMaillagePolygone2D.c \
grid_mesh.c

#cwipi_timer.c $(FORTRAN_IN_C)


