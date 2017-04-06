
# Add any debug flags

SUBDIRS = bft fvm .

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

AM_CPPFLAGS = $(DEBUG_CPPFLAGS) -I$(top_srcdir)/src/bft  -I$(top_srcdir)/src/fvm @MPI_CPPFLAGS@
AM_CFLAGS =
AM_LDFLAGS = @MPI_LDFLAGS@

# Public header files (to be installed)

cwipiincludedir = $(includedir)
cwipiinclude_HEADERS = \
cwipi.h 

if HAVE_FORTRAN
 cwipiinclude_HEADERS += cwipi.mod cwipi_f.f90
endif

noinst_HEADERS = \
cwipi_cf.h \
quickSort.h \
cwipi_timer.h \
solve_ax_b_4.h

# Library source files

if HAVE_FORTRAN
  lib_LTLIBRARIES = libcwipi.la libcwipi_f.la
else
  lib_LTLIBRARIES = libcwipi.la
endif

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
mesh.cxx \
conservativeMesh.cxx \
coupling.cxx \
couplingDataBase.cxx \
geomUtilities.cxx \
cwipi_timer.c \
creeMaillagePolygone2D.c \
grid_mesh.c

#cwipi_timer.c $(FORTRAN_IN_C)

if HAVE_FORTRAN
  libcwipi_f_la_LIBADD = \
  $(top_builddir)/src/bft/libbftc.la \
  $(top_builddir)/src/fvm/libfvmc.la \
  libcwipi.la \
  @MPI_LIBS@

  libcwipi_f_la_SOURCES = \
  cwipi_printfort.f90 \
  callfortinterpfct.f90 \
  printfort.f90 \
  cwipi_printfort_cf.cxx \
  cwipi_set_output_logical_unit.f90 \
  cwipi_cf.cxx \
  cwipi_f.f90
endif

