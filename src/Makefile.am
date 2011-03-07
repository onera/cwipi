
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
cwipi.mod \
cwipi_f.f90 \
cwipi.h 

noinst_HEADERS = \
cwipi_cf.h \
quickSort.h \
solve_ax_b_4.h

# Library source files

lib_LTLIBRARIES = libcwipi.la libcwipi_f.la
noinst_LTLIBRARIES =
libcwipi_la_LIBADD = \
$(top_builddir)/src/bft/libbfto.la \
$(top_builddir)/src/fvm/libfvmo.la

libcwipi_la_SOURCES = \
applicationProperties.cxx \
applicationPropertiesDataBase.cxx \
locationToLocalMesh.cxx \
locationToDistantMesh.cxx \
cwipi.cxx \
quickSort.c \
solve_ax_b_4.c \
mesh.cxx \
coupling.cxx \
couplingDataBase.cxx \
callfortinterpfct.f90

libcwipi_f_la_LIBADD = \
$(top_builddir)/src/bft/libbfto.la \
$(top_builddir)/src/fvm/libfvmo.la \
libcwipi.la

libcwipi_f_la_SOURCES = \
cwipi_cf.cxx \
cwipi_f.f90 \
printfort.f90


