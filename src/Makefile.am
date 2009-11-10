
# Add any debug flags

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

AM_CPPFLAGS = $(DEBUG_CPPFLAGS) @BFT_CPPFLAGS@ @FVM_CPPFLAGS@ @MPI_CPPFLAGS@
AM_CFLAGS =
AM_LDFLAGS = @BFT_LDFLAGS@ @FVM_LDFLAGS@ @MPI_LDFLAGS@

#LTFCCOMPILE = $(LIBTOOL) --mode=compile --tag=FC $(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_FCFLAGS) $(FCFLAGS) $(FCFLAGS_f90)
#FCCOMPILE = $(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(FCFLAGS_f90)
# Public header files (to be installed)

cwipiincludedir = $(includedir)
cwipiinclude_HEADERS = \
cwipi.mod \
cwipi.h 

noinst_HEADERS = \
coo_baryc.h \
cwipi_cf.h \
quickSort.h \
solve_ax_b_4.h

# Library source files

lib_LTLIBRARIES = libcwipi.la libcwipi_f.la
noinst_LTLIBRARIES =
libcwipi_la_LIBADD =

libcwipi_la_SOURCES = \
applicationProperties.cxx \
applicationPropertiesDataBase.cxx \
localLocation.cxx \
cwipi.cxx \
quickSort.c \
solve_ax_b_4.c \
coo_baryc.c \
mesh.cxx \
coupling.cxx \
couplingDataBase.cxx \
callfortinterpfct.f90

libcwipi_f_la_LIBADD = libcwipi.la

libcwipi_f_la_SOURCES = \
cwipi_cf.cxx \
cwipi_f.f90 \
printfort.f90


