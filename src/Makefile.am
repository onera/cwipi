
# Add any debug flags

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

AM_CPPFLAGS = $(DEBUG_CPPFLAGS) @BFT_CPPFLAGS@ @FVM_CPPFLAGS@ @MPI_CPPFLAGS@
AM_CFLAGS =
AM_LDFLAGS = @BFT_LDFLAGS@ @FVM_LDFLAGS@ @MPI_LDFLAGS@

# Public header files (to be installed)

couplingsincludedir = $(includedir)
couplingsinclude_HEADERS = \
couplings_f.h \
couplings.h

noinst_HEADERS = \
coo_baryc.h \
couplings_cf.h \
quickSort.h \
solve_ax_b_4.h

# Library source files

lib_LTLIBRARIES = libcouplings.la libcouplings_f.la
noinst_LTLIBRARIES =
libcouplings_la_LIBADD =

libcouplings_la_SOURCES = \
applicationProperties.cxx \
applicationPropertiesDataBase.cxx \
couplings.cxx \
quickSort.c \
solve_ax_b_4.c \
coo_baryc.c \
mesh.cxx \
coupling.cxx \
couplingDataBase.cxx

libcouplings_f_la_LIBADD = libcouplings.la

libcouplings_f_la_SOURCES = \
couplings_cf.cxx \
printfort.f \
couplings_f.f


