#-------------------------------------------------------------------------------
#   This file is part of the "Finite Volume Mesh" library, intended to provide
#   finite volume mesh and associated fields I/O and manipulation services.
#
#   Copyright (C) 2004-2008  EDF
#

# To disable dependency tracking
# AUTOMAKE_OPTIONS = no-dependencies

AUTOMAKE_OPTIONS = foreign

#include directories

#INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/src 
INCLUDES = -I$(includedir)

# Add any debug flags

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

# Main part

AM_CPPFLAGS = $(DEBUG_CPPFLAGS) @BFT_CPPFLAGS@ @FVM_CPPFLAGS@ @MPI_CPPFLAGS@ @METIS_CPPFLAGS@
AM_CFLAGS   =
AM_LDFLAGS  = @BFT_LDFLAGS@ @FVM_LDFLAGS@ @MPI_LDFLAGS@ @METIS_LDFLAGS@

LTFCCOMPILE = $(LIBTOOL) --mode=compile --tag=FC $(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(FCFLAGS_f90)
FCCOMPILE = $(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(FCFLAGS_f90)

LDADD = $(top_builddir)/src/libcouplings.la \
        @BFT_LIBS@ @FVM_LIBS@ @MPI_LIBS@ \
	-lm

testsdir = $(bindir)/cwipi_tests

tests_PROGRAMS = \
test2D_1_c1 \
test2D_1_c2 \
test3D_1_c1 \
test3D_1_c2 


test2D_1_c1_SOURCES = creeMaillagePolygone2D.c test2D_1_c1.c
test2D_1_c2_SOURCES = creeMaillagePolygone2D.c test2D_1_c2.f90
test3D_1_c1_SOURCES = test3D_1_c1.c
test3D_1_c2_SOURCES = test3D_1_c2.c


test2D_1_c1_LDADD = $(top_builddir)/src/libcouplings.la \
        @BFT_LIBS@ @FVM_LIBS@ @MPI_LIBS@ @METIS_LIBS@\
	-lm
test2D_1_c2_LDADD = $(top_builddir)/src/libcouplings.la \
        $(top_builddir)/src/libcouplings_f.la \
        @BFT_LIBS@ @FVM_LIBS@ @MPI_LIBS@ @METIS_LIBS@\
	-lm
test3D_1_c1_LDADD = $(top_builddir)/src/libcouplings.la \
        @BFT_LIBS@ @FVM_LIBS@ @MPI_LIBS@ @METIS_LIBS@\
	-lm
test3D_1_c2_LDADD = $(top_builddir)/src/libcouplings.la \
        @BFT_LIBS@ @FVM_LIBS@ @MPI_LIBS@ @METIS_LIBS@\
	-lm

# To link without libtool, replace libfvm.la above by .libs/libfvm.so
# and uncomment the following line:
# fvm_test_LINK    = $(CC) -o $@

install:
	@echo "Test programs not installed by default;" \
	 echo "use make install-test if preferred."

install-tests: install-testsPROGRAMS
	@cp $(top_srcdir)/tests/Readme $(testsdir)/
	@cp -r $(top_srcdir)/tests/meshes $(testsdir)/

uninstall-tests: uninstall-testsPROGRAMS

