#-----------------------------------------------------------------------------
# This file is part of the CWIPI library. 
#
# Copyright (C) 2011  ONERA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------
# To disable dependency tracking
# AUTOMAKE_OPTIONS = no-dependencies

AUTOMAKE_OPTIONS = foreign

# Add any debug flags

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

# Main part

AM_CPPFLAGS = \
-I$(top_builddir)/src \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
$(DEBUG_CPPFLAGS) \
@MPI_CPPFLAGS@

AM_CFLAGS   =
AM_LDFLAGS  = @MPI_LDFLAGS@

AM_FCFLAGS = -I$(top_builddir)/src -I$(top_srcdir)/src

LDADD = $(top_builddir)/src/bft/libbfto.la \
	$(top_builddir)/src/fvm/libfvmo.la \
	$(top_builddir)/src/libcwipi.la \
        @MPI_LIBS@ \
	-lm

check_PROGRAMS = \
c_linear_coupling \
c_surf_coupling_P1P1_micron_polygon \
c_surf_coupling_P1P1_polygon \
c_surf_coupling_P1P1 \
c_surf_coupling_P1P1_asynchronous \
c_mean_value_3D
#c_surf_coupling_P1P0_P0P1 \
#test3D_1_c1 \
#test3D_1_c2 

#test2D_1_c1 
#test2D_1_c2 

c_linear_coupling_SOURCES = c_linear_coupling.c
c_surf_coupling_P1P1_micron_polygon_SOURCES = creeMaillagePolygone2D.c c_surf_coupling_P1P1_micron_polygon.c
c_surf_coupling_P1P1_polygon_SOURCES = creeMaillagePolygone2D.c c_surf_coupling_P1P1_polygon.c
c_surf_coupling_P1P1_SOURCES = grid_mesh.c c_surf_coupling_P1P1.c
c_surf_coupling_P1P1_asynchronous_SOURCES = grid_mesh.c c_surf_coupling_P1P1_asynchronous.c
c_mean_value_3D_SOURCES = c_mean_value_3D.c
#test2D_1_c1_SOURCES = creeMaillagePolygone2D.c test2D_1_c1.c
#test2D_1_c2_SOURCES = creeMaillagePolygone2D.c test2D_1_c2.f90
#test3D_1_c1_SOURCES = test3D_1_c1.c
#test3D_1_c2_SOURCES = test3D_1_c2.c

c_mean_value_3D_LDADD = \
$(top_builddir)/src/bft/libbfto.la \
$(top_builddir)/src/fvm/libfvmo.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_linear_coupling_LDADD = \
$(top_builddir)/src/bft/libbfto.la \
$(top_builddir)/src/fvm/libfvmo.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_P1P1_micron_polygon_LDADD = \
$(top_builddir)/src/bft/libbfto.la \
$(top_builddir)/src/fvm/libfvmo.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_P1P1_polygon_LDADD = \
$(top_builddir)/src/bft/libbfto.la \
$(top_builddir)/src/fvm/libfvmo.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_P1P1_LDADD = \
$(top_builddir)/src/bft/libbfto.la \
$(top_builddir)/src/fvm/libfvmo.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

# test3D_1_c1_LDADD = \
# $(top_builddir)/src/bft/libbft.la \
# $(top_builddir)/src/fvm/libfvm.la \
# $(top_builddir)/src/libcwipi.la \
# @MPI_LIBS@ \
# -lm

# test3D_1_c2_LDADD = \
# $(top_builddir)/src/bft/libbft.la \
# $(top_builddir)/src/fvm/libfvm.la \
# $(top_builddir)/src/libcwipi.la \
# @MPI_LIBS@ \
# -lm

#test2D_1_c1_LDADD =  $(top_builddir)/src/libcwipi.la \
#        @FVM_LIBS@ @BFT_LIBS@ @MPI_LIBS@ \
#	-lm
#test2D_1_c2_LDADD =  $(top_builddir)/src/libcwipi.la \
#        $(top_builddir)/src/libcwipi_f.la \
#        @FVM_LIBS@ @BFT_LIBS@ @MPI_LIBS@ \
#	-lm

VALGRIND =
MPI_EXEC = mpiexec
MPI_NP = -np
NP = 2

TESTS_ENVIRONMENT = ${MPI_EXEC} ${MPI_NP} ${NP} ${VALGRIND}
TESTS = c_linear_coupling \
	c_surf_coupling_P1P1_micron_polygon \
	c_surf_coupling_P1P1_polygon \
	c_surf_coupling_P1P1 \
	c_surf_coupling_P1P1_asynchronous \
	c_mean_value_3D
if ENABLE_PYTHON 
  TESTS += python_api.py
endif

