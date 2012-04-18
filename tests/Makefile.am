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

c_vol_coupling_P1P1.o: copy-data

copy-data:
	@cp -r $(top_srcdir)/tests/meshes $(top_builddir)/tests

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

LDADD = $(top_builddir)/src/bft/libbftc.la \
	$(top_builddir)/src/fvm/libfvmc.la \
	$(top_builddir)/src/libcwipi.la \
        @MPI_LIBS@ \
	-lm

check_PROGRAMS = \
c_linear_coupling \
c_surf_coupling_P1P1_micron_polygon \
c_surf_coupling_P1P1_polygon \
c_surf_coupling_P1P1 \
c_surf_coupling_P1P1_asynchronous \
c_surf_coupling_P1P0_P0P1 \
c_surf_coupling_user_target_points \
c_surf_coupling_user_interpolation \
c_surf_coupling_user_interpolation2 \
c_surf_coupling_vector_exchange \
c_surf_coupling_output_error \
c_surf_coupling_simple_location \
c_surf_without_partitioning_coupling \
c_surf_without_partitioning_coupling2 \
c_surf_sequential_coupling \
c_vol_coupling_P1P1 \
c_vol_coupling_polyhedra_P1P1 \
c_vol_coupling_micro_polyhedra_P1P1 \
fortran_surf_P1P1

c_linear_coupling_SOURCES = c_linear_coupling.c
c_surf_coupling_P1P1_micron_polygon_SOURCES = creeMaillagePolygone2D.c c_surf_coupling_P1P1_micron_polygon.c
c_surf_coupling_P1P1_polygon_SOURCES = creeMaillagePolygone2D.c c_surf_coupling_P1P1_polygon.c
c_surf_coupling_P1P1_SOURCES = grid_mesh.c c_surf_coupling_P1P1.c
c_surf_coupling_P1P0_P0P1_SOURCES = grid_mesh.c c_surf_coupling_P1P0_P0P1.c
c_surf_coupling_P1P1_asynchronous_SOURCES = grid_mesh.c c_surf_coupling_P1P1_asynchronous.c
c_surf_coupling_user_target_points_SOURCES = grid_mesh.c c_surf_coupling_user_target_points.c
c_surf_coupling_user_interpolation_SOURCES = grid_mesh.c c_surf_coupling_user_interpolation.c
c_surf_coupling_user_interpolation2_SOURCES = grid_mesh.c c_surf_coupling_user_interpolation2.c
c_surf_coupling_vector_exchange_SOURCES = grid_mesh.c c_surf_coupling_vector_exchange.c
c_surf_coupling_output_error_SOURCES = grid_mesh.c c_surf_coupling_output_error.c
c_surf_coupling_simple_location_SOURCES = grid_mesh.c c_surf_coupling_simple_location.c
c_surf_without_partitioning_coupling_SOURCES = grid_mesh.c c_surf_without_partitioning_coupling.c
c_surf_without_partitioning_coupling2_SOURCES = grid_mesh.c c_surf_without_partitioning_coupling2.c
c_surf_sequential_coupling_SOURCES = grid_mesh.c c_surf_sequential_coupling.c
c_vol_coupling_P1P1_SOURCES = c_vol_coupling_P1P1.c
c_vol_coupling_micro_polyhedra_P1P1_SOURCES = c_vol_coupling_micro_polyhedra_P1P1.c
c_vol_coupling_polyhedra_P1P1_SOURCES = c_vol_coupling_polyhedra_P1P1.c
fortran_surf_P1P1_SOURCES = grid_mesh.c fortran_surf_P1P1.f90

c_linear_coupling_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_P1P1_micron_polygon_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_P1P1_polygon_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_P1P1_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_P1P0_P0P1_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_user_target_points_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_user_interpolation_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_user_interpolation2_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_vector_exchange_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_output_error_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_coupling_simple_location_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_without_partitioning_coupling_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_without_partitioning_coupling2_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_surf_sequential_coupling_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_vol_coupling_P1P1_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_vol_coupling_polyhedra_P1P1_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

c_vol_coupling_micro_polyhedra_P1P1_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
@MPI_LIBS@ \
-lm

fortran_surf_P1P1_LDADD = \
$(top_builddir)/src/bft/libbftc.la \
$(top_builddir)/src/fvm/libfvmc.la \
$(top_builddir)/src/libcwipi.la \
$(top_builddir)/src/libcwipi_f.la \
 @MPI_LIBS@ \
 -lm

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
	c_surf_coupling_P1P0_P0P1 \
	c_surf_coupling_user_target_points \
	c_surf_coupling_user_interpolation \
	c_surf_coupling_user_interpolation2 \
	c_surf_coupling_output_error \
	c_surf_coupling_simple_location \
	c_surf_without_partitioning_coupling \
	c_surf_without_partitioning_coupling2 \
	c_surf_sequential_coupling \
	c_vol_coupling_P1P1 \
	c_vol_coupling_polyhedra_P1P1 \
	c_vol_coupling_micro_polyhedra_P1P1 \
	fortran_surf_P1P1

if ENABLE_PYTHON 
  TESTS += python_api.py
endif

