#-------------------------------------------------------------------------------
#   This file is part of the "Finite Volume Mesh" library, intended to provide
#   finite volume mesh and associated fields I/O and manipulation services.
#
#   Copyright (C) 2004-2008  EDF
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# To disable dependency tracking
# AUTOMAKE_OPTIONS = no-dependencies

#include directories

#INCLUDES = -I$(top_srcdir)/include

# Add any debug flags

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

# Main part

AM_CPPFLAGS = $(DEBUG_CPPFLAGS) -I$(top_srcdir)/src/bft @MPI_CPPFLAGS@
AM_CFLAGS =
AM_LDFLAGS = @HDF5_LDFLAGS@ @MED_LDFLAGS@ @CGNS_LDFLAGS@ \
             @MPI_LDFLAGS@

# Conditionally compiled extensions

if HAVE_CGNS
libfvmc_cgns_la_CPPFLAGS = $(DEBUG_CPPFLAGS) @MPI_CPPFLAGS@ \
                         @CGNS_CPPFLAGS@
endif
if HAVE_MED
libfvmc_med_la_CPPFLAGS = $(DEBUG_CPPFLAGS)  @MPI_CPPFLAGS@ \
                         @HDF5_CPPFLAGS@ @MED_CPPFLAGS@ 
endif

# Public header files (to be installed)

fvmincludedir = $(includedir)
fvminclude_HEADERS =
#fvminclude_HEADERS = \
#fvmc_block_to_part.h \
#fvmc_coupling.h \
#fvmc_defs.h \
#fvmc_file.h \
#fvmc_group.h \
#fvmc_interface.h \
#fvmc_io_num.h \
#fvmc_locator.h \
#fvmc_neighborhood.h \
#fvmc_nodal.h \
#fvmc_nodal_append.h \
#fvmc_nodal_extrude.h \
#fvmc_nodal_from_desc.h \
#fvmc_nodal_order.h \
#fvmc_nodal_project.h \
#fvmc_nodal_triangulate.h \
#fvmc_nodal_extract.h \
#fvmc_order.h \
#fvmc_parall.h \
#fvmc_part_to_block.h \
#fvmc_periodicity.h \
#fvmc_selector.h \
#fvmc_triangulate.h \
#fvmc_point_location.h \
#fvmc_config_defs.h \
#fvmc_writer.h

noinst_HEADERS = \
fvmc_block_dist.h \
fvmc_box.h \
fvmc_box_priv.h \
fvmc_box_tree.h \
fvmc_convert_array.h \
fvmc_gather.h \
fvmc_morton.h \
fvmc_nodal_priv.h \
fvmc_selector_postfix.h \
fvmc_tesselation.h \
fvmc_to_cgns.h \
fvmc_to_med.h \
fvmc_to_ensight.h \
fvmc_to_ensight_case.h \
fvmc_to_text.h \
fvmc_writer_helper.h \
fvmc_writer_priv.h \
fvmc_block_to_part.h \
fvmc_coupling.h \
fvmc_defs.h \
fvmc_file.h \
fvmc_group.h \
fvmc_interface.h \
fvmc_io_num.h \
fvmc_locator.h \
fvmc_neighborhood.h \
fvmc_nodal.h \
fvmc_nodal_append.h \
fvmc_nodal_extrude.h \
fvmc_nodal_from_desc.h \
fvmc_nodal_order.h \
fvmc_nodal_project.h \
fvmc_nodal_triangulate.h \
fvmc_nodal_extract.h \
fvmc_order.h \
fvmc_parall.h \
fvmc_part_to_block.h \
fvmc_periodicity.h \
fvmc_selector.h \
fvmc_triangulate.h \
fvmc_point_location.h \
fvmc_config_defs.h \
fvmc_writer.h

# Library source files

if HAVE_MPI
lib_LTLIBRARIES = libfvmc.la
else
lib_LTLIBRARIES = libfvmc_scal.la
endif

noinst_LTLIBRARIES =

if HAVE_MPI
libfvmc_la_LIBADD = $(top_builddir)/src/bft/libbftc.la
else
libfvmc_scal_la_LIBADD = $(top_builddir)/src/bft/libbftc.la
endif

if HAVE_MPI
libfvmc_la_SOURCES = \
fvmc_block_dist.c \
fvmc_block_to_part.c \
fvmc_box.c \
fvmc_box_tree.c \
fvmc_convert_array.c \
fvmc_defs.c \
fvmc_file.c \
fvmc_gather.c \
fvmc_group.c \
fvmc_interface.c \
fvmc_io_num.c \
fvmc_locator.c \
fvmc_morton.c \
fvmc_neighborhood.c \
fvmc_nodal.c \
fvmc_nodal_append.c \
fvmc_nodal_extrude.c \
fvmc_nodal_from_desc.c \
fvmc_nodal_order.c \
fvmc_nodal_project.c \
fvmc_nodal_triangulate.c \
fvmc_nodal_extract.c \
fvmc_order.c \
fvmc_parall.c \
fvmc_part_to_block.c \
fvmc_periodicity.c \
fvmc_point_location.c \
fvmc_selector.c \
fvmc_selector_postfix.c \
fvmc_tesselation.c \
fvmc_to_ensight.c \
fvmc_to_ensight_case.c \
fvmc_to_text.c \
fvmc_triangulate.c \
fvmc_writer.c \
fvmc_writer_helper.c
else
libfvmc_scal_la_SOURCES = \
fvmc_block_dist.c \
fvmc_block_to_part.c \
fvmc_box.c \
fvmc_box_tree.c \
fvmc_convert_array.c \
fvmc_defs.c \
fvmc_file.c \
fvmc_gather.c \
fvmc_group.c \
fvmc_interface.c \
fvmc_io_num.c \
fvmc_locator.c \
fvmc_morton.c \
fvmc_neighborhood.c \
fvmc_nodal.c \
fvmc_nodal_append.c \
fvmc_nodal_extrude.c \
fvmc_nodal_from_desc.c \
fvmc_nodal_order.c \
fvmc_nodal_project.c \
fvmc_nodal_triangulate.c \
fvmc_nodal_extract.c \
fvmc_order.c \
fvmc_parall.c \
fvmc_part_to_block.c \
fvmc_periodicity.c \
fvmc_point_location.c \
fvmc_selector.c \
fvmc_selector_postfix.c \
fvmc_tesselation.c \
fvmc_to_ensight.c \
fvmc_to_ensight_case.c \
fvmc_to_text.c \
fvmc_triangulate.c \
fvmc_writer.c \
fvmc_writer_helper.c
endif

#if HAVE_MPI
#lib_LTLIBRARIES += libfvmc_coupl.la
#libfvmc_coupl_la_SOURCES = fvmc_coupling.c
#endif

if HAVE_CGNS
noinst_LTLIBRARIES += libfvmc_cgns.la
if HAVE_MPI
libfvmc_la_LIBADD += libfvmc_cgns.la
else
libfvmc_scal_la_LIBADD += libfvmc_cgns.la
endif
libfvmc_cgns_la_SOURCES = fvmc_to_cgns.c
libfvmc_cgns_la_LIBADD = $(top_builddir)/src/bft/libbftc.la
endif

if HAVE_MED
noinst_LTLIBRARIES += libfvmc_med.la
if HAVE_MPI
libfvmc_la_LIBADD += libfvmc_med.la $(top_builddir)/src/bft/libbfto.la
else
libfvmc_scal_la_LIBADD += libfvmc_med.la $(top_builddir)/src/bft/libbfto.la
endif
libfvmc_med_la_SOURCES = fvmc_to_med.c
libfvmc_med_la_LIBADD = $(top_builddir)/src/bft/libbftc.la
endif

if HAVE_MPI
libfvmc_la_LDFLAGS = -no-undefined -version-info 0:0:0
else
libfvmc_scal_la_LDFLAGS = -no-undefined -version-info 0:0:0
endif
#libfvmc_coupl_la_LDFLAGS = -no-undefined -version-info 0:0:0

