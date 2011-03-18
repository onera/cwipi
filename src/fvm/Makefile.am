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

CC = $(CXX)

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

libfvmo_cgns_la_CPPFLAGS = $(DEBUG_CPPFLAGS) @MPI_CPPFLAGS@ \
                         @CGNS_CPPFLAGS@
libfvmo_med_la_CPPFLAGS = $(DEBUG_CPPFLAGS)  @MPI_CPPFLAGS@ \
                         @HDF5_CPPFLAGS@ @MED_CPPFLAGS@ 

# Public header files (to be installed)

fvmincludedir = $(includedir)
fvminclude_HEADERS = \
fvm_block_to_part.h \
fvm_coupling.h \
fvm_defs.h \
fvm_file.h \
fvm_group.h \
fvm_interface.h \
fvm_io_num.h \
fvm_locator.h \
fvm_neighborhood.h \
fvm_nodal.h \
fvm_nodal_append.h \
fvm_nodal_extrude.h \
fvm_nodal_from_desc.h \
fvm_nodal_order.h \
fvm_nodal_project.h \
fvm_nodal_triangulate.h \
fvm_nodal_extract.h \
fvm_order.h \
fvm_parall.h \
fvm_part_to_block.h \
fvm_periodicity.h \
fvm_selector.h \
fvm_triangulate.h \
fvm_point_location.h \
fvm_config_defs.h \
fvm_writer.h

noinst_HEADERS = \
fvm_block_dist.h \
fvm_box.h \
fvm_box_priv.h \
fvm_box_tree.h \
fvm_convert_array.h \
fvm_gather.h \
fvm_morton.h \
fvm_nodal_priv.h \
fvm_selector_postfix.h \
fvm_tesselation.h \
fvm_to_cgns.h \
fvm_to_med.h \
fvm_to_ensight.h \
fvm_to_ensight_case.h \
fvm_to_text.h \
fvm_writer_helper.h \
fvm_writer_priv.h

# Library source files

if HAVE_MPI
lib_LTLIBRARIES = libfvmo.la
else
lib_LTLIBRARIES = libfvmo_scal.la
endif

noinst_LTLIBRARIES =

if HAVE_MPI
libfvmo_la_LIBADD = $(top_builddir)/src/bft/libbfto.la
else
libfvmo_scal_la_LIBADD = $(top_builddir)/src/bft/libbfto.la
endif

if HAVE_MPI
libfvmo_la_SOURCES = \
fvm_block_dist.cpp \
fvm_block_to_part.cpp \
fvm_box.cpp \
fvm_box_tree.cpp \
fvm_convert_array.cpp \
fvm_defs.cpp \
fvm_file.cpp \
fvm_gather.cpp \
fvm_group.cpp \
fvm_interface.cpp \
fvm_io_num.cpp \
fvm_locator.cpp \
fvm_morton.cpp \
fvm_neighborhood.cpp \
fvm_nodal.cpp \
fvm_nodal_append.cpp \
fvm_nodal_extrude.cpp \
fvm_nodal_from_desc.cpp \
fvm_nodal_order.cpp \
fvm_nodal_project.cpp \
fvm_nodal_triangulate.cpp \
fvm_nodal_extract.cpp \
fvm_order.cpp \
fvm_parall.cpp \
fvm_part_to_block.cpp \
fvm_periodicity.cpp \
fvm_point_location.cpp \
fvm_selector.cpp \
fvm_selector_postfix.cpp \
fvm_tesselation.cpp \
fvm_to_ensight.cpp \
fvm_to_ensight_case.cpp \
fvm_to_text.cpp \
fvm_triangulate.cpp \
fvm_writer.cpp \
fvm_writer_helper.cpp
else
libfvmo_scal_la_SOURCES = \
fvm_block_dist.cpp \
fvm_block_to_part.cpp \
fvm_box.cpp \
fvm_box_tree.cpp \
fvm_convert_array.cpp \
fvm_defs.cpp \
fvm_file.cpp \
fvm_gather.cpp \
fvm_group.cpp \
fvm_interface.cpp \
fvm_io_num.cpp \
fvm_locator.cpp \
fvm_morton.cpp \
fvm_neighborhood.cpp \
fvm_nodal.cpp \
fvm_nodal_append.cpp \
fvm_nodal_extrude.cpp \
fvm_nodal_from_desc.cpp \
fvm_nodal_order.cpp \
fvm_nodal_project.cpp \
fvm_nodal_triangulate.cpp \
fvm_nodal_extract.cpp \
fvm_order.cpp \
fvm_parall.cpp \
fvm_part_to_block.cpp \
fvm_periodicity.cpp \
fvm_point_location.cpp \
fvm_selector.cpp \
fvm_selector_postfix.cpp \
fvm_tesselation.cpp \
fvm_to_ensight.cpp \
fvm_to_ensight_case.cpp \
fvm_to_text.cpp \
fvm_triangulate.cpp \
fvm_writer.cpp \
fvm_writer_helper.cpp
endif

if HAVE_MPI
lib_LTLIBRARIES += libfvmo_coupl.la
libfvmo_coupl_la_SOURCES = fvm_coupling.cpp
endif

if HAVE_CGNS
noinst_LTLIBRARIES += libfvmo_cgns.la
if HAVE_MPI
libfvmo_la_LIBADD += libfvmo_cgns.la $(top_builddir)/src/bft/libbfto.la
else
libfvmo_scal_la_LIBADD += libfvmo_cgns.la $(top_builddir)/src/bft/libbfto.la
endif
libfvmo_cgns_la_SOURCES = fvm_to_cgns.cpp
endif

if HAVE_MED
noinst_LTLIBRARIES += libfvmo_med.la
if HAVE_MPI
libfvmo_la_LIBADD += libfvmo_med.la $(top_builddir)/src/bft/libbfto.la
else
libfvmo_scal_la_LIBADD += libfvmo_med.la $(top_builddir)/src/bft/libbfto.la
endif
libfvmo_med_la_SOURCES = fvm_to_med.cpp
endif

if HAVE_MPI
libfvmo_la_LDFLAGS = -no-undefined -version-info 0:0:0
else
libfvmo_scal_la_LDFLAGS = -no-undefined -version-info 0:0:0
endif
libfvmo_coupl_la_LDFLAGS = -no-undefined -version-info 0:0:0

