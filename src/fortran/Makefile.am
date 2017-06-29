#-------------------------------------------------------------------------------
#
#   Copyright (C) 2017       ONERA
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

# Recommended for VPATH support in autoconf documentation ; does not
# seem necessary with current Makefile functionality using
# autoconf 2.59 and automake 1.8.
# srcdir = @srcdir@
# VPATH = @srcdir@

# To disable dependency tracking
# AUTOMAKE_OPTIONS = no-dependencies

#include directories

# Add any debug flags

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

# Main part

AM_CPPFLAGS = $(DEBUG_CPPFLAGS) -I$(top_srcdir)/src/bft -I$(top_srcdir)/src/fvm -I$(top_srcdir)/src/ -I.
AM_CFLAGS   =
AM_LDFLAGS  =

# Public header files (to be installed)

cwipi_fincludedir = $(includedir) ../bft
if HAVE_FORTRAN
 cwipi_finclude_HEADERS = cwipi.mod cwipi_f.f90
endif

noinst_HEADERS = \
cwipi_cf.h

if HAVE_FORTRAN
  lib_LTLIBRARIES = libcwipi_f.la

  libcwipi_f_la_LIBADD = \
  $(top_builddir)/src/bft/libbftc.la \
  $(top_builddir)/src/fvm/libfvmc.la \
  $(top_builddir)/src/libcwipi.la \
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
