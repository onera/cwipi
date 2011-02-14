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

AM_CPPFLAGS = -I@PYTHON_INCLUDE@ -I@NUMPY_INCLUDE@ -I$(top_srcdir)/src $(DEBUG_CPPFLAGS) @BFT_CPPFLAGS@ @FVM_CPPFLAGS@ @MPI_CPPFLAGS@
AM_CFLAGS = -fno-strict-aliasing 
AM_LDFLAGS = @BFT_LDFLAGS@ @FVM_LDFLAGS@ @MPI_LDFLAGS@

pkgpython_PYTHON = __init__.py
pkgpyexec_LTLIBRARIES = CWIPI.la
CWIPI_la_SOURCES = CWIPI.c
CWIPI_la_LIBADD = $(top_builddir)/src/libcwipi.la @BFT_LIBS@ @FVM_LIBS@
CWIPI_la_LDFLAGS = @BFT_LDFLAGS@ @FVM_LDFLAGS@ _avoid_version -module

clean-local:
	rm -f core *.o *.bak *~ *.pyc *.pyo



