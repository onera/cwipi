#------------------------------------------------------------------------------
#   This file is part of the "Mathematical Expression Interpreter" library.
#
#   Copyright (C) 2010-2011  ONERA
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
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#	$(PYTHON) setup.py build_ext --include-dirs=/opt/intel/impi/3.2.1.009/include64:/home/equemera/opt/bft-1.1.3/Linux_x86_64/include:/home/equemera/opt/fvm-0.15.1-onera-5/Linux_x86_64/include:/home/equemera/opt/cwipi-0.2.10/Linux_x86_64/include \
#	--library-dirs=/opt/intel/impi/3.2.1.009/lib64:/home/equemera/opt/bft-1.1.3/Linux_x86_64/lib:/home/equemera/opt/fvm-0.15.1-onera-5/Linux_x86_64/lib:/home/equemera/opt/cwipi-0.2.10/Linux_x86_64/lib \
#	--rpath=/opt/intel/impi/3.2.1.009/lib64:/home/equemera/opt/bft-1.1.3/Linux_x86_64/lib:/home/equemera/opt/fvm-0.15.1-onera-5/Linux_x86_64/lib:/home/equemera/opt/cwipi-0.2.10/Linux_x86_64/lib ;\
#------------------------------------------------------------------------------

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



