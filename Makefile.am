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

AUX_DIST = $(ac_aux_dir)/config.guess \
           $(ac_aux_dir)/config.sub \
           $(ac_aux_dir)/depcomp \
           $(ac_aux_dir)/install-sh \
           $(ac_aux_dir)/missing \
           $(ac_aux_dir)/mkinstalldirs \
           $(ac_aux_dir)/ltmain.sh

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure config-h.in stamp-h.in $(AUX_DIST)

#SUBDIRS = src tests doc

if ENABLE_PYTHON 
 SUBDIRS = src cython tests
else
 SUBDIRS = src tests
endif

EXTRA_DIST = $(ac_aux_dir)/cwipi_auto_flags.sh \
             $(ac_aux_dir)/check_bft.m4 \
             $(ac_aux_dir)/check_fvm.m4 \
             $(ac_aux_dir)/cwipi_config_info.m4 \
             $(ac_aux_dir)/cwipi_config_publ.m4 \
             $(ac_aux_dir)/check_mpi.m4

bin_SCRIPTS = cwipi-config

configincludedir = $(includedir)
configinclude_HEADERS = cwipi_config.h 
#Avoir configinclude_HEADERS = cwipi_config_priv.h cwipi_config.h ???

uninstall-local:
	rmdir $(includedir)/cwipi

DISTCLEANFILES = cwipi-config cwipi_config.h cwipi_config_priv.h

ACLOCAL_AMFLAGS = -I config
