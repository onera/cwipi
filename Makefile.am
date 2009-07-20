AUX_DIST = $(ac_aux_dir)/config.guess \
           $(ac_aux_dir)/config.sub \
           $(ac_aux_dir)/depcomp \
           $(ac_aux_dir)/install-sh \
           $(ac_aux_dir)/missing \
           $(ac_aux_dir)/mkinstalldirs \
           $(ac_aux_dir)/ltmain.sh

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure config-h.in stamp-h.in $(AUX_DIST)

#SUBDIRS = src tests doc

SUBDIRS = src tests 

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

install-tests:
	@cd tests ; $(MAKE) install-tests

uninstall-tests:
	@cd tests ; $(MAKE) uninstall-tests
