dnl Copyright (C) 2005-2009 EDF
dnl
dnl This file is part of the FVM software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl FVM source distribution.

# FVMC_AC_TEST_HDF5
#----------------
# modifies or sets have_hdf5, HDF_CPPFLAGS, HDF_LDFLAGS, and HDF_LIBS
# depending on libraries found

AC_DEFUN([FVMC_AC_TEST_HDF5], [

AC_REQUIRE([FVMC_AC_CONFIG_PUBL_INIT])dnl

have_hdf5=no

AC_ARG_WITH(hdf5,
            [AS_HELP_STRING([--with-hdf5=PATH],
                            [specify prefix directory for HDF5])],
            [if test "x$withval" = "x"; then
               with_hdf5=yes
             fi],
            [with_hdf5=check])

AC_ARG_WITH(hdf5-include,
            [AS_HELP_STRING([--with-hdf5-include=PATH],
                            [specify directory for HDF5 include files])],
            [if test "x$with_hdf5" = "xcheck"; then
               with_hdf5=yes
             fi
             HDF5_CPPFLAGS="-I$with_hdf5_include"],
            [if test "x$with_hdf5" != "xno" -a "x$with_hdf5" != "xyes" \
	          -a "x$with_hdf5" != "xcheck"; then
               HDF5_CPPFLAGS="-I$with_hdf5/include"
             fi])

AC_ARG_WITH(hdf5-lib,
            [AS_HELP_STRING([--with-hdf5-lib=PATH],
                            [specify directory for HDF5 library])],
            [if test "x$with_hdf5" = "xcheck"; then
               with_hdf5=yes
             fi
             HDF5_LDFLAGS="-L$with_hdf5_lib"],
            [if test "x$with_hdf5" != "xno" -a "x$with_hdf5" != "xyes" \
	          -a "x$with_hdf5" != "xcheck"; then
               HDF5_LDFLAGS="-L$with_hdf5/lib"
             fi])


if test "x$with_hdf5" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  HDF5_LIBS="-lhdf5 $PTHREAD_LIBS"
  
  CPPFLAGS="${CPPFLAGS} ${HDF5_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${HDF5_LDFLAGS}"
  LIBS="${LIBS} ${HDF5_LIBS}"

  AC_CHECK_LIB(hdf5, H5Fopen, 
               [ AC_DEFINE([HAVE_HDF5], 1, [HDF5 file support])
                 have_hdf5=yes
               ], 
               [if test "x$with_hdf5" != "xcheck" ; then
                  AC_MSG_FAILURE([HDF5 support is requested, but test for HDF5 failed!])
                else
                  AC_MSG_WARN([no HDF5 file support])
                fi
               ],
               )

  if test "x$have_hdf5" = "xno"; then
    HDF5_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AC_SUBST(HDF5_CPPFLAGS)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(HDF5_LIBS)

])dnl

