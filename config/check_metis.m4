dnl  Copyright (C) 2009 ONERA
dnl
AC_DEFUN([CHECK_METIS],[

AC_CHECKING(for METIS Library)

AC_LANG_SAVE
AC_LANG_C

METIS_CPPFLAGS=""
METIS_LDFLAGS=""
METIS_LIBS=""

have_metis=no
metis_include_dir=""
metis_lib_dir=""

AC_CHECKING(for METIS location)
AC_ARG_WITH(metis,
            [  --with-metis=PATH          specify prefix directory for METIS]
)

AC_ARG_WITH(metis-include,
            [  --with-metis-include=PATH  specify directory for METIS include files]
)

AC_ARG_WITH(metis-lib,
            [  --with-metis-lib=PATH      specify directory for METIS library]
)

if test "x$with_metis_include" != "x" ; then
  METIS_CPPFLAGS="-I$with_metis_include"
  metis_include_dir=$with_metis_include
elif test "x$with_metis" != "x" ; then
  METIS_CPPFLAGS="-I$with_metis/Lib"
  metis_include_dir=$with_metis/Lib
fi

if test "x$with_metis_lib" != "x" ; then
  METIS_LDFLAGS="-L$with_metis_lib"
  metis_lib_dir=$with_metis_lib
elif test "x$with_metis" != "x" ; then
  METIS_LDFLAGS="-L$with_metis"
  metis_lib_dir=$with_metis           
fi

METIS_LIBS="-lmetis"

saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
saved_LIBS="$LIBS"

CPPFLAGS="${CPPFLAGS} ${METIS_CPPFLAGS}"
LDFLAGS="${LDFLAGS} ${METIS_LDFLAGS}"
LIBS="${LIBS} ${METIS_LIBS}"

dnl METIS headers
AC_CHECKING(for METIS headers)

if test "x${cross_compiling}" = "xno" ; then
 AC_CHECK_FILE(${metis_include_dir}/metis.h,
                have_metis=yes,
                have_metis=no)
fi

if test "x${have_metis}" = "xyes" ; then

  AC_TRY_COMPILE([#include <metis.h>],
                 [Change2CNumbering(0,0,0)],
                 have_metis=yes,
                 have_metis=no)
fi

AC_MSG_RESULT(for metis headers: $have_metis)

dnl METIS binaries
AC_CHECKING(for METIS binaries)

if test "x${cross_compiling}" = "xno" ; then
if test "x$have_metis" = "xyes" ; then

  AC_CHECK_FILE(${metis_lib_dir}/libmetis.a,
                  have_metis=yes,
                  have_metis=no)
fi
fi

if test "x$have_metis" = "xyes" ; then
  AC_TRY_LINK([#include <metis.h>],
              [Change2CNumbering(0,0,0)],
              have_metis=yes,
              have_metis=no)
fi

if test "x$have_metis" = "xno" ; then
  METIS_LIBS=""
fi

AC_MSG_RESULT(for metis binaries: $have_metis)

AM_CONDITIONAL(HAVE_METIS, test x$have_metis = xyes)

AC_MSG_RESULT(for metis: $have_metis)

CPPFLAGS="$saved_CPPFLAGS"
LDFLAGS="$saved_LDFLAGS"
LIBS="$saved_LIBS"

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(METIS_CPPFLAGS)
AC_SUBST(METIS_LDFLAGS)
AC_SUBST(METIS_LIBS)

AC_LANG_RESTORE

])dnl
