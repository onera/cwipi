
# CHECK_FVM(Minimal Release string, [Maximal Release string])
#------------------------------------------------------------------
# Check for FVM version ; defines FVMC_CPPFLAGS, FVMC_LDFLAGS, and FVMC_LIBS
# locally (i.e. as simple variables, not AC_SUBST)

AC_DEFUN([CHECK_FVM], [

AC_ARG_WITH(fvm,
  [  --with-fvm=PATH         specify prefix directory for FVM]
)

AC_ARG_WITH(fvm-exec,
  [  --with-fvm-exec=PATH    specify directory for FVM executables]
)

AC_ARG_WITH(fvm-include,
  [  --with-fvm-include=PATH  specify directory for FVM include files]
)

AC_ARG_WITH(fvm-lib,
  [  --with-fvm-lib=PATH     specify directory for FVM library]
)

if test "x$with_fvmc_exec" != "x" ; then
  fvmc_config="$with_fvmc_exec/fvm-config"
elif test "x$with_fvm" != "x" ; then
  fvmc_config="$with_fvm/bin/fvm-config"
else
  fvmc_config="fvm-config"
fi

if test "x$with_fvmc_include" != "x" ; then
  FVMC_CPPFLAGS="-I$with_fvmc_include"
elif test "x$with_fvm" != "x" ; then
  if test -f "$with_fvm/include/fvmc_config.h" ; then
    FVMC_CPPFLAGS="-I$with_fvm/include"
  elif test -f "$with_fvm/include/fvm/fvmc_config.h" ; then
    FVMC_CPPFLAGS="-I$with_fvm/include/fvm"
  fi
else
  fvmc_prefix=`fvmc_config --prefix`
  if test ! -f "$fvmc_prefix/include/fvmc_config.h" ; then
    if test -f "$fvmc_prefix/include/fvm/fvmc_config.h" ; then
      FVMC_CPPFLAGS="-I$fvmc_prefix/include/fvm"
    fi
  fi
fi

if test "x$with_fvmc_lib" != "x" ; then
  FVMC_LDFLAGS="-L$with_fvmc_lib"
elif test "x$with_fvm" != "x" ; then
  FVMC_LDFLAGS="-L$with_fvm/lib"
else
  FVMC_LDFLAGS=""
fi
FVMC_LIBS="-lfvm"

type "$fvmc_config" > /dev/null 2>&1
if test "$?" = "0" ; then
  FVMC_CPPFLAGS="$FVMC_CPPFLAGS `$fvmc_config --cppflags`"
  FVMC_LDFLAGS="$FVMC_LDFLAGS `$fvmc_config --ldflags`"
  FVMC_LIBS="$FVMC_LIBS `$fvmc_config --libs`"
fi

fvmc_version_min=$1
fvmc_version_max=$2

if test "x$fvmc_version_min" != "x" ; then
  if test "x$fvmc_version_max" != "x" ; then
    AC_MSG_CHECKING([for fvm version >= $1 and <= $2])
  else
    AC_MSG_CHECKING([for fvm version >= $1])
  fi
else
  fvmc_version_min="0.0.0"
fi

fvmc_version_major_min=`echo "$fvmc_version_min" | cut -f1 -d.`
fvmc_version_minor_min=`echo "$fvmc_version_min" | cut -f2 -d.`
fvmc_version_release_min=`echo "$fvmc_version_min" | cut -f3 -d.`
if test    "$fvmc_version_major_min" = "" \
        -o "$fvmc_version_minor_min" = "" \
        -o "$fvmc_version_release_min" = ""; then
  AC_MSG_FAILURE([bad FVM version definition in configure.ac: $fvmc_version_min])
fi

if test "x$fvmc_version_max" != "x" ; then
  fvmc_version_major_max=`echo "$fvmc_version_max" | cut -f1 -d.`
  fvmc_version_minor_max=`echo "$fvmc_version_max" | cut -f2 -d.`
  fvmc_version_release_max=`echo "$fvmc_version_max" | cut -f3 -d.`
  if test    "$fvmc_version_major_max" = "" \
          -o "$fvmc_version_minor_max" = "" \
          -o "$fvmc_version_release_max" = ""; then
    AC_MSG_FAILURE([bad FVM version definition in configure.ac: $fvmc_version_max])
  fi
else
  fvmc_version_major_max=99999999
  fvmc_version_minor_max=99999999
  fvmc_version_release_max=99999999
fi

saved_CPPFLAGS=$CPPFLAGS
saved_LDFLAGS=$LDFLAGS
saved_LIBS=$LIBS

CPPFLAGS="${CPPFLAGS} ${BFTC_CPPFLAGS} $FVMC_CPPFLAGS"
LDFLAGS="${LDFLAGS} ${BFTC_LDFLAGS} $FVMC_LDFLAGS"
LIBS="${LIBS} ${BFTC_LIBS} $FVMC_LIBS"

AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fvmc_config.h>
]],
[[#if FVMC_MAJOR_VERSION < $fvmc_version_major_min
#  error FVM major version < $fvmc_version_major_min
#elif FVMC_MAJOR_VERSION == $fvmc_version_major_min
#  if FVMC_MINOR_VERSION < $fvmc_version_minor_min
#    error FVM minor version < $fvmc_version_minor_min
#  elif FVMC_MINOR_VERSION == $fvmc_version_minor_min
#    if FVMC_RELEASE_VERSION < $fvmc_version_release_min
#      error FVM release version < $fvmc_version_release_min
#    endif
#  endif
#endif
#if FVMC_MAJOR_VERSION > $fvmc_version_major_max
#  error FVM major version > $fvmc_version_major_max
#elif FVMC_MAJOR_VERSION == $fvmc_version_major_max
#  if FVMC_MINOR_VERSION > $fvmc_version_minor_max
#    error FVM minor version > $fvmc_version_minor_max
#  elif FVMC_MINOR_VERSION == $fvmc_version_minor_max
#    if FVMC_RELEASE_VERSION > $fvmc_version_release_max
#      error FVM release version < $fvmc_version_release_max
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible fvm version found])],
               [AC_MSG_FAILURE([compatible fvm version not found])])

unset fvmc_version_major_min
unset fvmc_version_minor_min
unset fvmc_version_release_min
unset fvmc_version_major_max
unset fvmc_version_minor_max
unset fvmc_version_release_max

unset fvmc_version_min
unset fvmc_version_max
unset fvmc_config

# Restore old LIBS to add $FVMC_LIBS later, as other tests
# might otherwise not run if a shared library is not found

CPPFLAGS=$saved_CPPFLAGS
LDFLAGS=$saved_LDFLAGS
LIBS=$saved_LIBS

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(FVMC_CPPFLAGS)
AC_SUBST(FVMC_LDFLAGS)
AC_SUBST(FVMC_LIBS)

])dnl
