
# CHECK_BFT(Minimal Release string, [Maximal Release string])
#------------------------------------------------------------------
# Check for BFT version ; defines BFTC_CPPFLAGS, BFTC_LDFLAGS, and BFTC_LIBS
# locally (i.e. as simple variables, not AC_SUBST)

AC_DEFUN([CHECK_BFT], [

AC_ARG_WITH(bft,
  [  --with-bft=PATH         specify prefix directory for BFT]
)

AC_ARG_WITH(bft-exec,
  [  --with-bft-exec=PATH    specify directory for BFT executables]
)

AC_ARG_WITH(bft-include,
  [  --with-bft-include=PATH  specify directory for BFT include files]
)

AC_ARG_WITH(bft-lib,
  [  --with-bft-lib=PATH     specify directory for BFT library]
)

if test "x$with_bftc_exec" != "x" ; then
  bftc_config="$with_bftc_exec/bft-config"
elif test "x$with_bft" != "x" ; then
  bftc_config="$with_bft/bin/bft-config"
else
  bftc_config="bft-config"
fi

if test "x$with_bftc_include" != "x" ; then
  BFTC_CPPFLAGS="-I$with_bftc_include"
elif test "x$with_bft" != "x" ; then
  if test -f "$with_bft/include/bftc_config.h" ; then
    BFTC_CPPFLAGS="-I$with_bft/include"
  elif test -f "$with_bft/include/bft/bftc_config.h" ; then
    BFTC_CPPFLAGS="-I$with_bft/include/bft"
  fi
else
  bftc_prefix=`bftc_config --prefix`
  if test ! -f "$bftc_prefix/include/bftc_config.h" ; then
    if test -f "$bftc_prefix/include/bft/bftc_config.h" ; then
      BFTC_CPPFLAGS="-I$bftc_prefix/include/bft"
    fi
  fi
fi

if test "x$with_bftc_lib" != "x" ; then
  BFTC_LDFLAGS="-L$with_bftc_lib"
elif test "x$with_bft" != "x" ; then
  BFTC_LDFLAGS="-L$with_bft/lib"
else
  BFTC_LDFLAGS=""
fi
BFTC_LIBS="-lbft"

type "$bftc_config" > /dev/null 2>&1
if test "$?" = "0" ; then
  BFTC_CPPFLAGS="$BFTC_CPPFLAGS `$bftc_config --cppflags`"
  BFTC_LDFLAGS="$BFTC_LDFLAGS `$bftc_config --ldflags`"
  BFTC_LIBS="$BFTC_LIBS `$bftc_config --libs`"
fi

bftc_version_min=$1
bftc_version_max=$2

if test "x$bftc_version_min" != "x" ; then
  if test "x$bftc_version_max" != "x" ; then
    AC_MSG_CHECKING([for bft version >= $1 and <= $2])
  else
    AC_MSG_CHECKING([for bft version >= $1])
  fi
else
  bftc_version_min="0.0.0"
fi

bftc_version_major_min=`echo "$bftc_version_min" | cut -f1 -d.`
bftc_version_minor_min=`echo "$bftc_version_min" | cut -f2 -d.`
bftc_version_release_min=`echo "$bftc_version_min" | cut -f3 -d.`
if test    "$bftc_version_major_min" = "" \
        -o "$bftc_version_minor_min" = "" \
        -o "$bftc_version_release_min" = ""; then
  AC_MSG_FAILURE([bad BFT version definition in configure.ac: $bftc_version_min])
fi

if test "x$bftc_version_max" != "x" ; then
  bftc_version_major_max=`echo "$bftc_version_max" | cut -f1 -d.`
  bftc_version_minor_max=`echo "$bftc_version_max" | cut -f2 -d.`
  bftc_version_release_max=`echo "$bftc_version_max" | cut -f3 -d.`
  if test    "$bftc_version_major_max" = "" \
          -o "$bftc_version_minor_max" = "" \
          -o "$bftc_version_release_max" = ""; then
    AC_MSG_FAILURE([bad BFT version definition in configure.ac: $bftc_version_max])
  fi
else
  bftc_version_major_max=99999999
  bftc_version_minor_max=99999999
  bftc_version_release_max=99999999
fi

saved_CPPFLAGS=$CPPFLAGS
saved_LDFLAGS=$LDFLAGS
saved_LIBS=$LIBS

CPPFLAGS="${CPPFLAGS} $BFTC_CPPFLAGS"
LDFLAGS="${LDFLAGS} $BFTC_LDFLAGS"
LIBS="${LIBS} $BFTC_LIBS"

AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <bftc_config.h>
]],
[[#if BFTC_MAJOR_VERSION < $bftc_version_major_min
#  error BFT major version < $bftc_version_major_min
#elif BFTC_MAJOR_VERSION == $bftc_version_major_min
#  if BFTC_MINOR_VERSION < $bftc_version_minor_min
#    error BFT minor version < $bftc_version_minor_min
#  elif BFTC_MINOR_VERSION == $bftc_version_minor_min
#    if BFTC_RELEASE_VERSION < $bftc_version_release_min
#      error BFT release version < $bftc_version_release_min
#    endif
#  endif
#endif
#if BFTC_MAJOR_VERSION > $bftc_version_major_max
#  error BFT major version > $bftc_version_major_max
#elif BFTC_MAJOR_VERSION == $bftc_version_major_max
#  if BFTC_MINOR_VERSION > $bftc_version_minor_max
#    error BFT minor version > $bftc_version_minor_max
#  elif BFTC_MINOR_VERSION == $bftc_version_minor_max
#    if BFTC_RELEASE_VERSION > $bftc_version_release_max
#      error BFT release version < $bftc_version_release_max
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible bft version found])],
               [AC_MSG_FAILURE([compatible bft version not found])])

unset bftc_version_major_min
unset bftc_version_minor_min
unset bftc_version_release_min
unset bftc_version_major_max
unset bftc_version_minor_max
unset bftc_version_release_max

unset bftc_version_min
unset bftc_version_max
unset bftc_config

# Restore old LIBS to add $BFTC_LIBS later, as other tests
# might otherwise not run if a shared library is not found

CPPFLAGS=$saved_CPPFLAGS
LDFLAGS=$saved_LDFLAGS
LIBS=$saved_LIBS

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(BFTC_CPPFLAGS)
AC_SUBST(BFTC_LDFLAGS)
AC_SUBST(BFTC_LIBS)

])dnl
