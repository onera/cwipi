
# CHECK_FVM(Minimal Release string, [Maximal Release string])
#------------------------------------------------------------------
# Check for FVM version ; defines FVM_CPPFLAGS, FVM_LDFLAGS, and FVM_LIBS
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

if test "x$with_fvm_exec" != "x" ; then
  fvm_config="$with_fvm_exec/fvm-config"
elif test "x$with_fvm" != "x" ; then
  fvm_config="$with_fvm/bin/fvm-config"
else
  fvm_config="fvm-config"
fi

if test "x$with_fvm_include" != "x" ; then
  FVM_CPPFLAGS="-I$with_fvm_include"
elif test "x$with_fvm" != "x" ; then
  if test -f "$with_fvm/include/fvm_config.h" ; then
    FVM_CPPFLAGS="-I$with_fvm/include"
  elif test -f "$with_fvm/include/fvm/fvm_config.h" ; then
    FVM_CPPFLAGS="-I$with_fvm/include/fvm"
  fi
else
  fvm_prefix=`fvm_config --prefix`
  if test ! -f "$fvm_prefix/include/fvm_config.h" ; then
    if test -f "$fvm_prefix/include/fvm/fvm_config.h" ; then
      FVM_CPPFLAGS="-I$fvm_prefix/include/fvm"
    fi
  fi
fi

if test "x$with_fvm_lib" != "x" ; then
  FVM_LDFLAGS="-L$with_fvm_lib"
elif test "x$with_fvm" != "x" ; then
  FVM_LDFLAGS="-L$with_fvm/lib"
else
  FVM_LDFLAGS=""
fi
FVM_LIBS="-lfvm"

type "$fvm_config" > /dev/null 2>&1
if test "$?" = "0" ; then
  FVM_CPPFLAGS="$FVM_CPPFLAGS `$fvm_config --cppflags`"
  FVM_LDFLAGS="$FVM_LDFLAGS `$fvm_config --ldflags`"
  FVM_LIBS="$FVM_LIBS `$fvm_config --libs`"
fi

fvm_version_min=$1
fvm_version_max=$2

if test "x$fvm_version_min" != "x" ; then
  if test "x$fvm_version_max" != "x" ; then
    AC_MSG_CHECKING([for fvm version >= $1 and <= $2])
  else
    AC_MSG_CHECKING([for fvm version >= $1])
  fi
else
  fvm_version_min="0.0.0"
fi

fvm_version_major_min=`echo "$fvm_version_min" | cut -f1 -d.`
fvm_version_minor_min=`echo "$fvm_version_min" | cut -f2 -d.`
fvm_version_release_min=`echo "$fvm_version_min" | cut -f3 -d.`
if test    "$fvm_version_major_min" = "" \
        -o "$fvm_version_minor_min" = "" \
        -o "$fvm_version_release_min" = ""; then
  AC_MSG_FAILURE([bad FVM version definition in configure.ac: $fvm_version_min])
fi

if test "x$fvm_version_max" != "x" ; then
  fvm_version_major_max=`echo "$fvm_version_max" | cut -f1 -d.`
  fvm_version_minor_max=`echo "$fvm_version_max" | cut -f2 -d.`
  fvm_version_release_max=`echo "$fvm_version_max" | cut -f3 -d.`
  if test    "$fvm_version_major_max" = "" \
          -o "$fvm_version_minor_max" = "" \
          -o "$fvm_version_release_max" = ""; then
    AC_MSG_FAILURE([bad FVM version definition in configure.ac: $fvm_version_max])
  fi
else
  fvm_version_major_max=99999999
  fvm_version_minor_max=99999999
  fvm_version_release_max=99999999
fi

saved_CPPFLAGS=$CPPFLAGS
saved_LDFLAGS=$LDFLAGS
saved_LIBS=$LIBS

CPPFLAGS="${CPPFLAGS} ${BFT_CPPFLAGS} $FVM_CPPFLAGS"
LDFLAGS="${LDFLAGS} ${BFT_LDFLAGS} $FVM_LDFLAGS"
LIBS="${LIBS} ${BFT_LIBS} $FVM_LIBS"

AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fvm_config.h>
]],
[[#if FVM_MAJOR_VERSION < $fvm_version_major_min
#  error FVM major version < $fvm_version_major_min
#elif FVM_MAJOR_VERSION == $fvm_version_major_min
#  if FVM_MINOR_VERSION < $fvm_version_minor_min
#    error FVM minor version < $fvm_version_minor_min
#  elif FVM_MINOR_VERSION == $fvm_version_minor_min
#    if FVM_RELEASE_VERSION < $fvm_version_release_min
#      error FVM release version < $fvm_version_release_min
#    endif
#  endif
#endif
#if FVM_MAJOR_VERSION > $fvm_version_major_max
#  error FVM major version > $fvm_version_major_max
#elif FVM_MAJOR_VERSION == $fvm_version_major_max
#  if FVM_MINOR_VERSION > $fvm_version_minor_max
#    error FVM minor version > $fvm_version_minor_max
#  elif FVM_MINOR_VERSION == $fvm_version_minor_max
#    if FVM_RELEASE_VERSION > $fvm_version_release_max
#      error FVM release version < $fvm_version_release_max
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible fvm version found])],
               [AC_MSG_FAILURE([compatible fvm version not found])])

unset fvm_version_major_min
unset fvm_version_minor_min
unset fvm_version_release_min
unset fvm_version_major_max
unset fvm_version_minor_max
unset fvm_version_release_max

unset fvm_version_min
unset fvm_version_max
unset fvm_config

# Restore old LIBS to add $FVM_LIBS later, as other tests
# might otherwise not run if a shared library is not found

CPPFLAGS=$saved_CPPFLAGS
LDFLAGS=$saved_LDFLAGS
LIBS=$saved_LIBS

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(FVM_CPPFLAGS)
AC_SUBST(FVM_LDFLAGS)
AC_SUBST(FVM_LIBS)

])dnl
