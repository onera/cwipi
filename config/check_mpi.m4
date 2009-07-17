
# CHECK_MPI
#----------------
# optional MPI support (use CC=mpicc with configure if necessary)
# modifies or sets have_mpi, MPI_CPPFLAGS, MPI_LDFLAGS, and MPI_LIBS
# depending on libraries found

AC_DEFUN([CHECK_MPI], [

saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
saved_LIBS="$LIBS"

have_mpi=no
have_mpi_io=no
have_mpi_one_sided=no

AC_ARG_ENABLE(mpi,
  [  --disable-mpi           do not use MPI when available],
  [
    case "${enableval}" in
      yes) mpi=true ;;
      no)  mpi=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-mpi]) ;;
    esac
  ],
  [ mpi=true ]
)

AC_ARG_ENABLE(mpi-io,
  [  --disable-mpi-io        do not use MPI I/O when available],
  [
    case "${enableval}" in
      yes) mpi_io=true ;;
      no)  mpi_io=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-mpi-io]) ;;
    esac
  ],
  [ mpi_io=true ]
)

AC_ARG_WITH(mpi,
  [  --with-mpi=PATH         specify prefix directory for MPI]
)

AC_ARG_WITH(mpi-include,
  [  --with-mpi-include=PATH specify directory for MPI include files]
)

AC_ARG_WITH(mpi-lib,
  [  --with-mpi-lib=PATH     specify directory for MPI library]
)

if test "x$mpi" = "xtrue" ; then
  if test "x$with_mpi_include" != "x" ; then
    MPI_CPPFLAGS="$MPI_CPPFLAGS -I$with_mpi_include"
  elif test "x$with_mpi" != "x" ; then
    MPI_CPPFLAGS="$MPI_CPPFLAGS -I$with_mpi/include"
  fi
  if test "x$with_mpi_lib" != "x" ; then
    MPI_LDFLAGS="$MPI_LDFLAGS -L$with_mpi_lib"
  elif test "x$with_mpi" != "x" ; then
    MPI_LDFLAGS="$MPI_LDFLAGS -L$with_mpi/lib"
  fi
fi

# Just in case, remove excess whitespace from existing flag and libs variables.

if test "$MPI_CPPFLAGS" != "" ; then
  MPI_CPPFLAGS=`echo $MPI_CPPFLAGS | sed 's/^[ ]*//;s/[ ]*$//'`
fi
if test "$MPI_LDFLAGS" != "" ; then
  MPI_LDFLAGS=`echo $MPI_LDFLAGS | sed 's/^[ ]*//;s/[ ]*$//'`
fi
if test "$MPI_LIBS" != "" ; then
  MPI_LIBS=`echo $MPI_LIBS | sed 's/^[ ]*//;s/[ ]*$//'`
fi

# If we use mpicc or a variant, we have nothing to add

if test "x$mpi" = "xtrue" -a "x$CC" != "x" ; then
  CCNAME=`basename "$CC"`
  # Test for standard wrappers
  if test "$CCNAME" = "mpicc" -o "$CCNAME" = "mpiCC" \
                              -o "$CCNAME" = "mpic++" ; then
    have_mpi=yes
  # Else test for other known wrappers
  elif test "$CCNAME" = "mpcc" -o "$CCNAME" = "mpCC" \
                               -o "$CCNAME" = "mpixlc" ; then
    have_mpi=yes
  fi
  if test "x$have_mpi" = "xyes"; then
    if test "x$mpi_io" = "xtrue"; then
      AC_MSG_CHECKING([for MPI I/O])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_File_close((void *)0); ]])],
                     [have_mpi_io=yes],
                     [have_mpi_io=no])
      AC_MSG_RESULT($have_mpi_io)
    fi
    AC_MSG_CHECKING([for MPI2 one-sided communication])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Win_free((void *)0); ]])],
                   [have_mpi_one_sided=yes],
                   [have_mpi_one_sided=no])
    AC_MSG_RESULT($have_mpi_one_sided)
  fi
fi

# If we do not use an MPI compiler wrapper, we must add compilation
# and link flags; we try to detect the correct flags to add.

if test "x$mpi" = "xtrue" -a "x$have_mpi" = "xno" ; then

  # try several tests for MPI

  # Basic test
  AC_MSG_CHECKING([for MPI (basic test)])
  if test "$MPI_LIBS" = "" ; then
    MPI_LIBS="-lmpi $PTHREAD_LIBS"
  fi
  CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
  LIBS="$saved_LIBS $MPI_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                 [[ MPI_Init(0, (void *)0); ]])],
                 [have_mpi=yes],
                 [have_mpi=no])
  AC_MSG_RESULT($have_mpi)

  # If failed, test for mpich
  if test "x$have_mpi" = "xno"; then
    AC_MSG_CHECKING([for MPI (mpich test)])
    # First try (simplest)
    MPI_LIBS="-lmpich $PTHREAD_LIBS"
    LIBS="$saved_LIBS $MPI_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Init(0, (void *)0); ]])],
                   [have_mpi=yes],
                   [have_mpi=no])
    if test "x$have_mpi" = "xno"; then
      # Second try (with lpmpich)
      MPI_LIBS="-Wl,-lpmpich -Wl,-lmpich -Wl,-lpmpich -Wl,-lmpich"
      LIBS="$saved_LIBS $MPI_LIBS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_Init(0, (void *)0); ]])],
                     [have_mpi=yes],
                     [have_mpi=no])
    fi
    AC_MSG_RESULT($have_mpi)
  fi

  # If failed, test for lam-mpi
  if test "x$have_mpi" = "xno"; then
    AC_MSG_CHECKING([for MPI (lam-mpi test)])
    # First try (without MPI-IO)
    case $host_os in
      freebsd*)
        MPI_LIBS="-lmpi -llam $PTHREAD_LIBS";;
      *)
        MPI_LIBS="-lmpi -llam -lpthread";;
    esac
    LIBS="$saved_LIBS $MPI_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Init(0, (void *)0); ]])],
                   [have_mpi=yes],
                   [have_mpi=no])
    if test "x$have_mpi" = "xno"; then
      # Second try (with MPI-IO)
      case $host_os in
        freebsd*)
          MPI_LIBS="-lmpi -llam -lutil -ldl $PTHREAD_LIBS";;
        *)
          MPI_LIBS="-lmpi -llam -lutil -ldl -lpthread";;
      esac
      LIBS="$saved_LIBS $MPI_LIBS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_Init(0, (void *)0); ]])],
                     [have_mpi=yes],
                     [have_mpi=no])
    fi
    AC_MSG_RESULT($have_mpi)
  fi

  if test "x$have_mpi" = "xyes"; then
    if test "x$mpi_io" = "xtrue"; then
      AC_MSG_CHECKING([for MPI I/O])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_File_close((void *)0); ]])],
                     [have_mpi_io=yes],
                     [have_mpi_io=no])
      AC_MSG_RESULT($have_mpi_io)
    fi
    AC_MSG_CHECKING([for MPI2 one-sided communication])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Win_free((void *)0); ]])],
                   [have_mpi_one_sided=yes],
                   [have_mpi_one_sided=no])
    AC_MSG_RESULT($have_mpi_one_sided)
  else
    MPI_CPPFLAGS=""
    MPI_LDFLAGS=""
    MPI_LIBS=""
  fi

  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

fi

if test "x$have_mpi" = "xno"; then
 AC_MSG_FAILURE([No MPI library])
fi

AM_CONDITIONAL(HAVE_MPI, test x$have_mpi = xyes)
AM_CONDITIONAL(HAVE_MPI_IO, test x$have_mpi_io = xyes)
AM_CONDITIONAL(HAVE_MPI_ONE_SIDED, test x$have_mpi_one_sided = xyes)

CPPFLAGS="$saved_CPPFLAGS"
LDFLAGS="$saved_LDFLAGS"
LIBS="$saved_LIBS"

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(MPI_CPPFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)

])dnl

