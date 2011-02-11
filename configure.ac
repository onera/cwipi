# Process this file with autoconf to produce a configure script.

#------------------------------------------------------------------------------
# Initialization
#------------------------------------------------------------------------------

m4_define([licence_c_comment],
[/*

*/])

m4_define([cwipi_major_version], [0])
m4_define([cwipi_minor_version], [3])
m4_define([cwipi_release_version], [0])
m4_define([cwipi_version_extra], [])
m4_define([cwipi_version_string],
          [cwipi_major_version.cwipi_minor_version.cwipi_release_version@&t@cwipi_version_extra])

AC_INIT([cwipi],[cwipi_version_string],[eric.quemerais@onera.fr])
AC_CONFIG_SRCDIR([src/cwipi.h])

# Use the config directory for libtool stuff ... 
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_SUBST(ac_aux_dir)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([config_priv.h:config_priv.h.in])

AM_INIT_AUTOMAKE

AH_TOP([#ifndef __CONFIG_PRIV_H__
#define __CONFIG_PRIV_H__

/*
  This file is part of the CWIPI library, intended to interpolate and exchange fields between meshes distributed on many processus 
*/])

AH_BOTTOM([#endif /* __CONFIG_PRIV_H__ */])

#Define public configuration file

CWIPI_AC_CONFIG_PUBL_INIT([cwipi_config.h])
#CWIPI_AC_CONFIG_PUBL_VERBATIM([cwipi_licence_c_comment])

AM_MAINTAINER_MODE

# Version information

CWIPI_AC_CONFIG_PUBL_DEFINE_STRING(CWIPI_VERSION, cwipi_version_string,
                                 [Complete version number])

CWIPI_AC_CONFIG_PUBL_DEFINE(CWIPI_MAJOR_VERSION, cwipi_major_version,
                          [Major version number])
CWIPI_AC_CONFIG_PUBL_DEFINE(CWIPI_MINOR_VERSION, cwipi_minor_version,
                          [Minor version number])
CWIPI_AC_CONFIG_PUBL_DEFINE(CWIPI_RELEASE_VERSION, cwipi_release_version,
                          [Release version number])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_LDFLAGS=$LDFLAGS



AC_PROG_CC(cc gcc icc)

#LT_INIT
#LT_LANG([C++])

#LT_INIT
AC_PROG_CXX(g++ icpc)

#LT_INIT
#LT_LANG([Fortran])

#LT_INIT
AC_PROG_FC(gfortran ifort)

AC_PROG_LD(g++)

AC_PROG_INSTALL

AC_PROG_LN_S

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

# Debug or production compilation mode (debug by default) ?

AC_ARG_ENABLE(debug,
  [  --enable-debug          enable debugging (reduces optimization)],
  [
    case "${enableval}" in
      yes) debug=true ;;
      no)  debug=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
    esac
  ],
  [ debug=false ]
)
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

# Debug or production compilation mode (debug by default) ?

AC_ARG_ENABLE(protect,
  [  --enable-protect          enable use limitation (december 31, 2010)],
  [
    case "${enableval}" in
      yes) protect=yes ;;
      no)  protect=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-protect]) ;;
    esac
  ],
  [ protect=no ]
)

CWIPI_AC_CONFIG_PUBL_SET(PROTECT, $protect,
                       [cwipi library protected ?])

# Optionally deactivate automatic determination of flags on known systems

AC_ARG_ENABLE(auto-flags,
  [  --disable-auto-flags    Do not define *FLAGS on known systems],
  [
    case "${enableval}" in
      yes) auto_flags=true ;;
      no)  auto_flags=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-auto-flags]) ;;
    esac
  ],
  [ auto_flags=true ]
)

# Default compiler options (may be modified
# by defining CFLAGS in the environment)
#------------------------------------------

if test "x$auto_flags" = "xtrue" ; then

  # Source associated recommended compiler options

  if test -f "$srcdir/config/cwipi_auto_flags.sh" ; then
    AC_MSG_NOTICE([sourcing config/cwipi_auto_flags.sh (test for known compilers)])
    . "$srcdir/config/cwipi_auto_flags.sh"
  else
    AC_MSG_WARN([config/cwipi_auto_flags.sh default configuration file not found])
  fi

  if test "x$user_CPPFLAGS" = "x"; then
    CPPFLAGS="$cppflags_default"
  else
    AC_MSG_NOTICE([default CPPFLAGS overridden by user definition])
  fi


  if test "x$user_CXXFLAGS" = "x"; then

    CXXFLAGS="$cxxflags_default"

    if test "x$debug" = xtrue; then
      CXXFLAGS="${CXXFLAGS} $cxxflags_default_dbg"
    else
      CXXFLAGS="${CXXFLAGS} $cxxflags_default_opt"
    fi

  else

    AC_MSG_NOTICE([default CXXFLAGS overridden by user definition])

  fi

  if test "x$user_FCFLAGS" = "x"; then

    FCFLAGS="$fcflags_default"

    if test "x$debug" = xtrue; then
      FCFLAGS="${FCFLAGS} $fcflags_default_dbg"
    else
      FCFLAGS="${FCFLAGS} $fcflags_default_opt"
    fi

  else

    AC_MSG_NOTICE([default FCFLAGS overridden by user definition])

  fi

  if test "x$user_CFLAGS" = "x"; then

    CFLAGS="$cflags_default"

    if test "x$debug" = xtrue; then
      CFLAGS="${CFLAGS} $cflags_default_dbg"
    else
      CFLAGS="${CFLAGS} $cflags_default_opt"
    fi

  else

    AC_MSG_NOTICE([default CFLAGS overridden by user definition])

  fi

  if test "x$user_LDFLAGS" = "x"; then

    LDFLAGS="$ldflags_default"

    if test "x$debug" = xtrue; then
      LDFLAGS="${LDFLAGS} $ldflags_default_dbg"
    else
      LDFLAGS="${LDFLAGS} $ldflags_default_opt"
    fi

  else

    AC_MSG_NOTICE([default LDFLAGS overridden by user definition])

  fi

fi

#------------------------------------------------------------------------------
# Checks for optional Python bindings with SWIG support
#------------------------------------------------------------------------------

# Check for optional Python bindings with SWIG support

AC_ARG_ENABLE(python-bindings,
  [AS_HELP_STRING([--disable-python-bindings],[build without Python interface])],
  [
    case "${enableval}" in
      yes) python_bindings=true ;;
      no)  python_bindings=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-python-bindings]) ;;
    esac
  ],
  [ python_bindings=true ]
)
AM_CONDITIONAL(ENABLE_PYTHON, test x$python_bindings = xtrue)


AC_ARG_WITH(python-exec, [AS_HELP_STRING([--with-python-exec=PATH],
                                         [specify prefix directory for a Python interpreter])],
                         [PYTHON_EXE=${withval}],
                         [PYTHON_EXE=none])dnl

if test "x${python_bindings}" != "xfalse"; then

    #------------------
    # Checks for Python
    #------------------

    if test "x${PYTHON_EXE}" = "xnone" ; then
        AC_PATH_PROG([PYEXE], [python], ["no"])
    else
        AC_PATH_PROG([PYEXE], [python], ["no"], [${PYTHON_EXE}])
    fi

    PYTHON=${PYEXE}
    AM_PATH_PYTHON([2.6],, PYTHON="no")

    if test "x${PYTHON}" = "xno" ; then
        AC_MSG_WARN([])
        AC_MSG_WARN([Python interpreter not found])
        AC_MSG_WARN([Please install Python from http://www.python.org])
        AC_MSG_WARN([])
        AC_MSG_ERROR([Sorry, we can't continue without Python interpreter])
    fi

    PYTHON_INCLUDE=`${PYTHON} -c 'from distutils import sysconfig; print( sysconfig.get_python_inc() )'`

    if test "x${PYTHON_INCLUDE}" = "x" ; then
        AC_MSG_WARN([])
        AC_MSG_WARN([Python include not found])
        AC_MSG_WARN([Please install Python from http://www.python.org])
        AC_MSG_WARN([])
        AC_MSG_ERROR([Sorry, we can't continue without Python include])
    fi

    MPI4PY_INCLUDE=`${PYTHON} -c 'import mpi4py; print( mpi4py.get_include() )'`

    if test "x${MPI4PY_INCLUDE}" = "x" ; then
        AC_MSG_WARN([])
        AC_MSG_WARN([mpi4py package not found])
        AC_MSG_WARN([Please install mpi4py from http://code.google.com/p/mpi4py])
        AC_MSG_WARN([])
        AC_MSG_ERROR([Sorry, we can't continue without mpi4py])
    fi

    NUMPY_INCLUDE=`${PYTHON} -c 'import numpy; print( numpy.get_include() )'`

    if test "x${NUMPY_INCLUDE}" = "x" ; then
        AC_MSG_WARN([])
        AC_MSG_WARN([numpy package not found])
        AC_MSG_WARN([Please install mpi4py from http://numpy.scipy.org])
        AC_MSG_WARN([])
        AC_MSG_ERROR([Sorry, we can't continue without numpy])
    fi

    AC_SUBST(PYTHON_INCLUDE)

    AC_SUBST(MPI4PY_INCLUDE)

    AC_SUBST(NUMPY_INCLUDE)

else
    AM_CONDITIONAL(ENABLE_PYTHON, false)
fi

# Optional additional user flags
#-------------------------------

AC_ARG_VAR(CPP_ADD,  [Additional user C preprocessor flags])
AC_ARG_VAR(CC_ADD,   [Additional user C compiler flags])
AC_ARG_VAR(LD_ADD,   [Additional user linker flags, e.g. -L<libdir>])
AC_ARG_VAR(LIBS_ADD, [Additional user librairies, e.g. -l<lib>])

if test "x$CPP_ADD" != x ; then
  CPPFLAGS="${CPPFLAGS} $CPP_ADD"
fi

if test "x$CC_ADD" != x ; then
  CFLAGS="${CFLAGS} $CC_ADD"
fi

if test "x$LD_ADD" != x ; then
  LDFLAGS="${LDFLAGS} $LD_ADD"
fi

if test "x$LIBS_ADD" != x ; then
  LIBS="${LIBS} $LIBS_ADD"
fi

#------------------------------------------------------------------------------

# We may only turn on processing for libtool now that the basic compiler
# and linker flags are set (to avoid issues with linkers with different
# modes such as 32 and 64 bit which may be modified by compiler or
# linker flags).

# AC_DISABLE_SHARED may not be used inside a test, as it seems to
# take effect whether the result is true or not (maybe due to a bug),
# so we directly use autoconf's enable_shared variable.

if test "x$cwipi_disable_shared" = "xyes" ; then
  enable_shared=no
fi

LT_INIT


#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

# Check for BFT (modifies CPPFLAGS, LDFLAGS, and LIBS)
#-----------------------------------------------------

CHECK_BFT(["1.1.1"])

# Check for FVM (modifies CPPFLAGS, LDFLAGS, and LIBS)
#-----------------------------------------------------

CHECK_FVM(["0.15.1"])

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([sys/types.h stddef.h])

#------------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------------

AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_PID_T
AC_HEADER_TIME
AC_HEADER_STDBOOL

AC_CHECK_TYPES([long long, unsigned long long])
AC_CHECK_TYPES([int32_t])
AC_CHECK_TYPES([int64_t])
AC_CHECK_TYPES([uint32_t])
AC_CHECK_TYPES([uint64_t])
AC_CHECK_TYPES([ptrdiff_t])

#FVM_AC_CHECK_SIZEOF([void *], [fvm])
#FVM_AC_CHECK_SIZEOF([short], [fvm])
#FVM_AC_CHECK_SIZEOF([int], [fvm])
#FVM_AC_CHECK_SIZEOF([long], [fvm])
#FVM_AC_CHECK_SIZEOF([long long], [fvm])
#FVM_AC_CHECK_SIZEOF([float], [fvm])
#FVM_AC_CHECK_SIZEOF([double], [fvm])

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------

AC_CHECK_FUNCS([mkdir])
AC_CHECK_FUNCS([snprintf])

#------------------------------------------------------------------------------
# Checks for system services.
#------------------------------------------------------------------------------

#---------------------
# MPI support (use CC=mpicc with configure if necessary)
#---------------------

MPI_CPPFLAGS=""
MPI_LDFLAGS=""
MPI_LIBS=""

CHECK_MPI

# Whether we used mpicc or $CC, we now know if we have MPI support
CWIPI_AC_CONFIG_PUBL_SET(HAVE_MPI, $have_mpi,
                       [cwipi library configured with MPI support ?])
if test "x$have_mpi" = "xyes" ; then
  CWIPI_AC_CONFIG_PUBL_SET(HAVE_MPI_IO, $have_mpi_io,
                         [cwipi library configured with MPI IO support ?])
  CWIPI_AC_CONFIG_PUBL_SET(HAVE_MPI_ONE_SIDED, $have_mpi_one_sided,
                         [cwipi library configured with MPI One-sided support ?])
fi

#-----------------------------------------------------
# Check for METIS (modifies CPPFLAGS, LDFLAGS, and LIBS)
#-----------------------------------------------------

#
# A supprimer définitivement

#CHECK_METIS


#CWIPI_AC_CONFIG_PUBL_SET(HAVE_METIS, $have_metis,
#                       [cwipi library configured with METIS support ?])

#------------------------------------------------------------------------------
# Internationalization
#------------------------------------------------------------------------------

if test "${prefix}" = "NONE"; then
  localedir="/usr/local/share/locale"
else
  localedir="${prefix}/share/locale"
fi
AC_DEFINE_UNQUOTED(LOCALEDIR, "${localedir}", LOCALEDIR)

# Set of available languages now in po/LINGUAS.

# AM_GNU_GETTEXT

#------------------------------------------------------------------------------
# Compiler options and system characteristics
#------------------------------------------------------------------------------

AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------

# Determine CPPFLAGS, and LDFLAGS required for external code using FVM
#---------------------------------------------------------------------

for flag in [$CPPFLAGS] ; do
  case [$flag] in
    -I* | -isystem*)
    user_cppflags="${user_cppflags} $flag"
    ;;
  esac
done

for flag in [$LDFLAGS] ; do
  case [$flag] in
    -L* | -Wl,-L*)
    user_ldflags="${user_ldflags} $flag"
    ;;
  esac
done

# Output configuration info script
#---------------------------------

CWIPI_AC_CONFIG_INFO_INIT([cwipi-config])

CWIPI_AC_CONFIG_INFO_EXTRA([
# use debugging code:                         "$debug"
# MPI (Message Passing Interface) support:    "$have_mpi"
])

CWIPI_AC_CONFIG_INFO_CC([${CC}], [${cwipi_ac_cc_version}], [${cwipi_ac_cc_version_full}])

CWIPI_AC_CONFIG_INFO_VERSION([cwipi_version_string])

CWIPI_AC_CONFIG_INFO_DIRS([$prefix], [$exec_prefix], [$includedir], [$libdir])

CWIPI_AC_CONFIG_INFO_FLAGS([$user_cppflags], [], [$user_ldflags], [$LIBS])

CWIPI_AC_CONFIG_INFO_FLAGS([$CPPFLAGS], [$CFLAGS], [$LDFLAGS], [$LIBS],
                         [build])

CWIPI_AC_CONFIG_INFO_FLAGS([$BFT_CPPFLAGS], [], [$BFT_LDFLAGS], [$BFT_LIBS],
                         [bft])

CWIPI_AC_CONFIG_INFO_FLAGS([$FVM_CPPFLAGS], [], [$FVM_LDFLAGS], [$FVM_LIBS],
                         [fvm])

CWIPI_AC_CONFIG_INFO_FLAGS([$MPI_CPPFLAGS], [], [$MPI_LDFLAGS], [$MPI_LIBS],
                         [mpi])

CWIPI_AC_CONFIG_INFO_FINALIZE

#------------------------------------------------------------------------------

# If we use GNU gettext, we need intl/Makefile, po/Makefile, and m4/Makefile.

# AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile doc/Makefile
#                  intl/Makefile po/Makefile.in m4/Makefile])

#AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile doc/Makefile])
AC_CONFIG_FILES([Makefile src/Makefile cython/Makefile tests/Makefile])

#------------------------------------------------------------------------------

CWIPI_AC_CONFIG_PUBL_FINALIZE

AC_OUTPUT

#------------------------------------------------------------------------------

# Libtool configuration bug workaround (problem encountered with compilers
# which do not recognize "-soname", so that "-Wl,-soname" is needed to pass
# flags to the linker; Manually setting $wl to "-Wl" in libtool corrects
# the problem, and so does this workaround borrowed from
# the HDF5 configure.in, which does it automatically).

case $host_os in
  linux* | freebsd* )
    # If gcc is not used, need to set $wl to use "-Wl,"
    if $CC -v 2>&1 | grep '^gcc' > /dev/null ; then
      : using gcc
    else
      echo 'fixing $wl in' $ofile
ed - $ofile <<EOF 2> /dev/null
/^wl=""/s//wl="-Wl,"/
w
q
EOF
    fi
    ;;
esac

#case "$INSTALL" in
#  *install-sh*)
#    INSTALL='\${top_srcdir}/bin/install-sh -c'
#    ;;
#esac

#------------------------------------------------------------------------------

echo
echo "Configuration options:"
echo " use debugging code: "$debug""
echo " MPI (Message Passing Interface) support: "$have_mpi""
if (test x$have_mpi = xyes) ; then
  echo "   MPI I/O support: "$fvm_have_mpi_io""
  echo "   MPI2 one-sided communication support: "$have_mpi_one_sided""
fi
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the CWIPI package"
echo

