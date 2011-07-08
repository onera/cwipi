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

# Process this file with autoconf to produce a configure script.

#------------------------------------------------------------------------------
# Initialization
#------------------------------------------------------------------------------

m4_define([cwipi_licence_c_comment],
[/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/])

m4_define([bft_licence_c_comment],
[/*
  This file is part of the CWIPI library, derived from "Base Functions and Types"
  library, intended to simplify and enhance portability, memory 
  and I/O use for scientific codes.

  Copyright (C) 2004-2010  EDF   

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/])

m4_define([fvm_licence_c_comment],
[/*
  This file is part of the CWIPI library, derived from the "Finite Volume Mesh" 
  library, intended to provide finite volume mesh and associated 
  fields I/O and manipulation services.

  Copyright (C) 2011  ONERA

  Copyright (C) 2004-2010  EDF   

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/])

m4_define([cwipi_major_version], [0])
m4_define([cwipi_minor_version], [4])
m4_define([cwipi_release_version], [5])
m4_define([cwipi_version_extra], [])
m4_define([cwipi_version_string],
          [cwipi_major_version.cwipi_minor_version.cwipi_release_version@&t@cwipi_version_extra])

m4_define([fvm_major_version], [0])
m4_define([fvm_minor_version], [15])
m4_define([fvm_release_version], [1])
m4_define([fvm_version_extra], [-onera-6])
m4_define([fvm_version_string],
          [fvm_major_version.fvm_minor_version.fvm_release_version@&t@fvm_version_extra])

m4_define([bft_major_version], [1])
m4_define([bft_minor_version], [1])
m4_define([bft_release_version], [3])
m4_define([bft_version_string],
          [bft_major_version.bft_minor_version.bft_release_version])

AC_INIT([cwipi],[cwipi_version_string],[eric.quemerais@onera.fr])
AC_CONFIG_SRCDIR([src/cwipi.h])

# Use the config directory for libtool stuff ... 
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_SUBST(ac_aux_dir)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([config_priv.h:config_priv.h.in])
AC_CONFIG_HEADERS([fvm_config_priv.h:fvm_config_priv.h.in])
AC_CONFIG_HEADERS([bft_config_priv.h:bft_config_priv.h.in])

AM_INIT_AUTOMAKE

AH_TOP([#ifndef __CONFIG_PRIV_H__
#define __CONFIG_PRIV_H__

/*
  This file is part of the CWIPI library, intended to interpolate and exchange fields between meshes distributed on many processus 
*/])

AH_BOTTOM([#endif /* __CONFIG_PRIV_H__ */])

#Define public configuration file

CWIPI_AC_CONFIG_PUBL_INIT([cwipi_config.h])
CWIPI_AC_CONFIG_PUBL_VERBATIM([cwipi_licence_c_comment])

FVM_AC_CONFIG_PUBL_INIT([fvm_config.h])
FVM_AC_CONFIG_PUBL_VERBATIM([fvm_licence_c_comment])

BFT_AC_CONFIG_PUBL_INIT([bft_config.h])
BFT_AC_CONFIG_PUBL_VERBATIM([bft_licence_c_comment])


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

FVM_AC_CONFIG_PUBL_DEFINE_STRING(FVM_VERSION, fvm_version_string,
                                 [Complete version number])

FVM_AC_CONFIG_PUBL_DEFINE(FVM_MAJOR_VERSION, fvm_major_version,
                          [Major version number])
FVM_AC_CONFIG_PUBL_DEFINE(FVM_MINOR_VERSION, fvm_minor_version,
                          [Minor version number])
FVM_AC_CONFIG_PUBL_DEFINE(FVM_RELEASE_VERSION, fvm_release_version,
                          [Release version number])

BFT_AC_CONFIG_PUBL_DEFINE_STRING(BFT_VERSION, bft_version_string,
                                 [Complete version number])

BFT_AC_CONFIG_PUBL_DEFINE(BFT_MAJOR_VERSION, bft_major_version,
                          [Major version number])
BFT_AC_CONFIG_PUBL_DEFINE(BFT_MINOR_VERSION, bft_minor_version,
                          [Minor version number])
BFT_AC_CONFIG_PUBL_DEFINE(BFT_RELEASE_VERSION, bft_release_version,
                          [Release version number])

BFT_AC_CONFIG_PUBL_DEFINE(BFT_CPPCALLER, 1,
                          [Active namespace bft for cwipi])
FVM_AC_CONFIG_PUBL_DEFINE(FVM_CPPCALLER, 1,
                          [Active namespace fvm for cwipi])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_LDFLAGS=$LDFLAGS

AC_PROG_CC(cc gcc icc)

AC_PROG_CXX(g++ icpc)

AC_PROG_FC(gfortran ifort)

AC_PROG_LD(g++ icpc)

AC_PROG_INSTALL

AC_PROG_LN_S

#AC_LIBTOOL_DLOPEN
#AC_PROG_LIBTOOL

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
        AC_MSG_WARN([Sorry, we can't continue without Python interpreter])
        AM_CONDITIONAL(ENABLE_PYTHON, false)

    else
    
        PYTHON_INCLUDE=`${PYTHON} -c 'from distutils import sysconfig; print( sysconfig.get_python_inc() )'`

        if test "x${PYTHON_INCLUDE}" = "x" ; then
            AC_MSG_WARN([])
            AC_MSG_WARN([Python include not found])
            AC_MSG_WARN([Please install Python from http://www.python.org])
            AC_MSG_WARN([])
            AC_MSG_WARN([Sorry, we can't continue without Python include])
            PYTHON="no"        
        fi

        MPI4PY_INCLUDE=`${PYTHON} -c 'import mpi4py; print( mpi4py.get_include() )'`

        if test "x${MPI4PY_INCLUDE}" = "x" ; then
            AC_MSG_WARN([])
            AC_MSG_WARN([mpi4py package not found])
            AC_MSG_WARN([Please install mpi4py from http://code.google.com/p/mpi4py])
            AC_MSG_WARN([])
            AC_MSG_WARN([Sorry, we can't continue without mpi4py])
            PYTHON="no"        
        fi

        NUMPY_INCLUDE=`${PYTHON} -c 'import numpy; print( numpy.get_include() )'`

        if test "x${NUMPY_INCLUDE}" = "x" ; then
            AC_MSG_WARN([])
            AC_MSG_WARN([numpy package not found])
            AC_MSG_WARN([Please install mpi4py from http://numpy.scipy.org])
            AC_MSG_WARN([])
            AC_MSG_WARN([Sorry, we can't continue without numpy])
            PYTHON="no"        
        fi

        if test "x${PYTHON}" = "xno" ; then
            AM_CONDITIONAL(ENABLE_PYTHON, false)
        fi
        AC_SUBST(PYTHON_INCLUDE)
        AC_SUBST(MPI4PY_INCLUDE)
        AC_SUBST(NUMPY_INCLUDE)
    fi
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

# Add some preprocessor flags to include additional system features
#------------------------------------------------------------------

case "$host_os" in
  linux*)
    CPPFLAGS="${CPPFLAGS} -D_GNU_SOURCE"
    ;;
  *bsd* | solaris2.* | aix*)
    ;;
  darwin*)
    CPPFLAGS="${CPPFLAGS} -D_DARWIN_C_SOURCE"
    ;;
  osf*)
    CPPFLAGS="${CPPFLAGS} -D_OSF_SOURCE"
    ;;
  hpux*)
    CPPFLAGS="${CPPFLAGS} -D_HPUX_SOURCE"
    ;;
  *)
    ;;
esac

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

#--------------------------------------------------------------------------
# If the system is of a known type, declare it in the configuration header.
#--------------------------------------------------------------------------

case "$host_os" in
  linux*)
    AC_DEFINE(bft_OS_Linux, 1, bft_OS_Linux)
    ;;
  solaris2.*)
    AC_DEFINE(bft_OS_Solaris, 1, bft_OS_SunOS)
    ;;
  aix*)
    AC_DEFINE(bft_OS_AIX, 1, bft_OS_AIX)
    ;;
  osf*)
    AC_DEFINE(bft_OS_OSF1, 1, bft_OS_OSF1)
    ;;
  irix5.*|irix6.*)
    AC_DEFINE(bft_OS_IRIX64, 1, bft_OS_IRIX64)
    ;;
  hpux*)
    AC_DEFINE(bft_OS_HP_UX, 1, bft_OS_HP_UX)
    ;;
  uxpv*)
    AC_DEFINE(bft_OS_UXPV, 1, bft_OS_UNIX_System_V)
    ;;
esac

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

# Check for optionnal HDF support

FVM_AC_TEST_HDF5

# Check for optionnal CGNS support

FVM_AC_TEST_CGNS

# Check for optionnal MED support

FVM_AC_TEST_MED

# Optional Zlib Support

BFT_AC_ZLIB

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h unistd.h fcntl.h sys/types.h sys/signal.h])
AC_CHECK_HEADERS([sys/stat.h sys/procfs.h sys/sysinfo.h sys/utsname.h])
AC_CHECK_HEADERS([sys/resource.h stddef.h])

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

FVM_AC_CHECK_SIZEOF([void *], [fvm])
FVM_AC_CHECK_SIZEOF([short], [fvm])
FVM_AC_CHECK_SIZEOF([int], [fvm])
FVM_AC_CHECK_SIZEOF([long], [fvm])
FVM_AC_CHECK_SIZEOF([long long], [fvm])
FVM_AC_CHECK_SIZEOF([float], [fvm])
FVM_AC_CHECK_SIZEOF([double], [fvm])

BFT_AC_CHECK_SIZEOF([void *],    [bft])
BFT_AC_CHECK_SIZEOF([short],     [bft])
BFT_AC_CHECK_SIZEOF([int],       [bft])
BFT_AC_CHECK_SIZEOF([long],      [bft])
BFT_AC_CHECK_SIZEOF([long long], [bft])
BFT_AC_CHECK_SIZEOF([float],     [bft])
BFT_AC_CHECK_SIZEOF([double],    [bft])

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------

AC_CHECK_FUNCS([getrusage gettimeofday sbrk sysinfo])
AC_CHECK_FUNCS([mkdir])
AC_CHECK_FUNCS([snprintf])
AC_CHECK_FUNCS([posix_memalign])

#------------------------------------------------------------------------------
# Checks for system services.
#------------------------------------------------------------------------------

# Largefile support (may be disabled in case of conflicts)

AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Test for size of off_t once largefile support is activated

BFT_AC_CHECK_SIZEOF([off_t], [bft])

#------------------------------
# optional malloc hooks support (for glibc)
#------------------------------

have_malloc_hooks=no

AC_ARG_ENABLE(malloc_hooks,
  [  --enable-malloc-hooks  Use malloc hooks when available],
  [
    case "${enableval}" in
      yes) malloc_hooks=true ;;
      no)  malloc_hooks=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-malloc-hooks]) ;;
    esac
  ],
  [ malloc_hooks=false ]
)

saved_LIBS="$LIBS"

# If we use malloc hooks, we may need to add a -lmalloc link flag.

if test "x$malloc_hooks" = "xtrue" ; then

  # Basic test
  AC_MSG_CHECKING([for malloc hooks])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],
                 [[ int __malloc_hook(); __malloc_hook(); ]])],
                 [have_malloc_hooks=yes],
                 [have_malloc_hooks=no])
  AC_MSG_RESULT($have_malloc_hooks)

fi

if test "x$have_malloc_hooks" = "xyes"; then
  AC_DEFINE(HAVE_MALLOC_HOOKS, 1, HAVE_MALLOC_HOOKS)
fi

# Basic backtrace functionnality

have_glibc_backtrace=no

AC_MSG_CHECKING([for glibc backtrace])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#if defined(__GNUC__)
#define _GNU_SOURCE
#include <memory.h>
#include <execinfo.h>
#endif]],
               [[ void *tr_array[10]; backtrace(tr_array, 10); ]])],
               [have_glibc_backtrace=yes],
               [have_glibc_backtrace=no])
AC_MSG_RESULT($have_glibc_backtrace)

if test "x$have_glibc_backtrace" = "xyes"; then
  AC_DEFINE(HAVE_GLIBC_BACKTRACE, 1, HAVE_GLIBC_BACKTRACE)
fi

# C++ demangling

have_cplus_demangle=no

AC_MSG_CHECKING([for glibc cplus_demangle])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <demangle.h>]],
               [[ const char *s; cplus_demangle(s, auto_demangling); ]])],
               [have_cplus_demangle=yes],
               [have_cplus_demangle=no])
AC_MSG_RESULT($have_cplus_demangle)

if test "x$have_cplus_demangle" = "xyes"; then
  AC_DEFINE(HAVE_CPLUS_DEMANGLE, 1, HAVE_CPLUS_DEMANGLE)
fi

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

FVM_AC_CONFIG_PUBL_SET(FVM_HAVE_MPI, $have_mpi,
                       [FVM library configured with MPI support ?])
if test "x$have_mpi" = "xyes" ; then
  FVM_AC_CONFIG_PUBL_SET(FVM_HAVE_MPI_IO, $have_mpi_io,
                         [FVM library configured with MPI IO support ?])
  FVM_AC_CONFIG_PUBL_SET(FVM_HAVE_MPI_ONE_SIDED, $have_mpi_one_sided,
                         [FVM library configured with MPI One-sided support ?])
fi

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

#
# A voir si ca ne pose pas de probleme

AC_ARG_VAR(CPP_ADD,  [Additional user C compiler flags])
AC_ARG_VAR(CC_ADD,   [Additional user C compiler flags])
AC_ARG_VAR(LD_ADD,   [Additional user linker flags, e.g. -L<libdir>])
AC_ARG_VAR(LIBS_ADD, [Additional user librairies, e.g. -l<lib>])

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
BFT_AC_CONFIG_INFO_VERSION([bft_version_string])
FVM_AC_CONFIG_INFO_VERSION([fvm_version_string])

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
AC_CONFIG_FILES([Makefile src/Makefile src/bft/Makefile src/fvm/Makefile cython/Makefile tests/Makefile])

#------------------------------------------------------------------------------

CWIPI_AC_CONFIG_PUBL_FINALIZE
FVM_AC_CONFIG_PUBL_FINALIZE
BFT_AC_CONFIG_PUBL_FINALIZE

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

#------------------------------------------------------------------------------

echo
echo "Configuration options:"
echo " use debugging code: "$debug""
echo " MPI (Message Passing Interface) support: "$have_mpi""
if (test x$have_mpi = xyes) ; then
  echo "   MPI I/O support: "$have_mpi_io""
  echo "   MPI2 one-sided communication support: "$have_mpi_one_sided""
fi
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the CWIPI package"
echo

