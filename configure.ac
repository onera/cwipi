#-----------------------------------------------------------------------------
# This file is part of the CWIPI library. 
#
# Copyright (C) 2016  ONERA
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

  Copyright (C) 2011-2016  ONERA

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

m4_define([bftc_licence_c_comment],
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

m4_define([fvmc_licence_c_comment],
[/*
  This file is part of the CWIPI library, derived from the "Finite Volume Mesh" 
  library, intended to provide finite volume mesh and associated 
  fields I/O and manipulation services.

  Copyright (C) 2011-2016  ONERA

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
m4_define([cwipi_minor_version], [9])
m4_define([cwipi_release_version], [6])
m4_define([cwipi_version_extra], [])
m4_define([cwipi_version_string],
          [cwipi_major_version.cwipi_minor_version.cwipi_release_version@&t@cwipi_version_extra])

m4_define([fvmc_major_version], [0])
m4_define([fvmc_minor_version], [15])
m4_define([fvmc_release_version], [1])
m4_define([fvmc_version_extra], [-onera-6])
m4_define([fvmc_version_string],
          [fvmc_major_version.fvmc_minor_version.fvmc_release_version@&t@fvmc_version_extra])

m4_define([bftc_major_version], [1])
m4_define([bftc_minor_version], [1])
m4_define([bftc_release_version], [3])
m4_define([bftc_version_string],
          [bftc_major_version.bftc_minor_version.bftc_release_version])

AC_INIT([cwipi],[cwipi_version_string],[eric.quemerais@onera.fr])
AC_CONFIG_SRCDIR([src/cwipi.h])

# Use the config directory for libtool stuff ... 
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_SUBST(ac_aux_dir)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([config_priv.h:config_priv.h.in])
#AC_CONFIG_HEADERS([fvmc_config_priv.h:fvmc_config_priv.h.in])
#AC_CONFIG_HEADERS([bftc_config_priv.h:bftc_config_priv.h.in])

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

FVMC_AC_CONFIG_PUBL_INIT([fvmc_config.h])
FVMC_AC_CONFIG_PUBL_VERBATIM([fvmc_licence_c_comment])

BFTC_AC_CONFIG_PUBL_INIT([bftc_config.h])
BFTC_AC_CONFIG_PUBL_VERBATIM([bftc_licence_c_comment])


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

FVMC_AC_CONFIG_PUBL_DEFINE_STRING(FVMC_VERSION, fvmc_version_string,
                                 [Complete version number])

FVMC_AC_CONFIG_PUBL_DEFINE(FVMC_MAJOR_VERSION, fvmc_major_version,
                          [Major version number])
FVMC_AC_CONFIG_PUBL_DEFINE(FVMC_MINOR_VERSION, fvmc_minor_version,
                          [Minor version number])
FVMC_AC_CONFIG_PUBL_DEFINE(FVMC_RELEASE_VERSION, fvmc_release_version,
                          [Release version number])

BFTC_AC_CONFIG_PUBL_DEFINE_STRING(BFTC_VERSION, bftc_version_string,
                                 [Complete version number])

BFTC_AC_CONFIG_PUBL_DEFINE(BFTC_MAJOR_VERSION, bftc_major_version,
                          [Major version number])
BFTC_AC_CONFIG_PUBL_DEFINE(BFTC_MINOR_VERSION, bftc_minor_version,
                          [Minor version number])
BFTC_AC_CONFIG_PUBL_DEFINE(BFTC_RELEASE_VERSION, bftc_release_version,
                          [Release version number])

#BFTC_AC_CONFIG_PUBL_DEFINE(BFTC_CPPCALLER, 1,
#                          [Active namespace bft for cwipi])
#FVMC_AC_CONFIG_PUBL_DEFINE(FVMC_CPPCALLER, 1,
#                          [Active namespace fvm for cwipi])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_LDFLAGS=$LDFLAGS

LT_INIT

AC_PROG_CC(icc gcc)

AC_PROG_CXX(icpc g++)

AC_PROG_FC(ifort gfortran)

LT_PATH_LD([icpc g++])

AC_PROG_INSTALL

AC_PROG_LN_S

#_LT_SET_OPTION([LT_INIT],[dlopen])
AC_DIAGNOSE([obsolete],[_LT_SET_OPTION([LT_INIT],[dlopen])
AC_DIAGNOSE([obsolete],[AC_LIBTOOL_DLOPEN: Remove this warning and the call to _LT_SET_OPTION when you
put the 'dlopen' option into LT_INIT's first parameter.])
: Remove this warning and the call to _LT_SET_OPTION when you
put the `dlopen' option into LT_INIT's first parameter.])

#LT_INIT

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
# Checks for optional Python bindings 
#------------------------------------------------------------------------------

# Check for optional Python bindings 

AC_ARG_ENABLE(python-bindings,
  [AS_HELP_STRING([--enable-python-bindings],[build with Python interface])],
  [
    case "${enableval}" in
      yes) python_bindings=true ;;
      no)  python_bindings=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-python-bindings]) ;;
    esac
  ],
  [ python_bindings=false ]
)
AM_CONDITIONAL(ENABLE_PYTHON, test x$python_bindings = xtrue)

AC_ARG_ENABLE(rebuild-bindings,
  [AS_HELP_STRING([--enable-rebuild-bindings],[rebuild Cython bindings])],
  [
    case "${enableval}" in
      yes) rebuild_bindings=true ;;
      no)  rebuild_bindings=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-rebuild-bindings]) ;;
    esac
  ],
  [ rebuild_bindings=false ]
)

AM_CONDITIONAL(ENABLE_CYTHON, test x$rebuild_bindings = xtrue)

AC_ARG_WITH(python-exec, [AS_HELP_STRING([--with-python-exec=PATH],
                                         [specify prefix directory for a Python interpreter])],
                         [PYTHON_EXE=${withval}],
                         [PYTHON_EXE=none])dnl

AC_ARG_WITH(cython-exec, [AS_HELP_STRING([--with-cython-exec=PATH],
                                         [specify prefix directory for a Cython interpreter])],
                         [CYTHON_EXE=${withval}],
                         [CYTHON_EXE=none])dnl

if test "x${python_bindings}" != "xfalse"; then

    cp -f ${srcdir}/Cython/cwipi.c.default  ${srcdir}/Cython/cwipi.c

    #------------------
    # Checks for Python
    #------------------

    if test "x${PYTHON_EXE}" = "xnone" ; then
        AC_PATH_PROG([PYEXE], [python], ["no"])
    else
        AC_PATH_PROG([PYEXE], [python], ["no"], [${PYTHON_EXE}])
    fi

    python_version_min="2.7"
    mpi4py_version_min="1.2.2"
    numpy_version_min="1.5.1"
    cython_version_min="0.14.1"

    have_python="yes"
    have_python_binding="yes"

    have_python_include="no"
    have_numpy="no"
    have_mpi4py="no"
    have_cython="no"

    PYTHON=${PYEXE}
    if test "x${PYTHON}" = "xno" ; then
         have_python_binding="no"
         have_python="no"
    fi

    #
    # Check python version (for rebuild_bindings = false)

    if test "x${rebuild_bindings}" = "xfalse"; then    
       AM_PATH_PYTHON(,, PYTHON="no")
    else
       AM_PATH_PYTHON([${python_version_min}],, PYTHON="no")
    fi

    if test "x${PYTHON}" = "xno" ; then
       have_python_binding="no"
       have_python="no"
    fi

    if test "x${have_python}" = "xno" ; then
        AC_MSG_WARN([])
        AC_MSG_WARN([Python interpreter not found])
        AC_MSG_WARN([Please install Python from http://www.python.org])
        AC_MSG_WARN([])
        AC_MSG_WARN([Sorry, we can't continue without Python interpreter])
        AM_CONDITIONAL(ENABLE_PYTHON, false)

    else
    
        PYTHON_INCLUDE=`${PYTHON} -c 'from distutils import sysconfig; print( sysconfig.get_python_inc() )'`
        have_python_include="yes"        
        if test "x${PYTHON_INCLUDE}" = "x" ; then
            AC_MSG_WARN([])
            AC_MSG_WARN([Python include not found])
            AC_MSG_WARN([Please install Python from http://www.python.org])
            AC_MSG_WARN([])
            AC_MSG_WARN([Sorry, we can't continue without Python include])
            have_python_binding="no"
            have_python_include="no"        
        fi

        if test "x${have_python}" = "xyes" ; then

          MPI4PY_INCLUDE=`${PYTHON} -c 'import mpi4py; print( mpi4py.get_include() )'`
          have_mpi4py="yes"

          if test "x${MPI4PY_INCLUDE}" = "x" ; then
            AC_MSG_WARN([])
            AC_MSG_WARN([mpi4py package not found])
            AC_MSG_WARN([Please install mpi4py from http://code.google.com/p/mpi4py])
            AC_MSG_WARN([])
            AC_MSG_WARN([Sorry, we can't continue without mpi4py])
            have_python_binding="no"        
            have_mpi4py="no"
          else
           if test "x${rebuild_bindings}" = "xfalse"; then    
             MPI4PY_VERSION=`${PYTHON} -c 'import mpi4py; print( mpi4py.__version__ )'`

             mpi4py_version_major_min=`echo "$mpi4py_version_min" | cut -f1 -d.`
             mpi4py_version_minor_min=`echo "$mpi4py_version_min" | cut -f2 -d.`
             mpi4py_version_release_min=`echo "$mpi4py_version_min" | cut -f3 -d.`

             mpi4py_version_major=`echo "$MPI4PY_VERSION" | cut -f1 -d.`
             mpi4py_version_minor=`echo "$MPI4PY_VERSION" | cut -f2 -d.`
             mpi4py_version_release=`echo "$MPI4PY_VERSION" | cut -f3 -d.`
AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],
[[#if $mpi4py_version_major < $mpi4py_version_major_min
#  error mpi4py major version < $mpi4py_version_major_min
#elif $mpi4py_version_major == $mpi4py_version_major_min
#  if $mpi4py_version_minor < $mpi4py_version_minor_min
#    error mpi4py minor version < $mpi4py_version_minor_min
#  elif $mpi4py_version_minor == $mpi4py_version_minor_min
#    if $mpi4py_version_release < $mpi4py_version_release_min
#      error mpi4py release version < $mpi4py_version_release_min
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible mpi4py version found])],
               [AC_MSG_WARN([compatible mpi4py version not found (min version ${mpi4py_version_min})])
                have_python_binding="no"
                have_mpi4py="no"
                ])
           fi 
          fi
        fi

        if test "x${have_python}" = "xyes" ; then

          NUMPY_INCLUDE=`${PYTHON} -c 'import numpy; print( numpy.get_include() )'`
          have_numpy="yes"

          if test "x${NUMPY_INCLUDE}" = "x" ; then
            AC_MSG_WARN([])
            AC_MSG_WARN([numpy package not found])
            AC_MSG_WARN([Please install numpy from http://numpy.scipy.org])
            AC_MSG_WARN([])
            AC_MSG_WARN([Sorry, we can't continue without numpy])
            have_python_binding="no"        
            have_numpy="no"
          else
           if test "x${rebuild_bindings}" = "xfalse"; then    
             NUMPY_VERSION=`${PYTHON} -c 'import numpy; print( numpy.__version__ )'`

             numpy_version_major_min=`echo "$numpy_version_min" | cut -f1 -d.`
             numpy_version_minor_min=`echo "$numpy_version_min" | cut -f2 -d.`
             numpy_version_release_min=`echo "$numpy_version_min" | cut -f3 -d.`

             numpy_version_major=`echo "$NUMPY_VERSION" | cut -f1 -d.`
             numpy_version_minor=`echo "$NUMPY_VERSION" | cut -f2 -d.`
             numpy_version_release=`echo "$NUMPY_VERSION" | cut -f3 -d.`
AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],
[[#if $numpy_version_major < $numpy_version_major_min
#  error numpy major version < $numpy_version_major_min
#elif $numpy_version_major == $numpy_version_major_min
#  if $numpy_version_minor < $numpy_version_minor_min
#    error numpy minor version < $numpy_version_minor_min
#  elif $numpy_version_minor == $numpy_version_minor_min
#    if $numpy_version_release < $numpy_version_release_min
#      error numpy release version < $numpy_version_release_min
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible numpy version found])],
               [AC_MSG_WARN([compatible numpy version not found (min version ${numpy_version_min})])
                have_python_binding="no"
                have_numpy="no"
                ])
           fi
          fi
        fi

        if test "x${have_python_binding}" = "xno" ; then
            AM_CONDITIONAL(ENABLE_PYTHON, false)
        fi

        if test "x${have_python}" = "xyes" ; then
           if test "x${rebuild_bindings}" = "xtrue"; then    

             if test "x${CYTHON_EXE}" = "xnone" ; then
               AC_PATH_PROG([CYEXE], [cython], ["no"])
             else
               AC_PATH_PROG([CYEXE], [cython], ["no"], [${CYTHON_EXE}])
             fi

             if test "x${CYEXE}" = "xno" ; then
               AC_MSG_WARN([])
               AC_MSG_WARN([Cython package not found])
               AC_MSG_WARN([Please install Cython from http://http://www.cython.org or remove --enable-rebuild-bindings option])
               AC_MSG_WARN([])
               AC_MSG_WARN([Sorry, we can't continue without cython])
               have_python_binding="no"        
               have_cython="no"
             else
               have_cython="yes"
               CYTHON_VERSION=`${CYEXE} -V 2>&1 | awk '{print $3}'`

               cython_version_major_min=`echo "$cython_version_min" | cut -f1 -d.`
               cython_version_minor_min=`echo "$cython_version_min" | cut -f2 -d.`
               cython_version_release_min=`echo "$cython_version_min" | cut -f3 -d.`

               cython_version_major=`echo "$CYTHON_VERSION" | cut -f1 -d.`
               cython_version_minor=`echo "$CYTHON_VERSION" | cut -f2 -d.`
               cython_version_release=`echo "$CYTHON_VERSION" | cut -f3 -d.`
AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],
[[#if $cython_version_major < $cython_version_major_min
#  error cython major version < $cython_version_major_min
#elif $cython_version_major == $cython_version_major_min
#  if $cython_version_minor < $cython_version_minor_min
#    error cython minor version < $cython_version_minor_min
#  elif $cython_version_minor == $cython_version_minor_min
#    if $cython_version_release < $cython_version_release_min
#      error cython release version < $cython_version_release_min
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible cython version found])],
               [AC_MSG_WARN([compatible cython version not found (min version ${cython_version_min})])
                have_python_binding="no"
                have_cython="no"
                ])
                if test "x${have_cython}" = "xyes" ; then
                  echo "running cython : " ${CYEXE} -I${MPI4PY_INCLUDE} -I${srcdir}/src -I${NUMPY_INCLUDE} ${srcdir}/Cython/cwipi.pyx
                  cp ${srcdir}/Cython/cwipi.pyx cwipi.pyx
                  error=
                  error=`${CYEXE} -v -I${MPI4PY_INCLUDE} -I${srcdir}/src -I${NUMPY_INCLUDE} cwipi.pyx 2>&1 | grep -i Err`
                  if test "x${error}" != "x" ; then
                     echo $error
                     AC_MSG_ERROR([Error during cython binding])
                  fi
                  cp -f cwipi.c ${srcdir}/Cython/cwipi.c
                fi
             fi
           fi
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

#LT_INIT

#--------------------------------------------------------------------------
# If the system is of a known type, declare it in the configuration header.
#--------------------------------------------------------------------------

case "$host_os" in
  linux*)
    AC_DEFINE(bftc_OS_Linux, 1, bftc_OS_Linux)
    ;;
  solaris2.*)
    AC_DEFINE(bftc_OS_Solaris, 1, bftc_OS_SunOS)
    ;;
  aix*)
    AC_DEFINE(bftc_OS_AIX, 1, bftc_OS_AIX)
    ;;
  osf*)
    AC_DEFINE(bftc_OS_OSF1, 1, bftc_OS_OSF1)
    ;;
  irix5.*|irix6.*)
    AC_DEFINE(bftc_OS_IRIX64, 1, bftc_OS_IRIX64)
    ;;
  hpux*)
    AC_DEFINE(bftc_OS_HP_UX, 1, bftc_OS_HP_UX)
    ;;
  uxpv*)
    AC_DEFINE(bftc_OS_UXPV, 1, bftc_OS_UNIX_System_V)
    ;;
esac

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

# Check for optionnal HDF support

#FVMC_AC_TEST_HDF5

# Check for optionnal CGNS support

#FVMC_AC_TEST_CGNS

# Check for optionnal MED support

#FVMC_AC_TEST_MED

# Optional Zlib Support

#BFTC_AC_ZLIB

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

FVMC_AC_CHECK_SIZEOF([void *], [fvmc])
FVMC_AC_CHECK_SIZEOF([short], [fvmc])
FVMC_AC_CHECK_SIZEOF([int], [fvmc])
FVMC_AC_CHECK_SIZEOF([long], [fvmc])
FVMC_AC_CHECK_SIZEOF([long long], [fvmc])
FVMC_AC_CHECK_SIZEOF([float], [fvmc])
FVMC_AC_CHECK_SIZEOF([double], [fvmc])

BFTC_AC_CHECK_SIZEOF([void *],    [bftc])
BFTC_AC_CHECK_SIZEOF([short],     [bftc])
BFTC_AC_CHECK_SIZEOF([int],       [bftc])
BFTC_AC_CHECK_SIZEOF([long],      [bftc])
BFTC_AC_CHECK_SIZEOF([long long], [bftc])
BFTC_AC_CHECK_SIZEOF([float],     [bftc])
BFTC_AC_CHECK_SIZEOF([double],    [bftc])

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

BFTC_AC_CHECK_SIZEOF([off_t], [bftc])

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

FVMC_AC_CONFIG_PUBL_SET(FVMC_HAVE_MPI, $have_mpi,
                       [FVM library configured with MPI support ?])
if test "x$have_mpi" = "xyes" ; then
  FVMC_AC_CONFIG_PUBL_SET(FVMC_HAVE_MPI_IO, $have_mpi_io,
                         [FVM library configured with MPI IO support ?])
  FVMC_AC_CONFIG_PUBL_SET(FVMC_HAVE_MPI_ONE_SIDED, $have_mpi_one_sided,
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
BFTC_AC_CONFIG_INFO_VERSION([bftc_version_string])
FVMC_AC_CONFIG_INFO_VERSION([fvmc_version_string])

CWIPI_AC_CONFIG_INFO_DIRS([$prefix], [$exec_prefix], [$includedir], [$libdir])

CWIPI_AC_CONFIG_INFO_FLAGS([$user_cppflags], [], [$user_ldflags], [$LIBS])

CWIPI_AC_CONFIG_INFO_FLAGS([$CPPFLAGS], [$CFLAGS], [$LDFLAGS], [$LIBS],
                         [build])

CWIPI_AC_CONFIG_INFO_FLAGS([$BFTC_CPPFLAGS], [], [$BFTC_LDFLAGS], [$BFTC_LIBS],
                         [bftc])

CWIPI_AC_CONFIG_INFO_FLAGS([$FVMC_CPPFLAGS], [], [$FVMC_LDFLAGS], [$FVMC_LIBS],
                         [fvmc])

CWIPI_AC_CONFIG_INFO_FLAGS([$MPI_CPPFLAGS], [], [$MPI_LDFLAGS], [$MPI_LIBS],
                         [mpi])

CWIPI_AC_CONFIG_INFO_FINALIZE

#------------------------------------------------------------------------------

# If we use GNU gettext, we need intl/Makefile, po/Makefile, and m4/Makefile.

# AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile doc/Makefile
#                  intl/Makefile po/Makefile.in m4/Makefile])

AC_CONFIG_FILES([Makefile src/Makefile src/bft/Makefile src/fvm/Makefile Cython/Makefile tests/Makefile])

#------------------------------------------------------------------------------

CWIPI_AC_CONFIG_PUBL_FINALIZE
FVMC_AC_CONFIG_PUBL_FINALIZE
BFTC_AC_CONFIG_PUBL_FINALIZE

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
echo " enable python bindings : "${python_bindings}""
if test "x${python_bindings}" != "xfalse"; then
  if test "x${rebuild_bindings}" = "xtrue"; then    
    echo "   Rebuilding Python bindings with cython ("${CYTHON_VERSION}")"
    echo "     Python interpreter (${PYTHON_VERSION}): "$have_python""
    echo "     Python include             : "$have_python_include""
    echo "     mpi4py (>= ${mpi4py_version_min})          : "$have_mpi4py""
    echo "     numpy (>= ${numpy_version_min})           : "$have_numpy""
  else
    echo "   Python bindings support : "$have_python_binding""
    echo "     Python interpreter (>= ${python_version_min}): "$have_python""
    echo "     Python include             : "$have_python_include""
    echo "     mpi4py (>= ${mpi4py_version_min})          : "$have_mpi4py""
    echo "     numpy (>= ${numpy_version_min})           : "$have_numpy""
  fi
fi
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

