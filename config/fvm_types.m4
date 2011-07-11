dnl Copyright (C) 2004-2006 EDF
dnl
dnl This file is part of the FVM software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl FVM source distribution.

# FVMC_AC_CHECK_SIZEOF(TYPE, [PREFIX])
#------------------------------------
# get type sizes
# Optionnaly, the corresponding SIZEOF definition may be prefixed
# by the PREFIX variable

AC_DEFUN([FVMC_AC_CHECK_SIZEOF],[

if test "$1" = "" ; then
  AC_MSG_ERROR([configure test cannot be run])
fi

AC_REQUIRE([FVMC_AC_CONFIG_PUBL_INIT])dnl

fvmc_ac_lcname=`echo "$1" | sed y/' *'/'_p'/`
fvmc_ac_lower='abcdefghijklmnopqrstuvwxyz'
fvmc_ac_upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
if test "$2" != "" ; then
  fvmc_ac_szname_prefix=`echo $2 | sed y/$fvmc_ac_lower/$fvmc_ac_upper/`_
else
  fvmc_ac_szname_prefix=""
fi
fvmc_ac_szname_postfix=`echo $fvmc_ac_lcname | sed y/$fvmc_ac_lower/$fvmc_ac_upper/`
fvmc_ac_szname="${fvmc_ac_szname_prefix}SIZEOF_${fvmc_ac_szname_postfix}"
unset fvmc_ac_lower
unset fvmc_ac_upper
unset fvmc_ac_szname_prefix
unset fvmc_ac_szname_postfix

AC_CHECK_SIZEOF($1)
eval fvmc_ac_sizeof=\$ac_cv_sizeof_$fvmc_ac_lcname
if test "$fvmc_ac_sizeof" != "" -a "$fvmc_ac_sizeof" != "0"; then
  FVMC_AC_CONFIG_PUBL_DEFINE([$fvmc_ac_szname], [$fvmc_ac_sizeof],
                            [The size of a '$1', as computed by sizeof.])
else
  FVMC_AC_CONFIG_PUBL_SET([$fvmc_ac_szname], [no],
                         [The size of a '$1', as computed by sizeof.])
fi

unset fvmc_ac_lcname
unset fvmc_ac_szname
unset fvmc_ac_sizeof

/bin/rm -f conftest*])dnl

# FVMC_AC_SET_GNUM()
#------------------
# Select fvmc_gnum_t size

AC_DEFUN([FVMC_AC_SET_GNUM],[

AC_REQUIRE([FVMC_AC_CONFIG_PUBL_INIT])dnl

# Check if we use long global numbers

AC_ARG_ENABLE(long-gnum,
  [  --enable-long-gnum      use long global numbers],
  [
    case "${enableval}" in
      yes) fvmc_have_long_gnum=yes ;;
      no)  fvmc_have_long_gnum=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-long-gnum]) ;;
    esac
  ],
  [ fvmc_have_long_gnum=no ]
)

FVMC_AC_CONFIG_PUBL_SET(FVMC_HAVE_LONG_GNUM, [$fvmc_have_long_gnum],
                       [Use 64-bit type if available for fvmc_gnum_t.])

unset fvmc_have_long_gnum

])dnl

