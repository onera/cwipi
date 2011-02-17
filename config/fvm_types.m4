dnl Copyright (C) 2004-2006 EDF
dnl
dnl This file is part of the FVM software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl FVM source distribution.

# FVM_AC_CHECK_SIZEOF(TYPE, [PREFIX])
#------------------------------------
# get type sizes
# Optionnaly, the corresponding SIZEOF definition may be prefixed
# by the PREFIX variable

AC_DEFUN([FVM_AC_CHECK_SIZEOF],[

if test "$1" = "" ; then
  AC_MSG_ERROR([configure test cannot be run])
fi

AC_REQUIRE([FVM_AC_CONFIG_PUBL_INIT])dnl

fvm_ac_lcname=`echo "$1" | sed y/' *'/'_p'/`
fvm_ac_lower='abcdefghijklmnopqrstuvwxyz'
fvm_ac_upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
if test "$2" != "" ; then
  fvm_ac_szname_prefix=`echo $2 | sed y/$fvm_ac_lower/$fvm_ac_upper/`_
else
  fvm_ac_szname_prefix=""
fi
fvm_ac_szname_postfix=`echo $fvm_ac_lcname | sed y/$fvm_ac_lower/$fvm_ac_upper/`
fvm_ac_szname="${fvm_ac_szname_prefix}SIZEOF_${fvm_ac_szname_postfix}"
unset fvm_ac_lower
unset fvm_ac_upper
unset fvm_ac_szname_prefix
unset fvm_ac_szname_postfix

AC_CHECK_SIZEOF($1)
eval fvm_ac_sizeof=\$ac_cv_sizeof_$fvm_ac_lcname
if test "$fvm_ac_sizeof" != "" -a "$fvm_ac_sizeof" != "0"; then
  FVM_AC_CONFIG_PUBL_DEFINE([$fvm_ac_szname], [$fvm_ac_sizeof],
                            [The size of a '$1', as computed by sizeof.])
else
  FVM_AC_CONFIG_PUBL_SET([$fvm_ac_szname], [no],
                         [The size of a '$1', as computed by sizeof.])
fi

unset fvm_ac_lcname
unset fvm_ac_szname
unset fvm_ac_sizeof

/bin/rm -f conftest*])dnl

# FVM_AC_SET_GNUM()
#------------------
# Select fvm_gnum_t size

AC_DEFUN([FVM_AC_SET_GNUM],[

AC_REQUIRE([FVM_AC_CONFIG_PUBL_INIT])dnl

# Check if we use long global numbers

AC_ARG_ENABLE(long-gnum,
  [  --enable-long-gnum      use long global numbers],
  [
    case "${enableval}" in
      yes) fvm_have_long_gnum=yes ;;
      no)  fvm_have_long_gnum=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-long-gnum]) ;;
    esac
  ],
  [ fvm_have_long_gnum=no ]
)

FVM_AC_CONFIG_PUBL_SET(FVM_HAVE_LONG_GNUM, [$fvm_have_long_gnum],
                       [Use 64-bit type if available for fvm_gnum_t.])

unset fvm_have_long_gnum

])dnl

