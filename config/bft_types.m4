dnl Copyright (C) 2004 EDF
dnl
dnl This file is part of the BFT software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl BFT source distribution.

# BFTC_AC_CHECK_SIZEOF(TYPE, [PREFIX])
#------------------------------------
# get type sizes
# Optionnaly, the corresponding SIZEOF definition may be prefixed
# by the PREFIX variable

AC_DEFUN([BFTC_AC_CHECK_SIZEOF],[

if test "$1" = "" ; then
  AC_MSG_ERROR([configure test cannot be run])
fi

AC_REQUIRE([BFTC_AC_CONFIG_PUBL_INIT])dnl

bftc_ac_lcname=`echo "$1" | sed y/' *'/'_p'/`
bftc_ac_lower='abcdefghijklmnopqrstuvwxyz'
bftc_ac_upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
if test "$2" != "" ; then
  bftc_ac_szname_prefix=`echo $2 | sed y/$bftc_ac_lower/$bftc_ac_upper/`_
else
  bftc_ac_szname_prefix=""
fi
bftc_ac_szname_postfix=`echo $bftc_ac_lcname | sed y/$bftc_ac_lower/$bftc_ac_upper/`
bftc_ac_szname="${bftc_ac_szname_prefix}SIZEOF_${bftc_ac_szname_postfix}"
unset bftc_ac_lower
unset bftc_ac_upper
unset bftc_ac_szname_prefix
unset bftc_ac_szname_postfix

AC_CHECK_SIZEOF($1)
eval bftc_ac_sizeof=\$ac_cv_sizeof_$bftc_ac_lcname
if test "$bftc_ac_sizeof" != "" ; then
  BFTC_AC_CONFIG_PUBL_DEFINE([$bftc_ac_szname], [$bftc_ac_sizeof],
                            [The size of a '$1', as computed by sizeof.])
else
  BFTC_AC_CONFIG_PUBL_SET([$bftc_ac_szname], [no],
                         [The size of a '$1', as computed by sizeof.])
fi

unset bftc_ac_lcname
unset bftc_ac_szname
unset bftc_ac_sizeof

/bin/rm -f conftest*])dnl

