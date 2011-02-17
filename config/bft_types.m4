dnl Copyright (C) 2004 EDF
dnl
dnl This file is part of the BFT software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl BFT source distribution.

# BFT_AC_CHECK_SIZEOF(TYPE, [PREFIX])
#------------------------------------
# get type sizes
# Optionnaly, the corresponding SIZEOF definition may be prefixed
# by the PREFIX variable

AC_DEFUN([BFT_AC_CHECK_SIZEOF],[

if test "$1" = "" ; then
  AC_MSG_ERROR([configure test cannot be run])
fi

AC_REQUIRE([BFT_AC_CONFIG_PUBL_INIT])dnl

bft_ac_lcname=`echo "$1" | sed y/' *'/'_p'/`
bft_ac_lower='abcdefghijklmnopqrstuvwxyz'
bft_ac_upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
if test "$2" != "" ; then
  bft_ac_szname_prefix=`echo $2 | sed y/$bft_ac_lower/$bft_ac_upper/`_
else
  bft_ac_szname_prefix=""
fi
bft_ac_szname_postfix=`echo $bft_ac_lcname | sed y/$bft_ac_lower/$bft_ac_upper/`
bft_ac_szname="${bft_ac_szname_prefix}SIZEOF_${bft_ac_szname_postfix}"
unset bft_ac_lower
unset bft_ac_upper
unset bft_ac_szname_prefix
unset bft_ac_szname_postfix

AC_CHECK_SIZEOF($1)
eval bft_ac_sizeof=\$ac_cv_sizeof_$bft_ac_lcname
if test "$bft_ac_sizeof" != "" ; then
  BFT_AC_CONFIG_PUBL_DEFINE([$bft_ac_szname], [$bft_ac_sizeof],
                            [The size of a '$1', as computed by sizeof.])
else
  BFT_AC_CONFIG_PUBL_SET([$bft_ac_szname], [no],
                         [The size of a '$1', as computed by sizeof.])
fi

unset bft_ac_lcname
unset bft_ac_szname
unset bft_ac_sizeof

/bin/rm -f conftest*])dnl

