dnl Copyright (C) 2004 EDF
dnl
dnl This file is part of the BFT software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl BFT source distribution.

# BFTC_AC_CONFIG_PUBL_INIT([OUPUT FILE NAME])
#-------------------------------------------
# Initialize file

AC_DEFUN([BFTC_AC_CONFIG_PUBL_INIT],[

# First arg is output file name
if test "$1" = "" ; then
  bftc_ac_config_publ_h="config_publ.h"
else
  bftc_ac_config_publ_h=$1
fi
bftc_ac_lower='abcdefghijklmnopqrstuvwxyz.'
bftc_ac_upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ_'
bftc_ac_config_publ_upname=`echo $bftc_ac_config_publ_h | sed y/$bftc_ac_lower/$bftc_ac_upper/`
unset bftc_ac_lower
unset bftc_ac_upper
echo "#ifndef __"${bftc_ac_config_publ_upname}"__" >  "$bftc_ac_config_publ_h"-tmp
echo "#define __"${bftc_ac_config_publ_upname}"__" >> "$bftc_ac_config_publ_h"-tmp
echo >> "$bftc_ac_config_publ_h"-tmp
echo "/* $bftc_ac_config_publ_h. Generated by configure */" >> "$bftc_ac_config_publ_h"-tmp
echo >> "$bftc_ac_config_publ_h"-tmp

AC_MSG_NOTICE([initializing $bftc_ac_config_publ_h])
])dnl


# BFTC_AC_CONFIG_PUBL_FINALIZE
#----------------------------
# Finalize file

AC_DEFUN([BFTC_AC_CONFIG_PUBL_FINALIZE],[

AC_REQUIRE([BFTC_AC_CONFIG_PUBL_INIT])dnl

echo "#endif /* __"${bftc_ac_config_publ_upname}"__ */" >>  "$bftc_ac_config_publ_h"-tmp

AC_MSG_NOTICE([closing $bftc_ac_config_publ_h])

diff $bftc_ac_config_publ_h $bftc_ac_config_publ_h-tmp  > /dev/null 2>&1
if test $? -eq 0 ; then
  AC_MSG_NOTICE([$bftc_ac_config_publ_h  is unchanged])
  rm -f $bftc_ac_config_publ_h-tmp
else
  mv $bftc_ac_config_publ_h-tmp $bftc_ac_config_publ_h
fi

unset bftc_ac_config_publ_h
unset bftc_ac_config_publ_upname

])dnl


# BFTC_AC_CONFIG_PUBL_VERBATIM([VERBATIM TEXT])
#---------------------------------------------
# Add text to config file

AC_DEFUN([BFTC_AC_CONFIG_PUBL_VERBATIM],[

AC_REQUIRE([BFTC_AC_CONFIG_PUBL_INIT])dnl

echo "$1" >> "$bftc_ac_config_publ_h"-tmp
echo ""   >> "$bftc_ac_config_publ_h"-tmp

])dnl


# BFTC_AC_CONFIG_PUBL_DEFINE(VARIABLE NAME, VALUE, [COMMENT])
#-----------------------------------------------------------
# Define variable

AC_DEFUN([BFTC_AC_CONFIG_PUBL_DEFINE],[

AC_REQUIRE([BFTC_AC_CONFIG_PUBL_INIT])dnl

unset bftc_ac_lower
unset bftc_ac_upper
if test "$3" != "" ; then
  echo "/* $3 */" >> "$bftc_ac_config_publ_h"-tmp
fi
echo "#define $1 $2" >> "$bftc_ac_config_publ_h"-tmp
echo >> "$bftc_ac_config_publ_h"-tmp

])dnl


# BFTC_AC_CONFIG_PUBL_DEFINE_STRING(VARIABLE NAME, VALUE, [COMMENT])
#------------------------------------------------------------------
# Define string variable

AC_DEFUN([BFTC_AC_CONFIG_PUBL_DEFINE_STRING],[

AC_REQUIRE([BFTC_AC_CONFIG_PUBL_INIT])dnl

unset bftc_ac_lower
unset bftc_ac_upper
if test "$3" != "" ; then
  echo "/* $3 */" >> "$bftc_ac_config_publ_h"-tmp
fi
echo '#define $1 "$2"' >> "$bftc_ac_config_publ_h"-tmp
echo >> "$bftc_ac_config_publ_h"-tmp

])dnl


# BFTC_AC_CONFIG_PUBL_SET(VARIABLE NAME, YES OR NO, [COMMENT])
#------------------------------------------------------------
# Set variable (define to 1 if second arg is yes, undefine if no )

AC_DEFUN([BFTC_AC_CONFIG_PUBL_SET],[

# First arg is variable name, second arg is value, optional third arg is comment

AC_REQUIRE([BFTC_AC_CONFIG_PUBL_INIT])dnl

unset bftc_ac_lower
unset bftc_ac_upper
if test "$3" != "" ; then
  echo "/* $3 */" >> "$bftc_ac_config_publ_h"-tmp
fi
if test "$2" = "yes" ; then
  echo "#define $1 1" >> "$bftc_ac_config_publ_h"-tmp
else
  echo "#undef $1" >> "$bftc_ac_config_publ_h"-tmp
fi
echo >> "$bftc_ac_config_publ_h"-tmp

])dnl
