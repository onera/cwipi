

# CWIPI_AC_CONFIG_PUBL_INIT([OUPUT FILE NAME])
#-------------------------------------------
# Initialize file

AC_DEFUN([CWIPI_AC_CONFIG_PUBL_INIT],[

# First arg is output file name
if test "$1" = "" ; then
  cwipi_ac_config_publ_h="config_publ.h"
else
  cwipi_ac_config_publ_h=$1
fi
cwipi_ac_lower='abcdefghijklmnopqrstuvwxyz.'
cwipi_ac_upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ_'
cwipi_ac_config_publ_upname=`echo $cwipi_ac_config_publ_h | sed y/$cwipi_ac_lower/$cwipi_ac_upper/`
unset cwipi_ac_lower
unset cwipi_ac_upper
echo "#ifndef __"${cwipi_ac_config_publ_upname}"__" >  "$cwipi_ac_config_publ_h"-tmp
echo "#define __"${cwipi_ac_config_publ_upname}"__" >> "$cwipi_ac_config_publ_h"-tmp
echo >> "$cwipi_ac_config_publ_h"-tmp
echo "/* $cwipi_ac_config_publ_h. Generated by configure */" >> "$cwipi_ac_config_publ_h"-tmp
echo >> "$cwipi_ac_config_publ_h"-tmp

AC_MSG_NOTICE([initializing $cwipi_ac_config_publ_h])
])dnl


# CWIPI_AC_CONFIG_PUBL_FINALIZE
#----------------------------
# Finalize file

AC_DEFUN([CWIPI_AC_CONFIG_PUBL_FINALIZE],[

AC_REQUIRE([CWIPI_AC_CONFIG_PUBL_INIT])dnl

echo "#endif /* __"${cwipi_ac_config_publ_upname}"__ */" >>  "$cwipi_ac_config_publ_h"-tmp

AC_MSG_NOTICE([closing $cwipi_ac_config_publ_h])

diff $cwipi_ac_config_publ_h $cwipi_ac_config_publ_h-tmp  > /dev/null 2>&1
if test $? -eq 0 ; then
  AC_MSG_NOTICE([$cwipi_ac_config_publ_h  is unchanged])
  rm -f $cwipi_ac_config_publ_h-tmp
else
  mv $cwipi_ac_config_publ_h-tmp $cwipi_ac_config_publ_h
fi

unset cwipi_ac_config_publ_h
unset cwipi_ac_config_publ_upname

])dnl


# CWIPI_AC_CONFIG_PUBL_VERBATIM([VERBATIM TEXT])
#---------------------------------------------
# Add text to config file

AC_DEFUN([CWIPI_AC_CONFIG_PUBL_VERBATIM],[

AC_REQUIRE([CWIPI_AC_CONFIG_PUBL_INIT])dnl

echo "$1" >> "$cwipi_ac_config_publ_h"-tmp
echo ""   >> "$cwipi_ac_config_publ_h"-tmp

])dnl


# CWIPI_AC_CONFIG_PUBL_DEFINE(VARIABLE NAME, VALUE, [COMMENT])
#-----------------------------------------------------------
# Define variable

AC_DEFUN([CWIPI_AC_CONFIG_PUBL_DEFINE],[

AC_REQUIRE([CWIPI_AC_CONFIG_PUBL_INIT])dnl

unset cwipi_ac_lower
unset cwipi_ac_upper
if test "$3" != "" ; then
  echo "/* $3 */" >> "$cwipi_ac_config_publ_h"-tmp
fi
echo "#define $1 $2" >> "$cwipi_ac_config_publ_h"-tmp
echo >> "$cwipi_ac_config_publ_h"-tmp

])dnl


# CWIPI_AC_CONFIG_PUBL_DEFINE_STRING(VARIABLE NAME, VALUE, [COMMENT])
#------------------------------------------------------------------
# Define string variable

AC_DEFUN([CWIPI_AC_CONFIG_PUBL_DEFINE_STRING],[

AC_REQUIRE([CWIPI_AC_CONFIG_PUBL_INIT])dnl

unset cwipi_ac_lower
unset cwipi_ac_upper
if test "$3" != "" ; then
  echo "/* $3 */" >> "$cwipi_ac_config_publ_h"-tmp
fi
echo '#define $1 "$2"' >> "$cwipi_ac_config_publ_h"-tmp
echo >> "$cwipi_ac_config_publ_h"-tmp

])dnl


# CWIPI_AC_CONFIG_PUBL_SET(VARIABLE NAME, YES OR NO, [COMMENT])
#------------------------------------------------------------
# Set variable (define to 1 if second arg is yes, undefine if no )

AC_DEFUN([CWIPI_AC_CONFIG_PUBL_SET],[

# First arg is variable name, second arg is value, optional third arg is comment

AC_REQUIRE([CWIPI_AC_CONFIG_PUBL_INIT])dnl

unset cwipi_ac_lower
unset cwipi_ac_upper
if test "$3" != "" ; then
  echo "/* $3 */" >> "$cwipi_ac_config_publ_h"-tmp
fi
if test "$2" = "yes" ; then
  echo "#define $1 1" >> "$cwipi_ac_config_publ_h"-tmp
else
  echo "#undef $1" >> "$cwipi_ac_config_publ_h"-tmp
fi
echo >> "$cwipi_ac_config_publ_h"-tmp

])dnl
