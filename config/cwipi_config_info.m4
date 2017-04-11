
# CWIPI_AC_CONFIG_INFO_INIT([OUTPUT FILE NAME])
#------------------------
# Write main config info file header.

AC_DEFUN([CWIPI_AC_CONFIG_INFO_INIT], [

# First arg is output file name
if test "$1" = "" ; then
  cwipi_ac_config_info="cwipi-config"
else
  cwipi_ac_config_info=$1
fi

outfile="$cwipi_ac_config_info"-tmp

  cat > $outfile <<\_______EOF
#!/bin/sh

# This file is generated by the configure script.
_______EOF

AC_MSG_NOTICE([initializing $cwipi_ac_config_info])

])

# CWIPI_AC_CONFIG_INFO_EXTRA([extra config info strings])
#-------------------------
# Write extra info to config info file header.

AC_DEFUN([CWIPI_AC_CONFIG_INFO_EXTRA], [

AC_REQUIRE([CWIPI_AC_CONFIG_INFO_INIT])dnl

outfile="$cwipi_ac_config_info"-tmp

echo "$1" >> $outfile

])

# CWIPI_AC_CONFIG_INFO_CC([cc], [version], [version_full])
#----------------------
# Write available compiler info to config info file header.

AC_DEFUN([CWIPI_AC_CONFIG_INFO_CC], [

AC_REQUIRE([CWIPI_AC_CONFIG_INFO_INIT])dnl

outfile="$cwipi_ac_config_info"-tmp

if test "$2" != "" -o "$3" != "" ; then
  echo "" >> $outfile
  echo "# C compiler used for build: $2" >> $outfile
  echo "# --------------------------"    >> $outfile
  if test "$3" != "" ; then
    echo "$3" | sed 's/^/# /' >> $outfile
  fi
  echo "" >> $outfile
else
  AC_MSG_NOTICE([C compiler version info unavailable for configuration file])
fi

echo "cc=$1" >> $outfile
echo "" >> $outfile
])

# CWIPI_AC_CONFIG_INFO_VERSION([version])
#------------------------------
# Write version info.

AC_DEFUN([CWIPI_AC_CONFIG_INFO_VERSION], [

AC_REQUIRE([CWIPI_AC_CONFIG_INFO_INIT])dnl

outfile="$cwipi_ac_config_info"-tmp

echo ""                  >> $outfile
echo "# Package Version" >> $outfile
echo "#----------------" >> $outfile
echo ""                  >> $outfile

echo "version=\"$1\"" >> $outfile
echo "" >> $outfile

])

# CWIPI_AC_CONFIG_INFO_DIRS([prefix], [exec_prefix], [includedir], [libdir])
#------------------------
# Write install path info.

AC_DEFUN([CWIPI_AC_CONFIG_INFO_DIRS], [

AC_REQUIRE([CWIPI_AC_CONFIG_INFO_INIT])dnl

outfile="$cwipi_ac_config_info"-tmp

echo ""                           >> $outfile
echo "# Installation directories" >> $outfile
echo "#-------------------------" >> $outfile
echo ""                           >> $outfile

echo "prefix=\"$1\"" >> $outfile
if test "$2" = "NONE" ; then
  echo "exec_prefix=\"\${prefix}\"" >> $outfile
else
  echo "exec_prefix=\"$2\"" >> $outfile
fi
echo "includedir=\"$3\"" >> $outfile
echo "libdir=\"$4\"" >> $outfile
echo "" >> $outfile

])

# CWIPI_AC_CONFIG_INFO_FLAGS([cppflags], [cflags], [ldflags], [libs],
#                          [scope], [scope_comments])
#-------------------------
# Write compiler and linker flags for a given scope (default if scope is empty).

AC_DEFUN([CWIPI_AC_CONFIG_INFO_FLAGS], [

AC_REQUIRE([CWIPI_AC_CONFIG_INFO_INIT])dnl

outfile="$cwipi_ac_config_info"-tmp

if test "$5" != "" ; then
  cwipi_ac_config_info_scope="${cwipi_ac_config_info_scope} $5"
  scope_prefix="$5_"
else
  scope_prefix=""
fi

echo "${scope_prefix}cppflags=\"$1\"" >> $outfile
echo "${scope_prefix}cflags=\"$2\""   >> $outfile
echo "${scope_prefix}ldflags=\"$3\""  >> $outfile
echo "${scope_prefix}libs=\"$4\""     >> $outfile

echo "" >> $outfile

])

# CWIPI_AC_CONFIG_INFO_FINALIZE
#----------------------------
# Write main config info file body (script part).

AC_DEFUN([CWIPI_AC_CONFIG_INFO_FINALIZE], [

AC_REQUIRE([CWIPI_AC_CONFIG_INFO_INIT])dnl

outfile="$cwipi_ac_config_info"-tmp

cat >> $outfile <<\_______EOF

# Output info depending on command-line options
#----------------------------------------------

if test [$]# -eq 0; then
  show_help=yes

else
  for opt in [$]*
  do
    case [$]opt in
      --prefix)
        output="${output} $prefix"
        ;;
      --exec-prefix)
        if test "$exec_prefix" != "" ; then
          output="${output} $exec_prefix"
        else
          output="${output} $prefix/bin"
        fi
        ;;
      --includedir)
        if test "$includedir" != "" ; then
          output="${output} $includedir"
        else
          output="${output} $prefix/include"
        fi
        ;;
      --version)
        output="${output} $version"
        ;;
      --cppflags)
        echo_cppflags=yes
        ;;
      --cflags)
        echo_cflags=yes
        ;;
      --ldflags)
        echo_ldflags=yes
        ;;
      --libs)
        echo_libs=yes
        ;;
      --cc)
        output="${output} $cc"
        ;;
      -*)
        show_help=yes
        ;;
      *)
        scope=$opt
        ;;
    esac
  done
fi

if test "[$]show_help" = "yes"; then
  cat <<EOF
Usage: cwipi-config [options] [scope]

Options:
        --prefix            installation path prefix
        --exec-prefix       system-dependent path prefix
        --includedir        C header files path
        --version           library version

        --cppflags          C preprocessor flags (e.g. -D<macro>, ...)
        --cflags            C flags (e.g. -O, -g, ...)
        --ldflags           linker flags (e.g. -g, -L<path>, ...)
        --libs              librairies used (e.g. -l<libname>)

        --cc                C compiler used for build

Scope:
        use (default)       Options required for user code
_______EOF

  for scope in $cwipi_ac_config_info_scope ; do
    if test "$scope" = "build" ; then
      echo "        build               Options used for build" >> $outfile
    else
      echo "        $scope" >> $outfile
    fi
  done

cat >> $outfile <<\_______EOF
EOF
  exit 1
fi
_______EOF

if test "$cwipi_ac_config_info_scope" != "" ; then

  cat >> $outfile <<\_______EOF

case "[$]scope" in
_______EOF

  for scope in $cwipi_ac_config_info_scope ; do

    echo "  ${scope})"                      >> $outfile
    echo "    cppflags=\$${scope}_cppflags" >> $outfile
    echo "    cflags=\$${scope}_cflags"     >> $outfile
    echo "    ldflags=\$${scope}_ldflags"   >> $outfile
    echo "    libs=\$${scope}_libs"         >> $outfile
    echo "    ;;"                           >> $outfile

  done

  echo "esac" >> $outfile

fi

cat >> $outfile <<\_______EOF

if test "$echo_cppflags" = "yes" ; then
  output="${output} $cppflags"
fi
if test "$echo_cflags" = "yes" ; then
  output="${output} $cflags"
fi
if test "$echo_ldflags" = "yes" ; then
  output="${output} $ldflags"
fi
if test "$echo_libs" = "yes" ; then
  output="${output} $libs"
fi

if test "$output" != "" ; then
  echo $output
fi

_______EOF

AC_MSG_NOTICE([closing $cwipi_ac_config_info])

diff $outfile $cwipi_ac_config_info > /dev/null 2>&1
if test $? -eq 0 ; then
  AC_MSG_NOTICE([$cwipi_ac_config_info is unchanged])
  rm -f $outfile
else
  mv $outfile $cwipi_ac_config_info
  chmod +x $cwipi_ac_config_info
fi

])dnl
