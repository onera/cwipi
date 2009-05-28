# Shell script

# This file should be sourced by configure, and sets the following
# environment variables corresponding to the recommended settings for a
# given OS/CPU/compiler combination:
#
# cppflags_default       # Base CPPFLAGS                     (default: "")
#
# cflags_default         # Base CFLAGS                       (default: "")
# cflags_default_dbg     # Added to $CFLAGS for debugging    (default: "-g")
# cflags_default_opt     # Added to $CFLAGS for optimization (default: "-O")
# cflags_default_prf     # Added to $CFLAGS for profiling    (default: "")
#
# ldflags_default        # Base LDFLAGS                       (default: "")
# ldflags_default_dbg    # Added to $LDFLAGS for debugging    (default: "-g")
# ldflags_default_opt    # Added to $LDFLAGS for optimization (default: "-O")
# ldflags_default_prf    # Added to $LDFLAGS for profiling    (default: "")
#
# couplings_disable_shared     # Disable shared librairies          (default: "")

# Two other environment variable strings are defined, containing possibly
# more detailed compiler information:
#
# couplings_ac_cc_version      # Compiler version string, 1 line max.
# couplings_ac_cc_version_full # Compiler version string, 10 lines max.

# The sourcing approach and some tests are borrowed from the HDF5 configure
# environment.
#
# We choose to source this script rather than use a more classical m4 macro
# for this functionality, so that a user may more easily modify
# default compiler options or port to a new machine without requiring
# any advanced knowledge of autoconf or m4 macros, or installation of
# an autoconf environment on the target machine.

# Initialize local variables
#---------------------------

outfile=couplings_ac_cc_env-tmp

couplings_ac_cc_version=unknown

# Some compilers require a file to compile even for version info.

cat > conftest.c <<\_______EOF
int main()
{
  return 0;
}
_______EOF

# Compiler version info may be localization dependent (especially for gcc)

save_LANG=$LANG
unset LANG;

# Default pre-processor flags (not too dependent on compiler)
#----------------------------

case "$host_os" in
  *)
    cppflags_default=""
    ;;
esac

# Are we using gcc ?
#-------------------

couplings_gcc=no

if test "x$GCC" = "xyes"; then

  # Intel compiler passes as GCC but may be recognized by version string
  if test -n "`$CC --version | grep icc`" ; then
    couplings_gcc=icc
  else
    couplings_gcc=gcc
  fi

fi

if test "x$couplings_gcc" = "xgcc"; then

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
  couplings_ac_cc_version=`$CC --version 2>&1 | head -1`
  couplings_compiler_known=yes

  # Practical version info for option setting
  couplings_cc_version="`$CC -v 2>&1 |grep 'gcc version' |\
                  sed 's/.*gcc version \([-a-z0-9\.]*\).*/\1/'`"
  couplings_cc_vendor=`echo $couplings_cc_version |sed 's/\([a-z]*\).*/\1/'`
  couplings_cc_version=`echo $couplings_cc_version |sed 's/[-a-z]//g'`

  if test "x" = "x$couplings_cc_vendor" -a "x" != "x$couplings_cc_version"; then
    couplings_cc_vendor=gcc
  fi
  if test "-" != "$couplings_cc_vendor-$couplings_cc_version"; then
    echo "compiler '$CC' is GNU $couplings_cc_vendor-$couplings_cc_version"
  fi

  # Some version numbers
  couplings_cc_vers_major=`echo $couplings_cc_version | cut -f1 -d.`
  couplings_cc_vers_minor=`echo $couplings_cc_version | cut -f2 -d.`
  couplings_cc_vers_patch=`echo $couplings_cc_version | cut -f3 -d.`
  test -n "$couplings_cc_vers_major" || couplings_cc_vers_major=0
  test -n "$couplings_cc_vers_minor" || couplings_cc_vers_minor=0
  test -n "$couplings_cc_vers_patch" || couplings_cc_vers_patch=0

  # Default compiler flags
  cflags_default="-ansi -funsigned-char -pedantic -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused"
  cflags_default_dbg="-g"
  cflags_default_opt="-O2"
  cflags_default_prf="-pg"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in

    *i?86|*x86_64)
      cflags_default_opt="-funroll-loops -O2 -Wuninitialized"
      case "$host_cpu" in
        i686)
          case "$couplings_cc_vendor-$couplings_cc_version" in
            gcc-2.9[56]*|gcc-3*|gcc-4*)
              cflags_default_opt="$cflags_default_opt -march=i686"
            ;;
          esac
      esac
      ;;

    *alphaev6|*alphaev67|*alphaev68|*alphaev7)
      cflags_default_opt="-mcpu=ev6 -O"
      ;;

  esac

  # Modify default flags depending on gcc version (as older versions
  # may not handle all flags)

  case "$couplings_cc_vendor-$couplings_cc_version" in

    gcc-2.9[56]*)
      cflags_default="$cflags_default -Wno-long-long"
      ;;

    gcc-3.*|gcc-4.*)
      cflags_default="`echo $cflags_default | sed -e 's/-ansi/-std=c99/g'`"
      cflags_default="$cflags_default -Wfloat-equal -Wpadded"
      ;;

  esac

# Otherwise, are we using icc ?
#------------------------------

elif test "x$couplings_gcc" = "xicc"; then

  couplings_cc_version=`echo $CC --version | grep icc |sed 's/[a-zA-Z()]//g'`

  echo "compiler '$CC' is Intel ICC"

  # Version strings for logging purposes and known compiler flag
  $CC -V conftest.c > $outfile 2>&1
  couplings_ac_cc_version=`$CC --version 2>&1 | head -1`
  couplings_compiler_known=yes

  # Default compiler flags
  cflags_default="-strict-ansi -std=c99 -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0 -traceback"
  cflags_default_opt="-O2"
  cflags_default_prf="-p"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in
    *ia64)
      cflags_default_opt="-O2 -mcpu=itanium2-p9000"
      ;;
  esac

# Otherwise, are we using pgcc ?
#-------------------------------

else

  $CC -V 2>&1 | grep 'The Portland Group' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Portland Group pgcc"

    # Version strings for logging purposes and known compiler flag
    $CC -V conftest.c > $outfile 2>&1
    couplings_ac_cc_version=`grep -i pgcc $outfile`
    couplings_compiler_known=yes

    # Default compiler flags
    cflags_default="-Xa -fPIC"
    cflags_default_dbg="-g"
    cflags_default_opt="-fast -fastsse"
    cflags_default_prf="-Mprof=func,lines"

  fi

fi

# Compiler still not identified
#------------------------------

if test "x$couplings_compiler_known" != "xyes" ; then

  case "$host_os" in

    linux* | none)

      # IBM Blue Gene
      #--------------

      $CC -qversion 2>&1 | grep 'XL C' | grep 'Blue Gene' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is IBM XL C compiler for Blue Gene"

        # Version strings for logging purposes and known compiler flag
        $CC -qversion > $outfile 2>&1
        couplings_ac_cc_version=`grep 'XL C' $outfile`
        couplings_compiler_known=yes

        # Default compiler flags
        if test -d /bgl/BlueLight/ppcfloor/bglsys/include ; then
          cppflags_default="-I/bgl/BlueLight/ppcfloor/bglsys/include"
          cflags_default="-qmaxmem=-1 -qarch=440d -qtune=440"
        elif test -d /bgsys/drivers/ppcfloor/comm/include ; then
          cppflags_default="-I/bgsys/drivers/ppcfloor/comm/include"
          cflags_default="-qmaxmem=-1 -qarch=450d -qtune=450"
        else
          cppflags_default=""
        fi
        cflags_default_opt="-g -O3"
        cflags_default_prf="-pg"

        # Default linker flags
        if test -d /bgl/BlueLight/ppcfloor/bglsys ; then
          ldflags_default="-Wl,-allow-multiple-definition -L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts -lnss_files -lnss_dns -lresolv"
        elif test -d /bgsys/drivers/ppcfloor/comm/include ; then
          ldflags_default="-Wl,-allow-multiple-definition -L/bgsys/drivers/ppcfloor/comm/lib -lmpich.cnk -ldcmfcoll.cnk -ldcmf.cnk"
        else
          ldflags_default=""
        fi
        ldflags_default_opt="-g -O3"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-pg"
        couplings_linker_set=yes

        # Disable shared libraries in all cases
        couplings_disable_shared=yes

      fi
      ;;

    osf*)

      # Native Compaq Tru64 Unix C compiler
      #------------------------------------

      $CC -V 2>&1 | grep 'Compaq Tru64' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is Compaq Tru64 compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        couplings_ac_cc_version=`grep 'Compaq C' $outfile`
        couplings_compiler_known=yes

        # Default compiler flags
        case "$host_cpu" in
          alphaev6|alphaev67|alphaev68|alphaev7)
            cflags_default="-arch host -tune host -ansi_alias -std -check_bounds -trapuv -check -msg_enable alignment -msg_enable noansi -msg_enable performance -portable -msg_enable c_to_cxx"
            cflags_default_opt="-O"
            cflags_default_dbg="-g"
            cflags_default_prf="-pg"
          ;;
        esac

      fi
      ;;

    uxpv*)

      # Native Fujitsu vectorizing C compiler (tested on VPP5000)
      #---------------------------------------

      $CC -V 2>&1 | grep 'Fujitsu' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is Fujitsu compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        couplings_ac_cc_version=`grep ccom $outfile`
        couplings_compiler_known=yes
        couplings_linker_set=yes

        # Default compiler flags
        cflags_default="-KA64 -Kvp"
        cflags_default_opt="-O"
        cflags_default_dbg="-g -Kargchk -w4"
        cflags_default_prf="-pg"

        # Default linker flags
        ldflags_default="-Kvp -Kargchk -Wl,-S1:d"
        ldflags_default_opt="-O"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-pg"
        couplings_linker_set=yes

      fi
      ;;

    irix5.*|irix6.*)

      # Native SGI IRIX C compiler
      #---------------------------

      $CC -version 2>&1 | grep 'MIPSpro' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is MIPSpro compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -version > $outfile 2>&1
        couplings_ac_cc_version=`grep MIPSpro $outfile`
        couplings_compiler_known=yes

        # Default compiler flags
        cflags_default="-c99 -64"
        cflags_default_opt="-O2 -woff 1429,1521"
        cflags_default_dbg="-g -woff 1429,1521,1209 -fullwarn"
        cflags_default_prf="-O0"

      fi
      ;;

    hpux*)

      # Native HP-UX C compiler
      #------------------------

      $CC -V conftest.c 2>&1 | grep 'HP' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is HP compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        couplings_ac_cc_version=`grep ccom $outfile`
        couplings_compiler_known=yes
        couplings_linker_set=yes

        # Default compiler flags
        cflags_default="-Aa +e +DA2.0W"
        cflags_default_opt="+O2"
        cflags_default_dbg="-g"
        cflags_default_prf="-G"

        # Default linker flags
        ldflags_default="+DA2.0W"
        ldflags_default_opt="+O2"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-fbexe"

      fi
      ;;

    solaris2.*)

      # Sun Workshop compiler
      #----------------------

      $CC -V 2>&1 | grep 'WorkShop Compilers' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is Sun WorkShop Compilers"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        couplings_ac_cc_version=`grep cc $outfile`
        couplings_compiler_known=yes

        # Default compiler flags
        cflags_default="-Xa"
        cflags_default_opt="-xO2"
        cflags_default_dbg="-g"
        cflags_default_prf="-pg"

     fi
     ;;

    *)

      # Unknown
      #--------

      cflags_default=""
      cflags_default_opt="-O"
      cflags_default_dbg="-g"
      cflags_default_prf=""

      ;;

  esac

fi

# Default linker flags
#---------------------

if test "x$couplings_linker_set" != "xyes" ; then

  case "$host_os" in

    linux*)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ldflags_default_prf="-pg"
      ;;

    osf*)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g3"
      ldflags_default_prf="-pg"
      ;;

    irix5.*|irix6.*)
      ldflags_default="-64 -Wl,-woff,85"
      ldflags_default_opt=""
      ldflags_default_dbg="-g"
      ldflags_default_prf=""
      ;;

    solaris2.*)
      ldflags_default_opt=""
      ldflags_default_dbg="-g"
      ldflags_default_prf=""
      ;;

    *)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ldflags_default_prf="-pg"
      ;;

  esac

fi

# Finish

export LANG=$save_LANG

if test -f $outfile ; then 
  couplings_ac_cc_version_full=`sed -e '11,$d' $outfile`
fi

# Clean temporary files

rm -f conftest* a.out $outfile

