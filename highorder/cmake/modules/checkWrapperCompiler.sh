#!/bin/bash
# Check if compiler of mpi-wrapper is the same as classic compiler
# CMake standard : 0 for false, 1 for true
# Check intel, gnu, pgi, ibm compilers
# $1 defines kind of compilers
# $2 defines the wrapper
# $3 defines the reference compiler

if [ $1 = 'CC' ]
    then
    if   [ -n "$($2 --version 2>&1 | grep -i "^icc")" ];        then    CC=icc
    elif [ -n "$($2 --version 2>&1 | grep -i "^gcc")" ];        then	CC=gcc
    elif [ -n "$($2 -qversion 2>&1 | grep "XL C")" ];           then	CC=xlc
    elif [ -n "$($2 -V 2>&1 | grep "The Portland Group"  )" ];  then	CC=pgcc
    else exit 0
    fi
    
    if [ $CC = $(basename $3) ]; then exit 1; else exit 0; fi
fi


if [ $1 = 'CXX' ]
    then
    if   [ -n "$($2 --version 2>&1 | grep -i "^icpc")" ];       then    CXX=icpc
    elif [ -n "$($2 --version 2>&1 | grep -i "^g++")" ];        then	CXX=g++
    elif [ -n "$($2 -qversion 2>&1 | grep "XL C")" ];           then	CXX=xlc
    elif [ -n "$($2 -V 2>&1 | grep "The Portland Group"  )" ];  then	CXX=pgc++
    else exit 0
    fi

    if [ $CXX = $(basename $3) ]; then exit 1; else exit 0; fi
fi

if [ $1 = 'FC' ]
    then
    if   [ -n "$($2 --version 2>&1 | grep -i "^ifort")" ];       then   FC=ifort
    elif [ -n "$($2 --version 2>&1 | grep -i "^GNU Fortran")" ]; then	FC=gfortran
    elif [ -n "$($2 -qversion 2>&1 | grep "XL Fortran")" ];      then	FC=xlf
    elif [ -n "$($2 -V 2>&1 | grep "The Portland Group"  )" ];   then	FC=pgfortran
    else exit 0
    fi
    
    if [ $FC = $(basename $3) ]; then exit 1; else exit 0; fi
fi
