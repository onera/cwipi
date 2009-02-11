#!/bin/sh

MPI_RANK=`mpi_rank.sh`

if [ $MPI_RANK -ne 4 ] ; then 
  test1_f
else
  test1_c
fi
