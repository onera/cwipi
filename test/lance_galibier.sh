#!/bin/sh

MPI_RANK=`mpi_rank.sh`

if [ $MPI_RANK -eq 0 ] ; then 
  test1_f
else
  test1_c
fi
