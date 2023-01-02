#!/usr/bin/env python
#-----------------------------------------------------------------------------
# This file is part of the CWIPI library.
#
# Copyright (C) 2011  ONERA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

import mpi4py.MPI as MPI
import numpy as np
import sys

f=None

def runTest():
    """
    Run tests on Python interface of new API
    """
    global f
    comm = MPI.COMM_WORLD

    i_rank = comm.rank
    n_rank = comm.size

    if (i_rank == 0):
        print("\nSTART: python_api.py")

    if (n_rank != 2):
        if i_rank == 0:
            print("      Not executed : only available for 2 processes")
        return

    code_names = ["proc0","proc1"]

    try:
        from cwipi import cwipi
    except:
        if i_rank == 0:
            print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        sys.exit(1)

    # OUTPUT
    srank = '{0}'.format(rank)
    f=open("python_api_"+srank.zfill(4)+".txt",'w')
    cwipi.output_file_set(f)

    # INIT
    n_code = 1
    out = cwipi.init(comm,
                     n_code,
                     code_names[i_rank])
    f.write("cwipi.init:\n")
    f.write("  - is_active_rank : {param}\n".format(param=is_active_rank))
    f.write("  - time_init : {param}\n".format(param=time_init))
    f.write("  - intra_comms : {param}\n".format(param=intra_comms[0]))


    # PROPERTIES DUMP
    f.flush()
    cwipi.properties_dump()

    # END
    f.write("\nEnd.\n")

if __name__ == '__main__':
    runTest()
