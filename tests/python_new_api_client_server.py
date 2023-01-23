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

def runTest():
    """
    Run tests on Python interface of new API
    """
    global f
    comm = MPI.COMM_WORLD

    i_rank = comm.rank
    n_rank = comm.size

    if (i_rank == 0):
        print("\nSTART: python_new_api_client_server.py")

    if (n_rank != 2):
        if i_rank == 0:
            print("      Not executed : only available for 2 processes")
        return

    code_names = ["proc0","proc1"]

    if (i_rank == 0):
        code_name = "proc0"

    if (i_rank == 1):
        code_name = "proc1"

    try:
        from pycwpclt import pycwpclt
    except:
        if i_rank == 0:
            print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        sys.exit(1)

    # INIT
    print("pycwpclt.init:\n")
    config = "cwp_config_srv.txt"
    out = pycwpclt.init(comm,
                        config,
                        code_name)
    print("  - is_active_rank : {param}\n".format(param=out["is_active_rank"]))
    print("  - time_init : {param}\n".format(param=out["time_init"]))

    # TO DO

    # FINALIZE
    pycwpclt.finalize()

    # END
    print("\nEnd.\n")
    comm.Barrier()
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
