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
    config = "../server/cwp_config_srv.txt"
    is_active_rank = np.array([1], dtype=np.int32)
    time_init = np.array([0.], dtype=np.double)
    out = pycwpclt.init(comm,
                        config,
                        code_name,
                        is_active_rank,
                        time_init)

    # STATE UPDATE
    pycwpclt.state_update(code_names[i_rank], pycwpclt.STATE_IN_PROGRESS)
    print("pycwpclt.state_get:\n")
    state = pycwpclt.state_get(code_names[i_rank])
    print("  - state : {param}\n".format(param=state))
    pycwpclt.state_update(code_names[i_rank], pycwpclt.STATE_END)

    # PROPERTIES DUMP
    print("pycwpclt.properties_dump:\n")
    pycwpclt.properties_dump()

    # CODES
    print("pycwpclt.code:\n")
    n_code = pycwpclt.codes_nb_get()
    code   = pycwpclt.codes_list_get()
    n_loc_code = pycwpclt.loc_codes_nb_get()
    loc_code   = pycwpclt.loc_codes_list_get()
    print("  - n_code : {param}\n".format(param=n_code))
    for i in range(n_code):
        print("    --> {param}\n".format(param=code[i]))
    print("  - n_loc_code : {param}\n".format(param=n_loc_code))
    for i in range(n_loc_code):
        print("    --> {param}\n".format(param=loc_code[i]))

    # PARAM
    pycwpclt.param_lock(code_names[i_rank])
    pycwpclt.param_add_dbl(code_names[i_rank], "double", 0.5)
    pycwpclt.param_unlock(code_names[i_rank])

    pycwpclt.param_lock(code_names[i_rank])
    pycwpclt.param_add_str(code_names[i_rank], "str", "chat")
    pycwpclt.param_unlock(code_names[i_rank])

    if (i_rank == 0):
        pycwpclt.param_lock(code_names[i_rank])
        pycwpclt.param_add_int(code_names[i_rank], "entier", 1)
        pycwpclt.param_unlock(code_names[i_rank])

    if (i_rank == 1):
        pycwpclt.param_lock(code_names[i_rank])
        pycwpclt.param_add_int(code_names[i_rank], "entier", -1)
        pycwpclt.param_unlock(code_names[i_rank])

    comm.Barrier()

    print("cwp.param_get ({param}):\n".format(param=i_rank))
    value = pycwpclt.param_get(code_names[i_rank], "double", pycwpclt.DOUBLE)
    print("  - value (0): {param}\n".format(param=value))

    pycwpclt.param_lock(code_names[i_rank])
    pycwpclt.param_set_dbl(code_names[i_rank], "double", 0.25)
    pycwpclt.param_unlock(code_names[i_rank])

    pycwpclt.param_lock(code_names[i_rank])
    pycwpclt.param_set_str(code_names[i_rank], "str", "chien")
    pycwpclt.param_unlock(code_names[i_rank])

    if (i_rank == 0):
        pycwpclt.param_lock(code_names[i_rank])
        pycwpclt.param_set_int(code_names[i_rank], "entier", 2)
        pycwpclt.param_unlock(code_names[i_rank])

    comm.Barrier()

    print("pycwpclt.param_get ({param}):\n".format(param=i_rank))
    value = pycwpclt.param_get(code_names[i_rank], "double", pycwpclt.DOUBLE)
    print("  - value (1): {param}\n".format(param=value))

    pycwpclt.param_lock(code_names[i_rank])
    pycwpclt.param_del(code_names[i_rank], "str", pycwpclt.CHAR)
    pycwpclt.param_unlock(code_names[i_rank])

    comm.Barrier()

    print("pycwpclt.param_n_get:\n")
    n_param_str = pycwpclt.param_n_get(code_names[i_rank], pycwpclt.CHAR)
    n_param_int = pycwpclt.param_n_get(code_names[i_rank], pycwpclt.INT)
    print("  - n_param_str: {param}\n".format(param=n_param_str))
    print("  - n_param_int: {param}\n".format(param=n_param_int))

    print("pycwpclt.param_list_get:\n")
    str_param = pycwpclt.param_list_get(code_names[i_rank], pycwpclt.CHAR)
    for i in range(str_param['n_param']):
        print("    --> str_param: {param}\n".format(param=str_param['param_names'][i]))

    print("pycwpclt.param_is:\n")
    bool_int = pycwpclt.param_is(code_names[i_rank], "entier", pycwpclt.INT)
    print("  - bool_int 'entier': {param}\n".format(param=bool_int))
    bool_int = pycwpclt.param_is(code_names[i_rank], "chapeau", pycwpclt.INT)
    print("  - bool_int 'chapeau': {param}\n".format(param=bool_int))

    comm.Barrier()

    print("pycwpclt.param_list_get:\n")
    int_param = pycwpclt.param_list_get(code_names[i_rank], pycwpclt.INT)
    for i in range(int_param['n_param']):
        print("    --> int_param: {param}\n".format(param=int_param['param_names'][i]))

    print("pycwpclt.param_get ({param}):\n".format(param=i_rank))
    value = pycwpclt.param_get(code_names[i_rank], "entier", pycwpclt.INT)
    print("  - value int: {param}\n".format(param=value))

    comm.Barrier()

    print("pycwpclt.param_reduce:\n")
    result = pycwpclt.param_reduce(pycwpclt.OP_MIN, "entier",  pycwpclt.INT, 2, code_names)
    print("  - result: {param}\n".format(param=result))

    # Cpl
    print("pycwpclt.Coupling:\n")
    cpl = pycwpclt.Coupling(code_names[i_rank],
                            "test",
                            code_names[(i_rank+1)%2],
                            pycwpclt.INTERFACE_SURFACE,
                            pycwpclt.COMM_PAR_WITH_PART,
                            pycwpclt.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                            1,
                            pycwpclt.DYNAMIC_MESH_VARIABLE,
                            pycwpclt.TIME_EXCH_USER_CONTROLLED)

    # VISU
    cpl.visu_set(1,
                 pycwpclt.VISU_FORMAT_ENSIGHT,
                 "text")

    if (i_rank == 0):
        coord = np.array([0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0], dtype=np.double)
        connec_idx = np.array([0, 3, 6], dtype=np.int32)
        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

    if (i_rank == 1):
        coord = np.array([0, 1, 0, 0, 2, 0, 1, 1, 0, 1, 2, 0], dtype=np.double)
        connec_idx = np.array([0, 3, 6], dtype=np.int32)
        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

    print("cpl.mesh_interf_vtx_set:\n")
    cpl.mesh_interf_vtx_set(0,
                            4,
                            coord,
                            None)

    comm.Barrier()

    print("cpl.mesh_interf_block_add:\n")
    block_id = cpl.mesh_interf_block_add(pycwpclt.BLOCK_FACE_POLY)

    print("cpl.mesh_interf_f_poly_block_set ({param}):\n".format(param=i_rank))
    cpl.mesh_interf_f_poly_block_set(0,
                                     block_id,
                                     2,
                                     connec_idx,
                                     connec,
                                     None)

    print("cpl.mesh_interf_finalize:\n")
    cpl.mesh_interf_finalize()

    print("cpl.mesh_interf_f_poly_block_get:\n")
    out = cpl.mesh_interf_f_poly_block_get(0, block_id)
    print("  - n_elts : {param}\n".format(param=out["n_elts"]))
    print("  - connec_idx {param}\n".format(param=out["connec_idx"]))
    print("  - connec {param}\n".format(param=out["connec"]))
    print("  - global_num : {param}\n".format(param=out["global_num"]))

    # FINALIZE
    pycwpclt.finalize()

    # END
    print("\nEnd.\n")
    comm.Barrier()
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
