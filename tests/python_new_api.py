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

    # STATE UPDATE
    cwipi.state_update(code_names[i_rank], CWP_STATE_END)
    state = cwipi.state_get(code_names[i_rank])
    f.write("cwipi.state_get:\n")
    f.write("  - state : {param}\n".format(param=state))
    cwipi.state_update(code_names[i_rank], CWP_STATE_IN_PROGRESS)

    # TIME UPDATE
    cwipi.time_update(code_names[i_rank], 0.0)

    # PROPERTIES DUMP
    f.flush()
    cwipi.properties_dump()

    # CODES
    n_code = cwipi.codes_nb_get()
    code   = cwipi.codes_list_get()
    n_loc_code = cwipi.loc_codes_nb_get()
    loc_code   = cwipi.loc_codes_list_get()
    f.write("cwipi.code:\n")
    f.write("  - n_code : {param}\n".format(param=n_code))
    for i in range(n_code):
        f.write("    --> {param}\n".format(param=code[i]))
    f.write("  - n_loc_code : {param}\n".format(param=n_loc_code))
    for i in range(n_loc_code):
        f.write("    --> {param}\n".format(param=loc_code[i]))

    # PARAM
    cwipi.param_lock(code_names[i_rank])
    cwipi.param_add_dbl(code_names[i_rank], "double", 0.5)
    cwipi.param_unlock(code_names[i_rank])

    cwipi.param_lock(code_names[i_rank])
    cwipi.param_add_str(code_names[i_rank], "str", "chat")
    cwipi.param_unlock(code_names[i_rank])

    if (i_rank == 0):
        cwipi.param_lock(code_names[i_rank])
        cwipi.param_add_int(code_names[i_rank], "entier", 1)
        cwipi.param_unlock(code_names[i_rank])

    if (i_rank == 1):
        cwipi.param_lock(code_names[i_rank])
        cwipi.param_add_int(code_names[i_rank], "entier", -1)
        cwipi.param_unlock(code_names[i_rank])

    value = cwipi.param_get(code_names[i_rank], "double", CWP_DOUBLE)
    f.write("cwipi.param_get ({param}):\n".format(param=i_rank))
    f.write("  - value (0): {param}\n".format(param=value))

    cwipi.param_lock(code_names[i_rank])
    cwipi.param_set_dbl(code_names[i_rank], "double", 0.25)
    cwipi.param_unlock(code_names[i_rank])

    cwipi.param_lock(code_names[i_rank])
    cwipi.param_set_str(code_names[i_rank], "str", "chien")
    cwipi.param_unlock(code_names[i_rank])

    if (i_rank == 0):
        cwipi.param_lock(code_names[i_rank])
        cwipi.param_set_int(code_names[i_rank], "entier", 2)
        cwipi.param_unlock(code_names[i_rank])

    value = cwipi.param_get(code_names[i_rank], "double", CWP_DOUBLE)
    f.write("cwipi.param_get ({param}):\n".format(param=i_rank))
    f.write("  - value (1): {param}\n".format(param=value))

    cwipi.param_lock(code_names[i_rank])
    cwipi.param_del(code_names[i_rank], "str", CWP_CHAR)
    cwipi.param_unlock(code_names[i_rank])

    n_param_str = cwipi.param_n_get(code_names[i_rank], CWP_CHAR)
    n_param_int = cwipi.param_n_get(code_names[i_rank], CWP_INT)
    f.write("cwipi.param_n_get:\n")
    f.write("  - n_param_str: {param}\n".format(param=n_param_str))
    f.write("  - n_param_int: {param}\n".format(param=n_param_int))

    double_param = cwipi.param_list_get(code_names[i_rank], CWP_DOUBLE)
    f.write("cwipi.param_list_get:\n")
    f.write("  - double_param: {param}\n".format(param=double_param[0]))

    bool_int = cwipi.param_is(code_names[i_rank], "entier", CWP_INT)
    bool_int = cwipi.param_is(code_names[i_rank], "chapeau", CWP_INT)

    f.write("cwipi.param_is:\n")
    f.write("  - bool_int 'entier': {param}\n".format(param=bool_int))
    f.write("  - bool_int 'chapeau': {param}\n".format(param=bool_int))

    result = cwipi.param_reduce(CWP_OP_MIN, "entier",  CWP_INT, 2, code_names)
    f.write("cwipi.param_reduce:\n")
    f.write("  - result: {param}\n".format(param=result))

    # Cpl
    cpl = cwipi.Cpl(code_names[i_rank],
                    "test",
                    code_names[(i_rank+1)%2],
                    CWP_INTERFACE_VOLUME,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                    1,
                    CWP_DYNAMIC_MESH_STATIC,
                    CWP_TIME_EXCH_USER_CONTROLLED)

    # FIELD to do
    # cpl.field_create()
    # cpl.field_set()
    # cpl.field_get()
    # cpl.field_del()

    # MESH to do
    # cpl.mesh_interf_finalize()
    # cpl.mesh_interf_vtx_set()
    # cpl.mesh_interf_block_add()
    # cpl.mesh_interf_block_std_set()
    # cpl.mesh_interf_block_std_get()
    # cpl.mesh_interf_f_poly_block_set()
    # cpl.mesh_interf_f_poly_block_set()
    # cpl.mesh_interf_c_poly_block_set()
    # cpl.mesh_interf_c_poly_block_get()
    # cpl.mesh_interf_del()
    # cpl.mesh_interf_from_cellface_set()
    # cpl.mesh_interf_from_faceedge_set()

    # SPATIAL INTERPOLATION to do
    # cpl.spatial_interp_weights_compute()
    # cpl.spatial_interp_property_set()

    # VISU
    cpl.visu_set(1,
                 CWP_VISU_FORMAT_ENSIGHT,
                 "text")

    # USER TGT PTS to do
    # cpl.user_tgt_pts_set()

    # USER INTERPOLATION to do
    # cpl.interp_from_location_unset()
    # cpl.interp_from_location_set()

    # SEND/RECV to do
    # cwipi.field_issend()
    # cwipi.field_irecv()
    # cwipi.field_wait_issend()
    # cwipi.field_wait_irecv()

    # USER STRUCTURE to do
    # cwipi.user_structure_set()
    # cwipi.user_structure_get()

    del cpl

    # FINALIZE
    cwipi.finalize()

    # END
    f.write("\nEnd.\n")

if __name__ == '__main__':
    runTest()
