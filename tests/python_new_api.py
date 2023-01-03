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
        from pycwp import pycwp
    except:
        if i_rank == 0:
            print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        sys.exit(1)

    # OUTPUT
    srank = '{0}'.format(rank)
    f=open("python_api_"+srank.zfill(4)+".txt",'w')
    pycwp.output_file_set(f)

    # INIT
    n_code = 1
    out = pycwp.init(comm,
                   n_code,
                   code_names[i_rank])
    f.write("pycwp.init:\n")
    f.write("  - is_active_rank : {param}\n".format(param=out["is_active_rank"]))
    f.write("  - time_init : {param}\n".format(param=out["time_init"]))
    f.write("  - intra_comms : {param}\n".format(param=out["intra_comms"][0]))

    # STATE UPDATE
    pycwp.state_update(code_names[i_rank], pycwp.CWP_STATE_END)
    state = pycwp.state_get(code_names[i_rank])
    f.write("pycwp.state_get:\n")
    f.write("  - state : {param}\n".format(param=state))
    pycwp.state_update(code_names[i_rank], pycwp.CWP_STATE_IN_PROGRESS)

    # TIME UPDATE
    pycwp.time_update(code_names[i_rank], 0.0)

    # PROPERTIES DUMP
    f.flush()
    pycwp.properties_dump()

    # CODES
    n_code = pycwp.codes_nb_get()
    code   = pycwp.codes_list_get()
    n_loc_code = pycwp.loc_codes_nb_get()
    loc_code   = pycwp.loc_codes_list_get()
    f.write("pycwp.code:\n")
    f.write("  - n_code : {param}\n".format(param=n_code))
    for i in range(n_code):
        f.write("    --> {param}\n".format(param=code[i]))
    f.write("  - n_loc_code : {param}\n".format(param=n_loc_code))
    for i in range(n_loc_code):
        f.write("    --> {param}\n".format(param=loc_code[i]))

    # PARAM
    pycwp.param_lock(code_names[i_rank])
    pycwp.param_add_dbl(code_names[i_rank], "double", 0.5)
    pycwp.param_unlock(code_names[i_rank])

    pycwp.param_lock(code_names[i_rank])
    pycwp.param_add_str(code_names[i_rank], "str", "chat")
    pycwp.param_unlock(code_names[i_rank])

    if (i_rank == 0):
        pycwp.param_lock(code_names[i_rank])
        pycwp.param_add_int(code_names[i_rank], "entier", 1)
        pycwp.param_unlock(code_names[i_rank])

    if (i_rank == 1):
        pycwp.param_lock(code_names[i_rank])
        pycwp.param_add_int(code_names[i_rank], "entier", -1)
        pycwp.param_unlock(code_names[i_rank])

    value = pycwp.param_get(code_names[i_rank], "double", pycwp.CWP_DOUBLE)
    f.write("cwp.param_get ({param}):\n".format(param=i_rank))
    f.write("  - value (0): {param}\n".format(param=value))

    pycwp.param_lock(code_names[i_rank])
    pycwp.param_set_dbl(code_names[i_rank], "double", 0.25)
    pycwp.param_unlock(code_names[i_rank])

    pycwp.param_lock(code_names[i_rank])
    pycwp.param_set_str(code_names[i_rank], "str", "chien")
    pycwp.param_unlock(code_names[i_rank])

    if (i_rank == 0):
        pycwp.param_lock(code_names[i_rank])
        pycwp.param_set_int(code_names[i_rank], "entier", 2)
        pycwp.param_unlock(code_names[i_rank])

    value = pycwp.param_get(code_names[i_rank], "double", pycwp.CWP_DOUBLE)
    f.write("pycwp.param_get ({param}):\n".format(param=i_rank))
    f.write("  - value (1): {param}\n".format(param=value))

    pycwp.param_lock(code_names[i_rank])
    pycwp.param_del(code_names[i_rank], "str", pycwp.CWP_CHAR)
    pycwp.param_unlock(code_names[i_rank])

    n_param_str = pycwp.param_n_get(code_names[i_rank], pycwp.CWP_CHAR)
    n_param_int = pycwp.param_n_get(code_names[i_rank], pycwp.CWP_INT)
    f.write("pycwp.param_n_get:\n")
    f.write("  - n_param_str: {param}\n".format(param=n_param_str))
    f.write("  - n_param_int: {param}\n".format(param=n_param_int))

    double_param = pycwp.param_list_get(code_names[i_rank], pycwp.CWP_DOUBLE)
    f.write("pycwp.param_list_get:\n")
    f.write("  - double_param: {param}\n".format(param=double_param[0]))

    bool_int = pycwp.param_is(code_names[i_rank], "entier", pycwp.CWP_INT)
    bool_int = pycwp.param_is(code_names[i_rank], "chapeau", pycwp.CWP_INT)

    f.write("pycwp.param_is:\n")
    f.write("  - bool_int 'entier': {param}\n".format(param=bool_int))
    f.write("  - bool_int 'chapeau': {param}\n".format(param=bool_int))

    result = pycwp.param_reduce(pycwp.CWP_OP_MIN, "entier",  pycwp.CWP_INT, 2, code_names)
    f.write("pycwp.param_reduce:\n")
    f.write("  - result: {param}\n".format(param=result))

    # Cpl
    cpl = pycwp.Cpl(code_names[i_rank],
                    "test",
                    code_names[(i_rank+1)%2],
                    pycwp.CWP_INTERFACE_LINEAR,
                    pycwp.CWP_COMM_PAR_WITH_PART,
                    pycwp.CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                    1,
                    pycwp.CWP_DYNAMIC_MESH_STATIC,
                    pycwp.CWP_TIME_EXCH_USER_CONTROLLED)

    # MESH

    if (i_rank == 0):
        coord = np.array([0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0], dtype=np.double)
        connec_idx = np.array([0, 3, 6], dtype=np.int32)
        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

    if (i_rank == 1):
        coord = np.array([0, 1, 0, 0, 2, 0, 1, 1, 0, 1, 2, 0], dtype=np.double)
        connec_idx = np.array([0, 3, 6], dtype=np.int32)
        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

    cpl.mesh_interf_vtx_set(0,
                            4,
                            coord,
                            NULL)

    cpl.mesh_interf_finalize()

    # cpl.mesh_interf_block_add()
    # cpl.mesh_interf_block_std_set()
    # cpl.mesh_interf_block_std_get()

    cpl.mesh_interf_f_poly_block_set(0,
                                     0,
                                     2,
                                     connec_idx,
                                     connec,
                                     NULL)

    out = cpl.mesh_interf_f_poly_block_get(0, 0)
    f.write("mesh_interf_f_poly_block_get ({param}):\n".format(param=i_rank))
    f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
    f.write("  - connec_idx[0] : {param}\n".format(param=out["connec_idx"][0]))
    f.write("  - connec[0] : {param}\n".format(param=out["connec"][0]))
    f.write("  - global_num : {param}\n".format(param=out["global_num"]))

    # cpl.mesh_interf_c_poly_block_set()
    # cpl.mesh_interf_c_poly_block_get()
    # cpl.mesh_interf_from_cellface_set()
    # cpl.mesh_interf_from_faceedge_set()

    # FIELD to do

    sendField=np.array([0.0, 0.1, 0.2, 0.3], dtype=np.double)
    recvField=np.arange(4, dtype=np.double)

    cpl.field_create("champs",
                     pycwp.CWP_DOUBLE,
                     pycwp.CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     pycwp.CWP_DOF_LOCATION_NODE,
                     pycwp.CWP_FIELD_EXCH_SENDRECV,
                     pycwp.CWP_STATUS_OFF)

    cpl.field_set("champs",
                  0,
                  pycwp.CWP_FIELD_MAP_SOURCE,
                  sendField)

    out = cpl.field_get("champs")
    f.write("cpl.field_get ({param}):\n".format(param=i_rank))
    f.write("  - n_comp : {param}\n".format(param=out["n_comp"]))
    f.write("  - dof_loc : {param}\n".format(param=out["dof_loc"]))
    f.write("  - storage : {param}\n".format(param=out["storage"]))

    cpl.field_del("champs")

    # SPATIAL INTERPOLATION to do
    # cpl.spatial_interp_weights_compute()
    # cpl.spatial_interp_property_set()

    # VISU
    cpl.visu_set(1,
                 pycwp.CWP_VISU_FORMAT_ENSIGHT,
                 "text")

    # USER TGT PTS to do
    # cpl.user_tgt_pts_set()

    # USER INTERPOLATION to do
    # cpl.interp_from_location_unset()
    # cpl.interp_from_location_set()

    # SEND/RECV to do
    # pycwp.field_issend()
    # pycwp.field_irecv()
    # pycwp.field_wait_issend()
    # pycwp.field_wait_irecv()

    # USER STRUCTURE to do
    # pycwp.user_structure_set()
    # pycwp.user_structure_get()

    cpl.mesh_interf_del()
    del cpl

    # FINALIZE
    pycwp.finalize()

    # END
    f.write("\nEnd.\n")

if __name__ == '__main__':
    runTest()
