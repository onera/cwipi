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
import os
import time

def runTest():
    """
    Run tests on Python interface of new API
    """
    global f
    comm = MPI.COMM_WORLD

    i_rank = comm.rank
    n_rank = comm.size

    # launch server
    if (i_rank == 0):
        config = "python_new_api_client_server_o/code1/cwp_config_srv.txt"
    else:
        config = "python_new_api_client_server_o/code2/cwp_config_srv.txt"

    if (i_rank == 0):
        os.system("mkdir -p python_new_api_client_server_o/code1")
        os.system("mkdir -p python_new_api_client_server_o/code2")
        os.system("rm -f ./python_new_api_client_server_o/code1/cwp_config_srv.txt")
        os.system("rm -f ./python_new_api_client_server_o/code2/cwp_config_srv.txt")
        os.system("mpiexec -n 1 ../bin/cwp_server -cn code0 -p 49100 49100 -c \"python_new_api_client_server_o/code1/cwp_config_srv.txt\" : -n 1  ../bin/cwp_server -cn code1 -p 49101 49101 -c \"python_new_api_client_server_o/code2/cwp_config_srv.txt\" &")

    while (os.access(config, os.R_OK) != 0):
        time.sleep(1)

    time.sleep(5)

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

    intra_comm = comm.Split(i_rank)

    is_active_rank = 1
    time_init = 0.
    out = pycwpclt.init(intra_comm,
                        config,
                        code_name,
                        is_active_rank,
                        time_init)

    # # STATE UPDATE
    # pycwpclt.state_update(code_names[i_rank], pycwpclt.STATE_IN_PROGRESS)
    # print("pycwpclt.state_get:\n")
    # state = pycwpclt.state_get(code_names[i_rank])
    # print("  - state : {param}\n".format(param=state))
    # pycwpclt.state_update(code_names[i_rank], pycwpclt.STATE_END)

    # # PROPERTIES DUMP
    # print("pycwpclt.properties_dump:\n")
    # pycwpclt.properties_dump()

    # # CODES
    # print("pycwpclt.code:\n")
    # n_code = pycwpclt.codes_nb_get()
    # code   = pycwpclt.codes_list_get()
    # n_loc_code = pycwpclt.loc_codes_nb_get()
    # loc_code   = pycwpclt.loc_codes_list_get()
    # print("  - n_code : {param}\n".format(param=n_code))
    # for i in range(n_code):
    #     print("    --> {param}\n".format(param=code[i]))
    # print("  - n_loc_code : {param}\n".format(param=n_loc_code))
    # for i in range(n_loc_code):
    #     print("    --> {param}\n".format(param=loc_code[i]))

    # comm.Barrier()

    # # PARAM
    # pycwpclt.param_lock(code_names[i_rank])
    # pycwpclt.param_add_dbl(code_names[i_rank], "double", 0.5)
    # pycwpclt.param_unlock(code_names[i_rank])

    # pycwpclt.param_lock(code_names[i_rank])
    # pycwpclt.param_add_str(code_names[i_rank], "str", "chat")
    # pycwpclt.param_unlock(code_names[i_rank])

    # if (i_rank == 0):
    #     pycwpclt.param_lock(code_names[i_rank])
    #     pycwpclt.param_add_int(code_names[i_rank], "entier", 1)
    #     pycwpclt.param_unlock(code_names[i_rank])

    # if (i_rank == 1):
    #     pycwpclt.param_lock(code_names[i_rank])
    #     pycwpclt.param_add_int(code_names[i_rank], "entier", -1)
    #     pycwpclt.param_unlock(code_names[i_rank])

    # comm.Barrier()

    # print("cwp.param_get ({param}):\n".format(param=i_rank))
    # value = pycwpclt.param_get(code_names[i_rank], "double", pycwpclt.DOUBLE)
    # print("  - value (0): {param}\n".format(param=value))

    # pycwpclt.param_lock(code_names[i_rank])
    # pycwpclt.param_set_dbl(code_names[i_rank], "double", 0.25)
    # pycwpclt.param_unlock(code_names[i_rank])

    # pycwpclt.param_lock(code_names[i_rank])
    # pycwpclt.param_set_str(code_names[i_rank], "str", "chien")
    # pycwpclt.param_unlock(code_names[i_rank])

    # if (i_rank == 0):
    #     pycwpclt.param_lock(code_names[i_rank])
    #     pycwpclt.param_set_int(code_names[i_rank], "entier", 2)
    #     pycwpclt.param_unlock(code_names[i_rank])

    # comm.Barrier()

    # print("pycwpclt.param_get ({param}):\n".format(param=i_rank))
    # value = pycwpclt.param_get(code_names[i_rank], "double", pycwpclt.DOUBLE)
    # print("  - value (1): {param}\n".format(param=value))

    # pycwpclt.param_lock(code_names[i_rank])
    # pycwpclt.param_del(code_names[i_rank], "str", pycwpclt.CHAR)
    # pycwpclt.param_unlock(code_names[i_rank])

    # comm.Barrier()

    # print("pycwpclt.param_n_get:\n")
    # n_param_str = pycwpclt.param_n_get(code_names[i_rank], pycwpclt.CHAR)
    # n_param_int = pycwpclt.param_n_get(code_names[i_rank], pycwpclt.INT)
    # print("  - n_param_str: {param}\n".format(param=n_param_str))
    # print("  - n_param_int: {param}\n".format(param=n_param_int))

    # print("pycwpclt.param_list_get:\n")
    # str_param = pycwpclt.param_list_get(code_names[i_rank], pycwpclt.CHAR)
    # for i in range(str_param['n_param']):
    #     print("    --> str_param: {param}\n".format(param=str_param['param_names'][i]))

    # print("pycwpclt.param_is:\n")
    # bool_int = pycwpclt.param_is(code_names[i_rank], "entier", pycwpclt.INT)
    # print("  - bool_int 'entier': {param}\n".format(param=bool_int))
    # bool_int = pycwpclt.param_is(code_names[i_rank], "chapeau", pycwpclt.INT)
    # print("  - bool_int 'chapeau': {param}\n".format(param=bool_int))

    # comm.Barrier()

    # print("pycwpclt.param_list_get:\n")
    # int_param = pycwpclt.param_list_get(code_names[i_rank], pycwpclt.INT)
    # for i in range(int_param['n_param']):
    #     print("    --> int_param: {param}\n".format(param=int_param['param_names'][i]))

    # print("pycwpclt.param_get ({param}):\n".format(param=i_rank))
    # value = pycwpclt.param_get(code_names[i_rank], "entier", pycwpclt.INT)
    # print("  - value int: {param}\n".format(param=value))

    # comm.Barrier()

    # print("pycwpclt.param_reduce:\n")
    # result = pycwpclt.param_reduce(pycwpclt.OP_MIN, "entier",  pycwpclt.INT, 2, code_names)
    # print("  - result: {param}\n".format(param=result))

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

    # global data
    global_data_name = "poffertjes"

    send_data = np.array([13, 42, 99, 68], dtype=np.int32)
    if (i_rank == 0):
        cpl.global_data_issend(global_data_name,
                               2,
                               2,
                               send_data)

    recv_data = recv_global_data = np.zeros(4, dtype=np.int32)
    if (i_rank == 1):
        cpl.global_data_irecv(global_data_name,
                              2,
                              2,
                              recv_data)

    comm.Barrier()

    if (i_rank == 0):
        cpl.global_data_wait_issend(global_data_name)

    if (i_rank == 1):
        cpl.global_data_wait_irecv(global_data_name)


    if (i_rank == 0):
        print("send_global_data : {param}\n".format(param=send_data))
    if (i_rank == 1):
        print("recv_global_data : {param}\n".format(param=recv_data))

    # part data
    part_data_name = "ratchet and clank"
    gnum_elt = [np.array([1, 2, 3], dtype=np.int64)] # implicit that n_part = 1

    if (i_rank == 0):
        cpl.part_data_create(part_data_name,
                             pycwpclt.PARTDATA_SEND,
                             gnum_elt)

    if (i_rank == 1):
        cpl.part_data_create(part_data_name,
                             pycwpclt.PARTDATA_RECV,
                             gnum_elt)

    comm.Barrier()

    send_data = [np.array([10, 11, 20, 21, 30, 31], dtype=np.int32)]
    send_request = np.zeros(1, dtype=np.int32)
    if (i_rank == 0):
        cpl.part_data_issend(part_data_name,
                             2,
                             send_data,
                             send_request)

    recv_data = [np.zeros(6, dtype=np.int32)]
    recv_request = np.zeros(1, dtype=np.int32)
    if (i_rank == 1):
        cpl.part_data_irecv(part_data_name,
                            2,
                            recv_data,
                            recv_request)

    comm.Barrier()

    if (i_rank == 0):
        cpl.part_data_wait_issend(part_data_name,
                                  send_request[0])

    if (i_rank == 1):
        cpl.part_data_wait_irecv(part_data_name,
                                 recv_request[0])

    if (i_rank == 0):
        print("send_part_data : {param}\n".format(param=send_data))
    if (i_rank == 1):
        print("recv_part_data : {param}\n".format(param=recv_data))

    # MESH
    polygon = 0
    ho = 1

    # high order mesh
    if (ho):
        block_id = cpl.mesh_interf_block_add(pycwpclt.BLOCK_FACE_TRIAHO)
        face_vtx = np.array([1, 3, 6, 2, 5, 4], dtype=np.int32)
        cpl.mesh_interf_block_ho_set(0,
                                     block_id,
                                     1,
                                     2,
                                     face_vtx,
                                     None)

        vtx_coord = np.zeros(6*3, dtype=np.double)
        ijk_grid  = np.array([0, 0,
                             2, 0,
                             0, 2,
                             1, 0,
                             1, 1,
                             0, 1], dtype=np.int32)
        k = 0
        for j in range(3):
            for i in range(3):
                vtx_coord[k] = i
                k += 1
                vtx_coord[k] = j
                k += 1
        cpl.mesh_interf_ho_ordering_from_IJK_set(pycwpclt.BLOCK_FACE_TRIAHO,
                                                 2,
                                                 ijk_grid)

        out = cpl.mesh_interf_block_ho_get(0,
                                           block_id)

        print("mesh_interf_block_ho_get:\n")
        print("  - n_elts : {param}\n".format(param=out["n_elts"]))
        print("  - order {param}\n".format(param=out["order"]))
        print("  - connec {param}\n".format(param=out["connec"]))

    # std or polygon
    elif (polygon and not ho):
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


        # FIELD

        # USER INTERPOLATION to do has to be done
        # print("cpl.interp_from_location_set:\n")
        #
        # cpl.interp_from_location_set("champs", userInterp)

        sendField = np.array([0.0, 0.1, 0.2, 0.3], dtype=np.double)
        recvField = np.arange(4, dtype=np.double)

        if (i_rank == 0):
            cpl.field_create("champs",
                             pycwpclt.DOUBLE,
                             pycwpclt.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwpclt.DOF_LOCATION_NODE,
                             pycwpclt.FIELD_EXCH_SEND,
                             pycwpclt.STATUS_OFF)
            cpl.field_set("champs",
                          0,
                          pycwpclt.FIELD_MAP_SOURCE,
                          4,
                          sendField)

        if (i_rank == 1):
            cpl.field_create("champs",
                             pycwpclt.DOUBLE,
                             pycwpclt.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwpclt.DOF_LOCATION_NODE,
                             pycwpclt.FIELD_EXCH_RECV,
                             pycwpclt.STATUS_OFF)
            cpl.field_set("champs",
                          0,
                          pycwpclt.FIELD_MAP_TARGET,
                          4,
                          recvField)

        comm.Barrier()

        print("cpl.spatial_interp_property_set:\n")

        cpl.spatial_interp_property_set("tolerance", "double", "1e-2")

        comm.Barrier()

        print("cpl.spatial_interp_weights_compute:\n")

        cpl.spatial_interp_weights_compute()

        if (i_rank == 0):
            print("pycwpclt.field_issend (0):\n")

            cpl.field_issend("champs")

        if (i_rank == 1):
            print("cpl.field_irecv (1):\n")

            cpl.field_irecv("champs")

        if (i_rank == 0):
            print("cpl.field_wait_issend (0):\n")

            cpl.field_wait_issend("champs")

        if (i_rank == 1):
            print("cpl.field_wait_irecv (1):\n")

            cpl.field_wait_irecv("champs")

        comm.Barrier()

        print("cpl.field_del:\n")

        cpl.field_del("champs")

        # USER INTERPOLATION to do has to be done
        # print("cpl.interp_from_location_unset:\n")
        #
        # cpl.interp_from_location_unset("champs")

        # TIME UPDATE
        print("cpycwpclt.time_update:\n")

        pycwpclt.time_update(code_names[i_rank], 1.5)

        comm.Barrier()

        # USER TGT PTS
        coord = np.array([6, 7, 8, 9, 10, 11], dtype=np.double)
        print("cpl.user_tgt_pts_set:\n")

        cpl.user_tgt_pts_set(0,
                             2,
                             coord,
                             None)

        print("cpl.mesh_interf_del:\n")

        cpl.mesh_interf_del()

    else:
        # STD MESH
        print("cpl.mesh_interf_block_add:\n")

        block_id = cpl.mesh_interf_block_add(pycwpclt.BLOCK_FACE_TRIA3)
        print("cpl.mesh_interf_block_std_set:\n")

        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)
        cpl.mesh_interf_block_std_set(0,
                                      block_id,
                                      2,
                                      connec,
                                      None)

        print("cpl.mesh_interf_block_std_get:\n")

        out = cpl.mesh_interf_block_std_get(0,
                                            block_id)
        print("  - n_elts : {param}\n".format(param=out["n_elts"]))
        print("  - connec {param}\n".format(param=out["connec"]))
        print("  - global_num : {param}\n".format(param=out["global_num"]))


        print("cpl.mesh_interf_from_faceedge_set:\n")

        face_edge_idx = np.array([0, 3, 6], dtype=np.int32)
        face_edge = np.array([1, 2, 3, 5, 4, 3], dtype=np.int32)
        edge_vtx_idx = np.array([0, 2, 4, 6, 8, 10], dtype=np.int32)
        edge_vtx = np.array([3, 1, 1, 2, 2, 3, 4, 2, 3, 4], dtype=np.int32)
        cpl.mesh_interf_from_faceedge_set(0,
                                          2,
                                          face_edge_idx,
                                          face_edge,
                                          5,
                                          edge_vtx_idx,
                                          edge_vtx,
                                          None)

        print("cpl.mesh_interf_del:\n")

        cpl.mesh_interf_del()

    # Volumic Cpl
    print("pycwpclt.Coupling:\n")

    cpl2 = pycwpclt.Coupling(code_names[i_rank],
                         "test_vol",
                         code_names[(i_rank+1)%2],
                         pycwpclt.INTERFACE_VOLUME,
                         pycwpclt.COMM_PAR_WITH_PART,
                         pycwpclt.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                         1,
                         pycwpclt.DYNAMIC_MESH_STATIC,
                         pycwpclt.TIME_EXCH_USER_CONTROLLED)

    block_id = cpl2.mesh_interf_block_add(pycwpclt.BLOCK_CELL_POLY)

    connec_faces_idx = np.array([0, 3, 6, 9, 12], dtype=np.int32)
    connec_faces = np.array([1, 2, 3, 1, 2, 4, 2, 3, 4, 1, 3, 4], dtype=np.int32)
    connec_cells_idx = np.array([0, 4], dtype=np.int32)
    connec_cells = np.array([1, 2, 3, 4], dtype=np.int32)

    print("cpl2.mesh_interf_c_poly_block_set:\n")

    cpl2.mesh_interf_c_poly_block_set(0,
                                     block_id,
                                     1,
                                     4,
                                     connec_faces_idx,
                                     connec_faces,
                                     connec_cells_idx,
                                     connec_cells,
                                     None)

    print("cpl2.mesh_interf_c_poly_block_get:\n")

    out = cpl2.mesh_interf_c_poly_block_get(0,
                                           block_id)
    print("  - n_elts : {param}\n".format(param=out["n_elts"]))
    print("  - n_faces : {param}\n".format(param=out["n_faces"]))
    print("  - connec_faces_idx {param}\n".format(param=out["connec_faces_idx"]))
    print("  - connec_faces {param}\n".format(param=out["connec_faces"]))
    print("  - connec_cells_idx {param}\n".format(param=out["connec_cells_idx"]))
    print("  - connec_cells {param}\n".format(param=out["connec_cells"]))

    face_vtx_idx = np.array([0, 3, 6, 9, 12], dtype=np.int32)
    face_vtx = np.array([1, 2, 3, 1, 2, 4, 2, 3, 4, 1, 3, 4], dtype=np.int32)
    cell_face_idx = np.array([0, 4], dtype=np.int32)
    cell_face = np.array([1, 2, 3, 4], dtype=np.int32)

    print("cpl2.mesh_interf_from_cellface_set:\n")

    cpl2.mesh_interf_from_cellface_set(0,
                                      1,
                                      cell_face_idx,
                                      cell_face,
                                      4,
                                      face_vtx_idx,
                                      face_vtx,
                                      None)

    print("cpl2.mesh_interf_del:\n")

    cpl2.mesh_interf_del()

    comm.Barrier()

    # FINALIZE
    pycwpclt.finalize()

    # END
    print("\nEnd.\n")
    comm.Barrier()
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
