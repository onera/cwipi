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

def userInterp(interface_type,
               code_name,
               src_n_block,
               src_blocks_type,
               src_i_part,
               src_n_vtx,
               src_vtx_coords,
               src_vtx_global_num,
               src_n_elts,
               src_id_block,
               src_elt_in_block,
               src_elt_vtx_idx,
               src_elt_vtx,
               src_elts_global_num,
               tgt_n_pts,
               tgt_pts_elt_idx,
               tgt_pts_coords,
               tgt_pts_dist,
               tgt_pts_uvw,
               tgt_pts_weights_idx,
               tgt_pts_weights,
               stride,
               src_field_dof_location,
               src_field,
               tgt_field):
    """
    User interpolation method
    """
    global f

    f.write("in a user interpolation function !")

    if (distant_field is not None):
        distant_field = 1.234

def runTest():
    """
    Run tests on Python interface of new API
    """
    global f
    comm = MPI.COMM_WORLD

    i_rank = comm.rank
    n_rank = comm.size

    if (i_rank == 0):
        print("\nSTART: python_new_api.py")

    if (n_rank != 2):
        if i_rank == 0:
            print("      Not executed : only available for 2 processes")
        return

    code_names = ["proc0","proc1"]

    if (i_rank == 0):
        code_name = ["proc0"]

    if (i_rank == 1):
        code_name = ["proc1"]

    try:
        from pycwp import pycwp
    except:
        if i_rank == 0:
            print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        sys.exit(1)

    # OUTPUT
    srank = '{0}'.format(i_rank)
    f=open("python_new_api_"+srank.zfill(4)+".txt",'w')
    pycwp.output_file_set(f)

    # INIT
    f.write("pycwp.init:\n")
    n_code = 1
    is_active_rank = np.array([1], dtype=np.int32)
    time_init = np.array([0.], dtype=np.double)
    out = pycwp.init(comm,
                     n_code,
                     code_name,
                     is_active_rank,
                     time_init)
    f.write("  - intra_comms : {param}\n".format(param=out[0]))

    # STATE UPDATE
    pycwp.state_update(code_names[i_rank], pycwp.STATE_IN_PROGRESS)
    f.write("pycwp.state_get:\n")
    state = pycwp.state_get(code_names[i_rank])
    f.write("  - state : {param}\n".format(param=state))
    pycwp.state_update(code_names[i_rank], pycwp.STATE_END)
    f.flush()

    # PROPERTIES DUMP
    f.write("pycwp.properties_dump:\n")
    pycwp.properties_dump()

    # CODES
    f.write("pycwp.code:\n")
    n_code = pycwp.codes_nb_get()
    code   = pycwp.codes_list_get()
    n_loc_code = pycwp.loc_codes_nb_get()
    loc_code   = pycwp.loc_codes_list_get()
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

    comm.Barrier()

    f.write("cwp.param_get ({param}):\n".format(param=i_rank))
    value = pycwp.param_get(code_names[i_rank], "double", pycwp.DOUBLE)
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

    comm.Barrier()

    f.write("pycwp.param_get ({param}):\n".format(param=i_rank))
    value = pycwp.param_get(code_names[i_rank], "double", pycwp.DOUBLE)
    f.write("  - value (1): {param}\n".format(param=value))

    pycwp.param_lock(code_names[i_rank])
    pycwp.param_del(code_names[i_rank], "str", pycwp.CHAR)
    pycwp.param_unlock(code_names[i_rank])

    comm.Barrier()

    f.write("pycwp.param_n_get:\n")
    n_param_str = pycwp.param_n_get(code_names[i_rank], pycwp.CHAR)
    n_param_int = pycwp.param_n_get(code_names[i_rank], pycwp.INT)
    f.write("  - n_param_str: {param}\n".format(param=n_param_str))
    f.write("  - n_param_int: {param}\n".format(param=n_param_int))

    f.write("pycwp.param_list_get:\n")
    str_param = pycwp.param_list_get(code_names[i_rank], pycwp.CHAR)
    for i in range(str_param['n_param']):
        f.write("    --> str_param: {param}\n".format(param=str_param['param_names'][i]))

    f.write("pycwp.param_is:\n")
    bool_int = pycwp.param_is(code_names[i_rank], "entier", pycwp.INT)
    f.write("  - bool_int 'entier': {param}\n".format(param=bool_int))
    bool_int = pycwp.param_is(code_names[i_rank], "chapeau", pycwp.INT)
    f.write("  - bool_int 'chapeau': {param}\n".format(param=bool_int))

    comm.Barrier()

    f.write("pycwp.param_list_get:\n")
    int_param = pycwp.param_list_get(code_names[i_rank], pycwp.INT)
    for i in range(int_param['n_param']):
        f.write("    --> int_param: {param}\n".format(param=int_param['param_names'][i]))

    f.write("pycwp.param_get ({param}):\n".format(param=i_rank))
    value = pycwp.param_get(code_names[i_rank], "entier", pycwp.INT)
    f.write("  - value int: {param}\n".format(param=value))

    comm.Barrier()

    f.write("pycwp.param_reduce:\n")
    result = pycwp.param_reduce(pycwp.OP_MIN, "entier",  pycwp.INT, 2, code_names)
    f.write("  - result: {param}\n".format(param=result))
    f.flush()

    # Cpl
    f.write("pycwp.Coupling:\n")
    f.flush()
    cpl = pycwp.Coupling(code_names[i_rank],
                         "test",
                         code_names[(i_rank+1)%2],
                         pycwp.INTERFACE_SURFACE, #INTERFACE_LINEAR,
                         pycwp.COMM_PAR_WITH_PART,
                         pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                         1,
                         pycwp.DYNAMIC_MESH_VARIABLE,
                         pycwp.TIME_EXCH_USER_CONTROLLED)

    #if SPATIAL_INTERP_FROM_CLOSEST_POINT_LEAST_SQUARES, get munmap_chunk(): invalid pointer
    #if SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, blocked in weigths_compute

    # VISU
    cpl.visu_set(1,
                 pycwp.VISU_FORMAT_ENSIGHT,
                 "text")

    # std or polygon
    polygon = 1

    # MESH
    if (polygon):
        if (i_rank == 0):
            coord = np.array([0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0], dtype=np.double)
            connec_idx = np.array([0, 3, 6], dtype=np.int32)
            connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

        if (i_rank == 1):
            coord = np.array([0, 1, 0, 0, 2, 0, 1, 1, 0, 1, 2, 0], dtype=np.double)
            connec_idx = np.array([0, 3, 6], dtype=np.int32)
            connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

        f.write("cpl.mesh_interf_vtx_set:\n")
        f.flush()
        cpl.mesh_interf_vtx_set(0,
                                4,
                                coord,
                                None)

        f.write("cpl.mesh_interf_block_add:\n")
        f.flush()
        block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)
        f.write("cpl.mesh_interf_f_poly_block_set ({param}):\n".format(param=i_rank))
        f.flush()

        cpl.mesh_interf_f_poly_block_set(0,
                                         block_id,
                                         2,
                                         connec_idx,
                                         connec,
                                         None)

        f.write("cpl.mesh_interf_finalize:\n")
        f.flush()
        cpl.mesh_interf_finalize()

        f.write("cpl.mesh_interf_f_poly_block_get:\n")
        f.flush()
        out = cpl.mesh_interf_f_poly_block_get(0, block_id)
        f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
        f.write("  - connec_idx {param}\n".format(param=out["connec_idx"]))
        f.write("  - connec {param}\n".format(param=out["connec"]))
        f.write("  - global_num : {param}\n".format(param=out["global_num"]))
        f.flush()

        # FIELD

        # USER INTERPOLATION to do has to be done
        # f.write("cpl.interp_from_location_set:\n")
        # f.flush()
        # cpl.interp_from_location_set("champs", userInterp)

        sendField = np.array([0.0, 0.1, 0.2, 0.3], dtype=np.double)
        recvField = np.arange(4, dtype=np.double)

        if (i_rank == 0):
            cpl.field_create("champs",
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_NODE,
                             pycwp.FIELD_EXCH_SEND,
                             pycwp.STATUS_OFF)
            cpl.field_set("champs",
                          0,
                          pycwp.FIELD_MAP_SOURCE,
                          sendField)

        if (i_rank == 1):
            cpl.field_create("champs",
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_NODE,
                             pycwp.FIELD_EXCH_RECV,
                             pycwp.STATUS_OFF)
            cpl.field_set("champs",
                          0,
                          pycwp.FIELD_MAP_TARGET,
                          recvField)

        comm.Barrier()

        f.write("cpl.spatial_interp_property_set:\n")
        f.flush()
        cpl.spatial_interp_property_set("tolerance", "double", "1e-2")

        comm.Barrier()

        f.write("cpl.spatial_interp_weights_compute:\n")
        f.flush()
        cpl.spatial_interp_weights_compute()

        if (i_rank == 0):
            f.write("pycwp.field_issend (0):\n")
            f.flush()
            pycwp.field_issend(code_names[i_rank], "test", "champs")

        if (i_rank == 1):
            f.write("pycwp.field_irecv (1):\n")
            f.flush()
            pycwp.field_irecv(code_names[i_rank], "test", "champs")

        if (i_rank == 0):
            f.write("pycwp.field_wait_issend (0):\n")
            f.flush()
            pycwp.field_wait_issend(code_names[i_rank], "test", "champs")

        if (i_rank == 1):
            f.write("pycwp.field_wait_irecv (1):\n")
            f.flush()
            pycwp.field_wait_irecv(code_names[i_rank], "test", "champs")

        comm.Barrier()

        f.write("cpl.field_del:\n")
        f.flush()
        cpl.field_del("champs")

        # USER INTERPOLATION to do has to be done
        # f.write("cpl.interp_from_location_unset:\n")
        # f.flush()
        # cpl.interp_from_location_unset("champs")

        # TIME UPDATE
        f.write("cpycwp.time_update:\n")
        f.flush()
        pycwp.time_update(code_names[i_rank], 1.5)

        comm.Barrier()

        # USER TGT PTS
        coord = np.array([6, 7, 8, 9, 10, 11], dtype=np.double)
        f.write("cpl.user_tgt_pts_set:\n")
        f.flush()
        cpl.user_tgt_pts_set(0,
                             2,
                             coord,
                             None)

        # USER STRUCTURE
        class userClass:
            animal = "chat"
            aliment = "aligot"
            ville = "Toulouse"

        userObj = userClass()

        f.write("pycwp.user_structure_set:\n")
        f.flush()
        pycwp.user_structure_set(code_names[i_rank],
                                 userObj)
        f.write("pycwp.user_structure_get:\n")
        f.flush()
        user_structure = pycwp.user_structure_get(code_names[i_rank])

        print(user_structure.animal)
        print(user_structure.aliment)
        print(user_structure.ville)

        f.write("cpl.mesh_interf_del:\n")
        f.flush()
        cpl.mesh_interf_del()

    else:
        # STD MESH
        f.write("cpl.mesh_interf_block_add:\n")
        f.flush()
        block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_TRIA3)
        f.write("cpl.mesh_interf_block_std_set:\n")
        f.flush()
        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)
        cpl.mesh_interf_block_std_set(0,
                                      block_id,
                                      2,
                                      connec,
                                      None)

        f.write("cpl.mesh_interf_block_std_get:\n")
        f.flush()
        out = cpl.mesh_interf_block_std_get(0,
                                            block_id)
        f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
        f.write("  - connec {param}\n".format(param=out["connec"]))
        f.write("  - global_num : {param}\n".format(param=out["global_num"]))
        f.flush()

        f.write("cpl.mesh_interf_from_faceedge_set:\n")
        f.flush()
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

        f.write("cpl.mesh_interf_del:\n")
        f.flush()
        cpl.mesh_interf_del()

    # Volumic Cpl
    f.write("pycwp.Coupling:\n")
    f.flush()
    cpl2 = pycwp.Coupling(code_names[i_rank],
                         "test_vol",
                         code_names[(i_rank+1)%2],
                         pycwp.INTERFACE_VOLUME,
                         pycwp.COMM_PAR_WITH_PART,
                         pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                         1,
                         pycwp.DYNAMIC_MESH_STATIC,
                         pycwp.TIME_EXCH_USER_CONTROLLED)

    block_id = cpl2.mesh_interf_block_add(pycwp.BLOCK_CELL_POLY)

    connec_faces_idx = np.array([0, 3, 6, 9, 12], dtype=np.int32)
    connec_faces = np.array([1, 2, 3, 1, 2, 4, 2, 3, 4, 1, 3, 4], dtype=np.int32)
    connec_cells_idx = np.array([0, 4], dtype=np.int32)
    connec_cells = np.array([1, 2, 3, 4], dtype=np.int32)

    f.write("cpl2.mesh_interf_c_poly_block_set:\n")
    f.flush()
    cpl2.mesh_interf_c_poly_block_set(0,
                                     block_id,
                                     1,
                                     4,
                                     connec_faces_idx,
                                     connec_faces,
                                     connec_cells_idx,
                                     connec_cells,
                                     None)

    f.write("cpl2.mesh_interf_c_poly_block_get:\n")
    f.flush()
    out = cpl2.mesh_interf_c_poly_block_get(0,
                                           block_id)
    f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
    f.write("  - n_faces : {param}\n".format(param=out["n_faces"]))
    f.write("  - connec_faces_idx {param}\n".format(param=out["connec_faces_idx"]))
    f.write("  - connec_faces {param}\n".format(param=out["connec_faces"]))
    f.write("  - connec_cells_idx {param}\n".format(param=out["connec_cells_idx"]))
    f.write("  - connec_cells {param}\n".format(param=out["connec_cells"]))

    face_vtx_idx = np.array([0, 3, 6, 9, 12], dtype=np.int32)
    face_vtx = np.array([1, 2, 3, 1, 2, 4, 2, 3, 4, 1, 3, 4], dtype=np.int32)
    cell_face_idx = np.array([0, 4], dtype=np.int32)
    cell_face = np.array([1, 2, 3, 4], dtype=np.int32)

    f.write("cpl2.mesh_interf_from_cellface_set:\n")
    f.flush()
    cpl2.mesh_interf_from_cellface_set(0,
                                      1,
                                      cell_face_idx,
                                      cell_face,
                                      4,
                                      face_vtx_idx,
                                      face_vtx,
                                      None)

    f.write("cpl2.mesh_interf_del:\n")
    f.flush()
    cpl2.mesh_interf_del()

    # FINALIZE
    pycwp.finalize()

    # END
    f.write("\nEnd.\n")
    f.close()
    comm.Barrier()
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
