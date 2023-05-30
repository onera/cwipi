#!/usr/bin/env python
#-----------------------------------------------------------------------------
# This file is part of the CWIPI library.
#
# Copyright (C) 2023  ONERA
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

def first_interpolation(local_code_name,
                        cpl_id,
                        field_id,
                        i_part,
                        spatial_interp_algorithm,
                        storage,
                        buffer_in,
                        buffer_out):

  try:
    from pycwp import pycwp
  except:
    if i_rank == 0:
      print("      Error : CWIPI module not found (update PYTHONPATH variable)")
      print(f"cwp : {pycwp.__file__}")
      sys.exit(1)

  n_comp = pycwp.interp_field_n_components_get(local_code_name,
                                               cpl_id,
                                               field_id)

  if spatial_interp_algorithm == pycwp.SPATIAL_INTERP_FROM_CLOSEST_SOURCES_LEAST_SQUARES:
    tgt_data = pycwp.interp_tgt_data_get(local_code_name,
                                         cpl_id,
                                         field_id,
                                         i_part)
    n_tgt          = tgt_data["n_elt_tgt"]
    ref_tgt        = tgt_data["referenced_tgt"]
    tgt_to_src_idx = tgt_data["tgt_come_from_src_idx"]

    distance2 = pycwp.interp_closest_points_distances_get(local_code_name,
                                                          cpl_id,
                                                          field_id,
                                                          i_part)

    for i, jtgt in enumerate(ref_tgt):
      itgt = jtgt-1
      buffer_out[n_comp*itgt:n_comp*(itgt+1)] = 0
      sum_w = 0
      for isrc in range(tgt_to_src_idx[i], tgt_to_src_idx[i+1]):
        w = 1./max(1.e-24, distance2[isrc])
        sum_w += w
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] + \
        w*buffer_in[n_comp*isrc:n_comp*(isrc+1)]

      buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
      buffer_out[n_comp*itgt:n_comp*(itgt+1)] / sum_w

  elif spatial_interp_algorithm == pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE:
    src_data = pycwp.interp_src_data_get(local_code_name,
                                         cpl_id,
                                         field_id,
                                         i_part)
    n_src = src_data["n_elt_src"]
    src_to_tgt_idx = src_data["src_to_tgt_idx"]

    weight = pycwp.interp_location_weights_get(local_code_name,
                                               cpl_id,
                                               field_id,
                                               i_part)

    cell_data = pycwp.interp_location_internal_cell_vtx_get(local_code_name,
                                                            cpl_id,
                                                            field_id,
                                                            i_part)
    cell_vtx_idx = cell_data["cell_vtx_idx"]
    cell_vtx     = cell_data["cell_vtx"]

    for isrc in range(n_src):
      for itgt in range(src_to_tgt_idx[isrc], src_to_tgt_idx[isrc+1]):
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] = 0
        for i in range(cell_vtx_idx[isrc], cell_vtx_idx[isrc+1]):
          ivtx = cell_vtx[i] - 1
          w = weight[ivtx]
          buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
          buffer_out[n_comp*itgt:n_comp*(itgt+1)] + \
          w*buffer_in[n_comp*ivtx:n_comp*(ivtx+1)]

  else:
    print(f"      Error : wrong spatial_interp_algorithm ({spatial_interp_algorithm})")
    sys.exit(1)


def second_interpolation(local_code_name,
                         cpl_id,
                         field_id,
                         i_part,
                         spatial_interp_algorithm,
                         storage,
                         buffer_in,
                         buffer_out):

  try:
    from pycwp import pycwp
  except:
    if i_rank == 0:
      print("      Error : CWIPI module not found (update PYTHONPATH variable)")
      print(f"cwp : {pycwp.__file__}")
      sys.exit(1)

  n_comp = pycwp.interp_field_n_components_get(local_code_name,
                                               cpl_id,
                                               field_id)

  if spatial_interp_algorithm == pycwp.SPATIAL_INTERP_FROM_CLOSEST_SOURCES_LEAST_SQUARES:
    tgt_data = pycwp.interp_tgt_data_get(local_code_name,
                                         cpl_id,
                                         field_id,
                                         i_part)
    n_tgt          = tgt_data["n_elt_tgt"]
    ref_tgt        = tgt_data["referenced_tgt"]
    tgt_to_src_idx = tgt_data["tgt_come_from_src_idx"]

    distance2 = pycwp.interp_closest_points_distances_get(local_code_name,
                                                          cpl_id,
                                                          field_id,
                                                          i_part)

    for i, jtgt in enumerate(ref_tgt):
      itgt = jtgt-1
      buffer_out[n_comp*itgt:n_comp*(itgt+1)] = 0
      sum_w = 0
      for isrc in range(tgt_to_src_idx[i], tgt_to_src_idx[i+1]):
        w = 1./max(1.e-24, distance2[isrc])
        sum_w += w
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] + \
        w*buffer_in[n_comp*isrc:n_comp*(isrc+1)]

      buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
      buffer_out[n_comp*itgt:n_comp*(itgt+1)] / sum_w

  elif spatial_interp_algorithm == pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE:
    src_data = pycwp.interp_src_data_get(local_code_name,
                                         cpl_id,
                                         field_id,
                                         i_part)
    n_src = src_data["n_elt_src"]
    src_to_tgt_idx = src_data["src_to_tgt_idx"]

    weight = pycwp.interp_location_weights_get(local_code_name,
                                               cpl_id,
                                               field_id,
                                               i_part)

    cell_data = pycwp.interp_location_internal_cell_vtx_get(local_code_name,
                                                            cpl_id,
                                                            field_id,
                                                            i_part)
    cell_vtx_idx = cell_data["cell_vtx_idx"]
    cell_vtx     = cell_data["cell_vtx"]

    for isrc in range(n_src):
      for itgt in range(src_to_tgt_idx[isrc], src_to_tgt_idx[isrc+1]):
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] = 0
        for i in range(cell_vtx_idx[isrc], cell_vtx_idx[isrc+1]):
          ivtx = cell_vtx[i] - 1
          w = weight[ivtx]
          buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
          buffer_out[n_comp*itgt:n_comp*(itgt+1)] + \
          w*buffer_in[n_comp*ivtx:n_comp*(ivtx+1)]

  else:
    print(f"      Error : wrong spatial_interp_algorithm ({spatial_interp_algorithm})")
    sys.exit(1)


def runTest():
    """
    Run tests on Python interface of new API
    """
    global f
    comm = MPI.COMM_WORLD

    i_rank = comm.rank
    n_rank = comm.size

    if (i_rank == 0):
        print("\nSTART: python_new_api_multiple_field_with_callback.py")

    # Load Python PDM module
    try:
      from Pypdm import Pypdm
    except:
      if i_rank == 0:
        print("      Error : PDM module not found (update PYTHONPATH variable)")
        print(f"pdm : {Pypdm.__file__}")
        sys.exit(1)

    # Load Python CWIPI module
    try:
      from pycwp import pycwp
    except:
      if i_rank == 0:
        print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        print(f"cwp : {pycwp.__file__}")
        sys.exit(1)

    # CODE
    proc0 = 0
    proc1 = 0
    if (i_rank < n_rank // 2):
        proc0 = 1
    else:
        proc1 = 1

    # INIT
    if (proc0):
        code_name         = "proc0"
        coupled_code_name = "proc1"

    if (proc1):
        code_name         = "proc1"
        coupled_code_name = "proc0"

    n_code = 1
    is_active_rank = [True for i in range(n_code)]
    intra_comm = pycwp.init(comm,
                            [code_name],
                            is_active_rank)

    # COUPLING
    cpl = pycwp.Coupling(code_name,
                         "coupling_multiple_field_with_callback",
                         coupled_code_name,
                         pycwp.INTERFACE_SURFACE,
                         pycwp.COMM_PAR_WITH_PART,
                         pycwp.SPATIAL_INTERP_FROM_CLOSEST_SOURCES_LEAST_SQUARES,
                         1,
                         pycwp.DYNAMIC_MESH_STATIC,
                         pycwp.TIME_EXCH_USER_CONTROLLED)

    # MESH
    mesh = Pypdm.generate_mesh_rectangle_simplified(intra_comm[0],
                                                    5)

    cpl.mesh_interf_vtx_set(0,
                            mesh["n_vtx"],
                            mesh["coords"],
                            None)

    block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

    cpl.mesh_interf_f_poly_block_set(0,
                                     block_id,
                                     mesh["n_elt"],
                                     mesh["elt_vtx_idx"],
                                     mesh["elt_vtx"],
                                     None)

    cpl.mesh_interf_finalize()

    # FIELD 1 - x
    field1_name = "Field 1"

    send_field1_data = np.zeros(mesh["n_vtx"], dtype=np.double)
    for i in range(mesh["n_vtx"]):
        send_field1_data[i] = mesh["coords"][3*i]
    recv_field1_data = np.zeros(mesh["n_vtx"], dtype=np.double)

    if (proc0) :
        exchange_type = pycwp.FIELD_EXCH_SEND
    else :
        exchange_type = pycwp.FIELD_EXCH_RECV

    cpl.field_create(field1_name,
                     pycwp.DOUBLE,
                     pycwp.FIELD_STORAGE_INTERLACED,
                     1,
                     pycwp.DOF_LOCATION_NODE,
                     exchange_type,
                     pycwp.STATUS_OFF)

    if (proc0) :
        cpl.field_set(field1_name,
                      0,
                      pycwp.FIELD_MAP_SOURCE,
                      send_field1_data)
    else :
        cpl.field_set(field1_name,
                      0,
                      pycwp.FIELD_MAP_TARGET,
                      recv_field1_data)

        # USER FUNCTION
        cpl.field_interp_function_set(field1_name,
                                      first_interpolation)

    # FIELD 2 - y
    field2_name = "Field 2"

    send_field2_data = np.zeros(mesh["n_vtx"], dtype=np.double)
    for i in range(mesh["n_vtx"]):
        send_field2_data[i] = mesh["coords"][3*i+1]
    recv_field2_data = np.zeros(mesh["n_vtx"], dtype=np.double)

    if (proc0) :
        exchange_type = pycwp.FIELD_EXCH_SEND
    else :
        exchange_type = pycwp.FIELD_EXCH_RECV

    cpl.field_create(field2_name,
                     pycwp.DOUBLE,
                     pycwp.FIELD_STORAGE_INTERLACED,
                     1,
                     pycwp.DOF_LOCATION_NODE,
                     exchange_type,
                     pycwp.STATUS_OFF)

    if (proc0) :
        cpl.field_set(field2_name,
                      0,
                      pycwp.FIELD_MAP_SOURCE,
                      send_field2_data)
    else :
        cpl.field_set(field2_name,
                      0,
                      pycwp.FIELD_MAP_TARGET,
                      recv_field2_data)

        # USER FUNCTION
        cpl.field_interp_function_set(field2_name,
                                      second_interpolation)

    # INTERPOLATION
    cpl.spatial_interp_property_set("n_closest_pts",
                                    "int",
                                    "1")

    cpl.spatial_interp_weights_compute()

    # EXCHANGE
    if (proc0) :
        cpl.field_issend(field1_name)
        cpl.field_issend(field2_name)
    else :
        cpl.field_irecv (field1_name)
        cpl.field_irecv (field2_name)

    if (proc0) :
        cpl.field_wait_issend(field1_name)
        cpl.field_wait_issend(field2_name)
    else :
        cpl.field_wait_irecv (field1_name)
        cpl.field_wait_irecv (field2_name)

    # CHECK
    for i in range(mesh["n_vtx"]):
        if (proc1) :
            egal1 = (recv_field1_data[i] == mesh["coords"][3*i])
            if (egal1 == 0) :
                sys.exit(1)

            egal2 = (recv_field2_data[i] == mesh["coords"][3*i+1])
            if (egal2 == 0) :
                sys.exit(1)

    # FINALIZE
    cpl.mesh_interf_del()
    pycwp.finalize()

    # END
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
