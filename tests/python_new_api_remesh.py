#!/usr/bin/env python
#-----------------------------------------------------------------------------
# This file is part of the CWIPI library.
#
# Copyright (C) 2024  ONERA
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


def eval_field1(coord, time):
  return np.cos( np.sqrt(coord[0::3]**2 + coord[1::3]**2) - 0.1*time )

def eval_field2(elt_vtx_idx, elt_vtx, coord, time):
  vtx_field_value = np.cos(coord[0::3] + 0.1*time)

  n_elt = len(elt_vtx_idx) - 1
  elt_field_value = np.empty(n_elt, dtype=np.double)
  for i_elt in range(n_elt):
    l_vtx = elt_vtx[elt_vtx_idx[i_elt]:elt_vtx_idx[i_elt+1]] - 1
    elt_field_value[i_elt] = np.sum(vtx_field_value[l_vtx])

  elt_field_value /= np.diff(elt_vtx_idx)

  return elt_field_value


def run_test():
  # Initialize MPI
  comm = MPI.COMM_WORLD

  # Load Python CWIPI module
  try:
    from pycwp import pycwp
  except:
    if comm.rank == 0:
      print("      Error : CWIPI module not found (update PYTHONPATH variable)")
      print(f"cwp : {pycwp.__file__}")
      sys.exit(1)

  # Load Python CWIPI Test module
  try:
    from pycwpt import pycwpt
  except:
    if comm.rank == 0:
      print("      Error : CWIPI Test module not found (update PYTHONPATH variable)")
      print(f"cwpt : {pycwpt.__file__}")
      sys.exit(1)


  # Even ranks run code0, odd ranks run code1
  i_code = comm.rank % 2
  code_name = f"code{i_code}"

  # Coupled with the other code
  coupled_code_name = f"code{(i_code+1)%2}"

  # Initialize CWIPI
  is_active_rank = True

  intra_comm = pycwp.init(comm,
                          [code_name],
                          is_active_rank)

  # Create coupling environment
  cpl = pycwp.Coupling(code_name,
                       "python_new_api_remesh",
                       coupled_code_name,
                       pycwp.INTERFACE_SURFACE,
                       pycwp.COMM_PAR_WITH_PART,
                       pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                       1,
                       pycwp.DYNAMIC_MESH_VARIABLE,
                       pycwp.TIME_EXCH_USER_CONTROLLED)

  # Set visu status
  cpl.visu_set(1,
               pycwp.VISU_FORMAT_ENSIGHT,
               "text")

  # Initial interface mesh
  n_vtx_seg = 10 - 3*i_code

  mesh = pycwpt.generate_mesh_rectangle_simplified(intra_comm[0], n_vtx_seg)

  cpl.mesh_interf_vtx_set(0,
                          mesh["coords"],
                          None)

  block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

  cpl.mesh_interf_f_poly_block_set(0,
                                   block_id,
                                   mesh["elt_vtx_idx"],
                                   mesh["elt_vtx"],
                                   None)

  cpl.mesh_interf_finalize()

  # Create fields
  field_name1 = "field1"
  field_name2 = "field2"

  if i_code == 0:
    # Code0 is node-centered
    field1 = cpl.field_create(field_name1,
                              pycwp.DOUBLE,
                              pycwp.FIELD_STORAGE_INTERLACED,
                              1,
                              pycwp.DOF_LOCATION_NODE,
                              pycwp.FIELD_EXCH_SEND,
                              pycwp.STATUS_ON)

    field2 = cpl.field_create(field_name2,
                              pycwp.DOUBLE,
                              pycwp.FIELD_STORAGE_INTERLACED,
                              1,
                              pycwp.DOF_LOCATION_NODE,
                              pycwp.FIELD_EXCH_RECV,
                              pycwp.STATUS_ON)

    n_vtx = len(mesh["coords"]) // 3
    send_data1 = np.empty(n_vtx, dtype=np.double)
    recv_data2 = np.empty(n_vtx, dtype=np.double)

  else:
    # Code1 is cell-centered
    field1 = cpl.field_create(field_name1,
                              pycwp.DOUBLE,
                              pycwp.FIELD_STORAGE_INTERLACED,
                              1,
                              pycwp.DOF_LOCATION_CELL_CENTER,
                              pycwp.FIELD_EXCH_RECV,
                              pycwp.STATUS_ON)

    field2 = cpl.field_create(field_name2,
                              pycwp.DOUBLE,
                              pycwp.FIELD_STORAGE_INTERLACED,
                              1,
                              pycwp.DOF_LOCATION_CELL_CENTER,
                              pycwp.FIELD_EXCH_SEND,
                              pycwp.STATUS_ON)

    n_elt = len(mesh["elt_vtx_idx"]) - 1
    recv_data1 = np.empty(n_elt, dtype=np.double)
    send_data2 = np.empty(n_elt, dtype=np.double)


  # Time loop
  n_time_steps = 30
  time = 0.0
  freq_remesh = 5 + 3*i_code


  for i_step in range(n_time_steps):

    if comm.rank == 0:
      print(f"Step {i_step}", flush=True)

    # Begin time step
    pycwp.time_step_beg(code_name,
                        time)

    # Remesh if required
    if (i_step+1)%freq_remesh == 0:
      if intra_comm[0].rank == 0:
        print(f"  Remesh {code_name}", flush=True)


      #  Delete current interface mesh
      cpl.mesh_interf_del()


      #  Generate new interface mesh
      n_vtx_seg += 2 + 2*i_code

      mesh = pycwpt.generate_mesh_rectangle_simplified(intra_comm[0], n_vtx_seg)


      #  Set new interface mesh
      cpl.mesh_interf_vtx_set(0,
                              mesh["coords"],
                              None)

      block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

      cpl.mesh_interf_f_poly_block_set(0,
                                       block_id,
                                       mesh["elt_vtx_idx"],
                                       mesh["elt_vtx"],
                                       None)

      cpl.mesh_interf_finalize()


      #  Reallocate array of field values
      if i_code == 0:
        n_vtx = len(mesh["coords"]) // 3
        send_data1 = np.empty(n_vtx, dtype=np.double)
        recv_data2 = np.empty(n_vtx, dtype=np.double)
      else:
        n_elt = len(mesh["elt_vtx_idx"]) - 1
        recv_data1 = np.empty(n_elt, dtype=np.double)
        send_data2 = np.empty(n_elt, dtype=np.double)


    # Set pointers to array of field values
    if i_code == 0:
      send_data1[:] = eval_field1(mesh["coords"], time)
      field1.data_set(0,
                      pycwp.FIELD_MAP_SOURCE,
                      send_data1)
      field2.data_set(0,
                      pycwp.FIELD_MAP_TARGET,
                      recv_data2)
    else:
      send_data2[:] = eval_field2(mesh["elt_vtx_idx"], mesh["elt_vtx"], mesh["coords"], time)
      field1.data_set(0,
                      pycwp.FIELD_MAP_TARGET,
                      recv_data1)
      field2.data_set(0,
                      pycwp.FIELD_MAP_SOURCE,
                      send_data2)


    # Compute spatial interpolation weights
    cpl.spatial_interp_weights_compute()


    # Exchange fields
    if i_code == 0:
      field1.issend()
      field2.irecv()
    else:
      field1.irecv()
      field2.issend()

    # Overlap exchanges with computations if possible...

    # Finalize field exchange
    if i_code == 0:
      field1.wait_issend()
      field2.wait_irecv()
    else:
      field1.wait_irecv()
      field2.wait_issend()


    # End time step
    pycwp.time_step_end(code_name)

    time += 1.0

  # Finalize CWIPI
  pycwp.finalize()


if __name__ == '__main__':
  run_test()

