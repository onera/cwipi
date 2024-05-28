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

def generate_mesh(intra_comm, n_vtx_seg):
  if intra_comm.rank != 0:
    return {
    "n_vtx"        : 0,
    "n_elt"        : 0,
    "vtx_coord"    : np.empty(0, dtype=np.double),
    "vtx_ln_to_gn" : None,
    "connec"       : np.empty(0, dtype=np.int32),
    "elt_ln_to_gn" : None
    }

  n_vtx = n_vtx_seg**2
  n_elt = (n_vtx_seg-1)**2

  # Vertices
  vtx_coord = np.empty(3*n_vtx, dtype=np.double)
  vtx_ln_to_gn = 1 + np.arange(n_vtx, dtype=np.int64)

  step = 1./(n_vtx_seg - 1)

  k = 0
  for j in range(n_vtx_seg):
    for i in range(n_vtx_seg):
      vtx_coord[3*k  ] = i*step
      vtx_coord[3*k+1] = j*step
      vtx_coord[3*k+2] = 0
      k += 1

  # Elements
  connec = np.empty(4*n_elt, dtype=np.int32)
  elt_ln_to_gn = 1 + np.arange(n_elt, dtype=np.int64)

  k = 0
  for j in range(n_vtx_seg-1):
    for i in range(n_vtx_seg-1):
      connec[4*k  ] = n_vtx_seg*j     + i     + 1
      connec[4*k+1] = n_vtx_seg*j     + (i+1) + 1
      connec[4*k+2] = n_vtx_seg*(j+1) + (i+1) + 1
      connec[4*k+3] = n_vtx_seg*(j+1) + i     + 1
      k += 1

  return {
    "n_vtx"        : n_vtx,
    "n_elt"        : n_elt,
    "vtx_coord"    : vtx_coord,
    "vtx_ln_to_gn" : vtx_ln_to_gn,
    "connec"       : connec,
    "elt_ln_to_gn" : elt_ln_to_gn
    }


def run():
  try:
    from pycwp import pycwp
  except:
    if i_rank == 0:
      print("      Error : CWIPI module not found (update PYTHONPATH variable)")
      sys.exit(1)

  # Initialize MPI
  comm = MPI.COMM_WORLD

  is_code1 = (comm.rank%2 == 0)

  if is_code1:
    code_name         = "code1"
    coupled_code_name = "code2"
  else:
    code_name         = "code2"
    coupled_code_name = "code1"

  # Initialize CWIPI
  is_active_rank = True
  intra_comm = pycwp.init(comm,
                          [code_name],
                          is_active_rank)

  # Create coupling environment
  coupling_name = "python_new_api_empty_parts"

  cpl = pycwp.Coupling(code_name,
                       coupling_name,
                       coupled_code_name,
                       pycwp.INTERFACE_SURFACE,
                       pycwp.COMM_PAR_WITH_PART,
                       pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                       1,
                       pycwp.DYNAMIC_MESH_STATIC,
                       pycwp.TIME_EXCH_USER_CONTROLLED)

  cpl.visu_set(1,
               pycwp.VISU_FORMAT_ENSIGHT,
               "text")

  # Define interface mesh
  mesh = generate_mesh(intra_comm[0],
                       4 if is_code1 else 3)

  print(f"rank {comm.rank} has {mesh['n_vtx']} vtx and {mesh['n_elt']} elt", flush=True)

  cpl.mesh_interf_vtx_set(0,
                          mesh["vtx_coord"],
                          mesh["vtx_ln_to_gn"])

  block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_QUAD4)

  cpl.mesh_interf_block_std_set(0,
                                block_id,
                                mesh["connec"],
                                mesh["elt_ln_to_gn"])

  cpl.mesh_interf_finalize()


  # Define fields
  field_name = "field"

  if is_code1:
    # send_data = np.array(mesh["elt_ln_to_gn"], dtype=np.double)
    send_data = np.empty(mesh["n_elt"], dtype=np.double)
    if mesh["n_elt"] > 0:
      send_data[:] = mesh["elt_ln_to_gn"][:]

    field = cpl.field_create(field_name,
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_CELL_CENTER,
                             pycwp.FIELD_EXCH_SEND,
                             pycwp.STATUS_ON)
    field.data_set(0,
                   pycwp.FIELD_MAP_SOURCE,
                   send_data)
  else:
    recv_data = np.empty(mesh["n_vtx"], dtype=np.double)

    field = cpl.field_create(field_name,
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_NODE,
                             pycwp.FIELD_EXCH_RECV,
                             pycwp.STATUS_ON)
    field.data_set(0,
                   pycwp.FIELD_MAP_TARGET,
                   recv_data)

  pycwp.time_step_beg(code_name,
                      0.0)

  # Compute spatial interpolation weights
  cpl.spatial_interp_weights_compute()

  # Exchange interpolated fields
  if is_code1:
    field.issend()
  else:
    field.irecv()

  if is_code1:
    field.wait_issend()
  else:
    field.wait_irecv()


  pycwp.time_step_end(code_name)


  if comm.rank == 0:
    print("End :)", flush=True)


  # Finalize CWIPI
  pycwp.finalize()

  # Finalize MPI
  MPI.Finalize()



if __name__ == '__main__':
  run()

