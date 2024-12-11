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


  # Create coupling
  cpl = pycwp.Coupling(code_name,
                       "python_new_api_vol_from_cell_vtx",
                       coupled_code_name,
                       pycwp.INTERFACE_VOLUME,
                       pycwp.COMM_PAR_WITH_PART,
                       pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                       1,
                       pycwp.DYNAMIC_MESH_STATIC,
                       pycwp.TIME_EXCH_USER_CONTROLLED)

  # Set visu status
  cpl.visu_set(1,
               pycwp.VISU_FORMAT_ENSIGHT,
               "text")

  # Define interface mesh
  vtx_coord = np.empty(3*27, dtype=np.double)
  i_vtx = 0
  for k in range(3):
    for j in range(3):
      for i in range(3):
        vtx_coord[3*i_vtx  ] = i - 1
        vtx_coord[3*i_vtx+1] = j - 1
        vtx_coord[3*i_vtx+2] = k - 1
        i_vtx += 1
    
  if i_code == 0:
    # Rotate 90° around z-axis
    x = np.array(vtx_coord[0::3])
    vtx_coord[0::3] = -vtx_coord[1::3]
    vtx_coord[1::3] = x

  cell_vtx_idx = np.array([0, 8, 16, 22, 28, 34, 40, 45, 50, 55, 60, 65, 70, 74, 78, 82, 86, 90, 94, 98, 102, 106, 110], dtype=np.int32)
  cell_vtx = np.array([
    1, 2, 5, 4, 10, 11, 14, 13,
    4, 5, 8, 7, 13, 14, 17, 16,
    2, 3, 5, 11, 12, 14,
    3, 6, 5, 12, 15, 14,
    5, 6, 9, 14, 15, 18,
    5, 9, 8, 14, 18, 17,
    10, 11, 14, 13, 20,
    13, 14, 23, 22, 20,
    10, 13, 22, 19, 20,
    13, 14, 17, 16, 26,
    13, 22, 23, 14, 26,
    13, 16, 25, 22, 26,
    11, 12, 14, 20,
    14, 23, 20, 24,
    14, 15, 24, 12,
    14, 12, 24, 20,
    12, 24, 20, 21,
    14, 15, 18, 24,
    14, 23, 24, 26,
    14, 18, 17, 26,
    14, 18, 26, 24,
    18, 26, 24, 27
  ], dtype=np.int32)

  cpl.mesh_interf_vtx_set(0,
                          vtx_coord,
                          None)

  cpl.mesh_interf_from_cellvtx_set(0,
                                   cell_vtx_idx,
                                   cell_vtx,
                                   None)

  cpl.mesh_interf_finalize()
  

  # Define field
  n_vtx = vtx_coord.size//3

  send_val = np.array(vtx_coord[0::3] + vtx_coord[1::3] + vtx_coord[2::3], dtype=np.double)
  recv_val = np.empty(n_vtx, dtype=np.double)

  field = cpl.field_create("field",
                           pycwp.DOUBLE,
                           pycwp.FIELD_STORAGE_INTERLACED,
                           1,
                           pycwp.DOF_LOCATION_NODE,
                           pycwp.FIELD_EXCH_SENDRECV,
                           pycwp.STATUS_ON)

  pycwp.time_step_beg(code_name, 0.)

  
  field.data_set(0,
                 pycwp.FIELD_MAP_SOURCE,
                 send_val)

  
  field.data_set(0,
                 pycwp.FIELD_MAP_TARGET,
                 recv_val)


  # Compute spatial interpolation weights
  cpl.spatial_interp_weights_compute()

  # Exchange interpolated field
  field.issend()
  field.irecv ()

  field.wait_issend()
  field.wait_irecv ()

  pycwp.time_step_end(code_name)

  # Check interpolated field
  max_err = max(np.absolute(send_val - recv_val))

  if max_err > 1e-14:
    print(f"Interpolation error is too high {max_err}")
    sys.exit(1)

  # Finalize CWIPI
  pycwp.finalize()

  if comm.rank == 0:
    print("The End")


if __name__ == '__main__':
  run_test()

