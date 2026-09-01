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

def grid(n, z=0):
  coord = np.empty(3*n**2, dtype=np.double)
  step = 1./float(n-1)
  k = 0
  for j in range(n):
    for i in range(n):
      coord[3*k  ] = i*step
      coord[3*k+1] = j*step
      coord[3*k+2] = z
      k += 1

  connec_idx = 4*np.arange((n-1)**2 + 1, dtype=np.int32)
  connec = np.empty(4*(n-1)**2, dtype=np.int32)
  k = 0
  for j in range(n-1):
      for i in range(n-1):
        connec[4*k  ] = n*j + i + 1
        connec[4*k+1] = connec[4*k] + 1
        connec[4*k+2] = connec[4*k] + n + 1
        connec[4*k+3] = connec[4*k] + n
        k += 1

  return coord, connec_idx, connec

def run():
  # Initialize MPI
  comm = MPI.COMM_WORLD
  i_rank = comm.rank

  # Load Python CWIPI module
  try:
    from pycwp import pycwp
  except:
    if i_rank == 0:
      print("      Error : CWIPI module not found (update PYTHONPATH variable)")
      print(f"cwp : {pycwp.__file__}")
      sys.exit(1)

  # Even ranks run code0, odd ranks run code1
  i_code = i_rank % 2
  code_name = f"code{i_code}"

  # Coupled with the other code
  coupled_code_name = f"code{(i_code+1)%2}"

  # Initialize CWIPI
  is_active_rank = True

  intra_comm = pycwp.init(comm,
                          [code_name],
                          is_active_rank)

  # Create coupling
  print(f"I am rank {i_rank} and I run {code_name}, I am coupled with {coupled_code_name}", flush=True)
  cpl = pycwp.Coupling(code_name,
                       "python_new_api_conditioned_exchange",
                       coupled_code_name,
                       pycwp.INTERFACE_SURFACE,
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
  coord, connec_idx, connec = grid(3+i_code)

  cpl.mesh_interf_vtx_set(0,
                          coord,
                          None)

  block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

  cpl.mesh_interf_f_poly_block_set(0,
                                   block_id,
                                   connec_idx,
                                   connec,
                                   None)

  cpl.mesh_interf_finalize()

  # Create control parameter to trigger conditional field exchanges
  # if i_code == 0:
  pycwp.param_lock(code_name)
  pycwp.param_add_int(code_name, "do_exchange_field", 0)
  pycwp.param_unlock(code_name)

  # Define field
  n_vtx = len(coord)//3
  if i_code == 0:
    field = cpl.field_create("field",
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_NODE,
                             pycwp.FIELD_EXCH_SEND,
                             pycwp.STATUS_ON)

    send_data    = np.empty(n_vtx, dtype=np.double)
    send_data[:] = coord[0::3]

    field.data_set(0,
                   pycwp.FIELD_MAP_SOURCE,
                   send_data)
  else:
    field = cpl.field_create("field",
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_NODE,
                             pycwp.FIELD_EXCH_RECV,
                             pycwp.STATUS_ON)

    recv_data = np.empty(n_vtx, dtype=np.double)

    field.data_set(0,
                   pycwp.FIELD_MAP_TARGET,
                   recv_data)

  # Time loop
  n_time_steps = 10
  time = 0.0

  # Begin empty time step
  pycwp.time_step_beg(code_name,
                      time)
  pycwp.time_step_end(code_name)
  time += 1.0

  for step in range(n_time_steps):
    print(f"begin step {step}", flush=True)

    # Begin time step
    pycwp.time_step_beg(code_name,
                        time)

    # Compute spatial interpolation weights
    if step == 0:
      cpl.spatial_interp_weights_compute()

    cpl.barrier() # Seems mandatory
    # Update field values in code0
    if i_code == 0:
      send_data += 1

    if i_code == 0:
      # Update 'do_exchange_field' and trigger send if condition is satisfied
      if intra_comm[0].rank == 0:
        r = np.random.rand()
      else:
        r = None
      r = intra_comm[0].bcast(r, root=0)

      do_exchange_field_code0 = int(r > -1)
      print(f"do_exchange_field_code0 = {do_exchange_field_code0}")
      pycwp.param_lock("code0")
      pycwp.param_set_int("code0", "do_exchange_field", do_exchange_field_code0)
      pycwp.param_unlock("code0")

    cpl.barrier() # Seems mandatory

    if i_code == 0:
      if do_exchange_field_code0:
        field.issend()

    if i_code == 1:
      # Check condition for exchange
      do_exchange_field_code1 = pycwp.param_get("code0", "do_exchange_field", pycwp.INT)
      print(f"do_exchange_field_code1 = {do_exchange_field_code1}")

      if do_exchange_field_code1:
        field.irecv()

    cpl.barrier() # Seems mandatory

    if i_code == 0:
      if do_exchange_field_code0:
        field.wait_issend()
    else:
      if do_exchange_field_code1:
        field.wait_irecv()

    # End time step
    print(f"end step {step}", flush=True)
    pycwp.time_step_end(code_name)

    cpl.barrier() # Seems mandatory

    time += 1.0

  # Finalize CWIPI
  del field
  cpl.mesh_interf_del()
  del cpl
  pycwp.finalize()

if __name__ == '__main__':
  run()
