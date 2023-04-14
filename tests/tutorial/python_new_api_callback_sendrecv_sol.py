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

def my_interpolation(local_code_name,
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

  print("in a user interpolation function !")

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
      # TO DO : some tgt have zero src???
      for isrc in range(tgt_to_src_idx[itgt], tgt_to_src_idx[itgt+1]):
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



def run_coupling():
  # Initialize MPI
  comm = MPI.COMM_WORLD
  i_rank = comm.rank
  n_rank = comm.size

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

  print(f"Python : {i_rank}/{n_rank} je suis là")


  # Initialize CWIPI
  n_code = 1

  code_name = ["codePython"]

  is_active_rank = np.ones (1, dtype=np.int32)
  time_init      = np.zeros(1, dtype=np.double)

  intra_comm = pycwp.init(comm,
                          n_code,
                          code_name,
                          is_active_rank,
                          time_init)

  print(f"Python : {i_rank}/{n_rank} CWP_Init OK")

  n_part = 1

  # Generate mesh
  mesh = Pypdm.generate_mesh_rectangle_simplified(intra_comm[0],
                                                  20)
  print(f"Python : {i_rank}/{n_rank} generate_mesh_rectangle_simplified OK")

  # Create first coupling C <-> Python
  cpl_CP = pycwp.Coupling(code_name[0],
                          "coupling_C_Python",
                          "codeC",
                          pycwp.INTERFACE_SURFACE,
                          pycwp.COMM_PAR_WITH_PART,
                          pycwp.SPATIAL_INTERP_FROM_CLOSEST_SOURCES_LEAST_SQUARES,
                          n_part,
                          pycwp.DYNAMIC_MESH_STATIC,
                          pycwp.TIME_EXCH_USER_CONTROLLED)

  print(f"Python : {i_rank}/{n_rank} Coupling OK")


  # Set coupling visualisation
  cpl_CP.visu_set(1,
                  pycwp.VISU_FORMAT_ENSIGHT,
                  "text")

  # Set the mesh vertices
  cpl_CP.mesh_interf_vtx_set(0,
                             mesh["n_vtx"],
                             mesh["coords"],
                             None)

  # Set the mesh elements
  block_id = cpl_CP.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

  cpl_CP.mesh_interf_f_poly_block_set(0,
                                      block_id,
                                      mesh["n_elt"],
                                      mesh["elt_vtx_idx"],
                                      mesh["elt_vtx"],
                                      None)

  # Finalize mesh
  cpl_CP.mesh_interf_finalize()

  print(f"Python : {i_rank}/{n_rank} mesh_interf_finalize OK")

  # Define field
  field_name = "coord_x"

  cpl_CP.field_create(field_name,
                      pycwp.DOUBLE,
                      pycwp.FIELD_STORAGE_INTERLACED,
                      1,
                      pycwp.DOF_LOCATION_NODE,
                      pycwp.FIELD_EXCH_SENDRECV,
                      pycwp.STATUS_ON)


  send_field_data = mesh["coords"][0::3]
  recv_field_data = np.zeros(mesh["n_vtx"], dtype=np.double)

  cpl_CP.field_set(field_name,
                   0,
                   pycwp.FIELD_MAP_SOURCE,
                   send_field_data)

  cpl_CP.field_set(field_name,
                   0,
                   pycwp.FIELD_MAP_TARGET,
                   recv_field_data)

  # Set user-defined interpolation function
  cpl_CP.interp_function_set(field_name,
                             my_interpolation)


  # Spatial interpolation
  cpl_CP.spatial_interp_property_set("n_closest_pts",
                                     "int",
                                     "3")

  cpl_CP.spatial_interp_weights_compute()

  # Exchange interpolated fields
  cpl_CP.field_issend(field_name)
  cpl_CP.field_irecv (field_name)

  cpl_CP.field_wait_issend(field_name)
  cpl_CP.field_wait_irecv (field_name)

  # Delete Mesh
  cpl_CP.mesh_interf_del()




  # # Create second coupling Python <-> Fortran
  # cpl_PF = pycwp.Coupling(code_name[0],
  #                         "coupling_Python_Fortran",
  #                         ["codeC"],
  #                         pycwp.INTERFACE_SURFACE,
  #                         pycwp.COMM_PAR_WITH_PART,
  #                         pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
  #                         n_part,
  #                         pycwp.DYNAMIC_MESH_STATIC,
  #                         pycwp.TIME_EXCH_USER_CONTROLLED)


  # # Set coupling visualisation
  # cpl_PF.visu_set(1,
  #                 pycwp.VISU_FORMAT_ENSIGHT,
  #                 "text")

  # # Set the mesh vertices
  # cpl_PF.mesh_interf_vtx_set(0,
  #                            mesh["n_vtx"],
  #                            mesh["coords"],
  #                            None)

  # # Set the mesh elements
  # block_id = cpl_PF.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

  # cpl_PF.mesh_interf_f_poly_block_set(0,
  #                                     block_id,
  #                                     mesh["n_elt"],
  #                                     mesh["elt_vtx_idx"],
  #                                     mesh["elt_vtx"],
  #                                     None)

  # # Finalize mesh
  # cpl_PF.mesh_interf_finalize()


  # # Define field
  # field_name = "coord_x"

  # cpl_PF.field_create(field_name,
  #                     pycwp.DOUBLE,
  #                     pycwp.FIELD_STORAGE_INTERLACED,
  #                     1,
  #                     pycwp.DOF_LOCATION_NODE,
  #                     pycwp.FIELD_EXCH_SENDRECV,
  #                     pycwp.STATUS_ON)

  # cpl_PF.field_set(field_name,
  #                  0,
  #                  pycwp.FIELD_MAP_SOURCE,
  #                  send_field_data)

  # cpl_PF.field_set(field_name,
  #                  0,
  #                  pycwp.FIELD_MAP_TARGET,
  #                  recv_field_data)

  # # Set user-defined interpolation function
  # cpl_PF.interp_function_set(field_name,
  #                            my_interpolation)


  # # Spatial interpolation
  # cpl_PF.spatial_interp_property_set("tolerance",
  #                                    "double",
  #                                    "0.1")

  # cpl_PF.spatial_interp_weights_compute()

  # # Exchange interpolated fields
  # cpl_PF.field_issend(field_name)
  # cpl_PF.field_irecv (field_name)

  # cpl_PF.field_wait_issend(field_name)
  # cpl_PF.field_wait_irecv (field_name)

  # # Delete Mesh
  # cpl_PF.mesh_interf_del()

  # Finalize CWIPI
  pycwp.finalize()

  print(f"Python rank {i_rank} FINISHED :D")

  # Finalize MPI
  MPI.Finalize()


if __name__ == '__main__':
  print("C'est parti!")
  run_coupling()

