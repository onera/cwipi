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
import Pypdm.Pypdm as PDM
import numpy as np
import sys
import argparse



def gen_mesh(comm, n_part, n, center, radius, part_method):
  dmn_capsule = PDM.sphere_surf_icosphere_gen_nodal(comm,
                                                    n,
                                                    center[0],
                                                    center[1],
                                                    center[2],
                                                    radius)

  mpart = PDM.MultiPart(1,
                        np.array([n_part]).astype(np.int32),
                        0,
                        part_method,
                        1,
                        np.ones(1).astype(np.double),
                        comm)

  mpart.multipart_register_dmesh_nodal(0, dmn_capsule)

  mpart.multipart_run_ppart()


  pvtx_coord    = [mpart.multipart_vtx_coord_get(i, 0)["np_vtx_coord"]       for i in range(n_part)]
  pvtx_ln_to_gn = [mpart.multipart_ln_to_gn_get(i, 0) ["np_entity_ln_to_gn"] for i in range(n_part)]

  pface_vtx_idx = []
  pface_vtx     = []
  for i in range(n_part):
    edges = mpart.multipart_connectivity_get(i, 0, PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)
    faces = mpart.multipart_connectivity_get(i, 0, PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)

    face_edge_idx = faces["np_entity1_entity2_idx"]
    face_edge     = faces["np_entity1_entity2"]

    face_vtx = PDM.compute_face_vtx_from_face_and_edge(face_edge_idx,
                                                       face_edge,
                                                       edge_vtx)
    pface_vtx_idx.append(face_edge_idx)
    pface_vtx.append(face_vtx)


  return {
    "pvtx_coord"     : pvtx_coord,
    "pvtx_ln_to_gn"  : pvtx_ln_to_gn,
    "pface_vtx_idx"  : pface_vtx_idx,
    "pface_vtx"      : pface_vtx,
    "pface_ln_to_gn" : pface_ln_to_gn
  }




def runTest():
  try:
    from pycwp import pycwp
  except:
    print("Error : CWIPI module not found (update PYTHONPATH variable)")
    sys.exit(1)

  print("file :", PDM.__file__)

  parser = argparse.ArgumentParser()

  parser.add_argument("-n",          "--n_subdiv",            default=1)
  parser.add_argument("-n_part1",    "--n_part1",             default=1)
  parser.add_argument("-n_part2",    "--n_part2",             default=1)
  parser.add_argument("-disjoint",   "--disjoint_comm",       action="store_true")
  parser.add_argument("-v",          "--verbose",             action="store_true")
  parser.add_argument("-swap_codes", "--swap_codes",          action="store_true")
  parser.add_argument("-algo",       "--spatial_interp_algo", default=pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE)


  args = parser.parse_args()

  n_subdiv            = int(args.n_subdiv)
  n_part1             = int(args.n_part1)
  n_part2             = int(args.n_part2)
  disjoint_comms      = args.disjoint_comm
  verbose             = args.verbose
  swap_codes          = args.swap_codes
  spatial_interp_algo = int(args.spatial_interp_algo)

  part_method = 3

  global f
  comm = MPI.COMM_WORLD

  i_rank = comm.rank
  n_rank = comm.size

  if i_rank == 0:
    print("\nSTART: python_new_api_surf.py")


  if disjoint_comms:
    has_code = [
      i_rank <  n_rank//2,
      i_rank >= n_rank//2,
    ]
  else:
    has_code = [
      i_rank <  (2*n_rank)//3,
      i_rank >= n_rank//3,
    ]

  if swap_codes:
    has_code = has_code[::-1]


  all_code_name = ["code%d" % (i+1) for i in range(2)]
  all_n_part    = [n_part1, n_part2]


  code_name         = []
  coupled_code_name = []
  n_part            = []
  for i in range(2):
    if has_code[i]:
      code_name.append(all_code_name[i])
      coupled_code_name.append(all_code_name[(i+1)%2])
      n_part.append(all_n_part[i])


  # OUTPUT
  f = None
  if verbose:
    srank = '{0}'.format(i_rank)
    f = open("python_new_api_surf" + srank.zfill(4) + ".txt",'w')
    pycwp.output_file_set(f)


  # Initialize CWIPI
  if verbose:
    f.write("pycwp.init:\n")
  n_code = len(code_name)
  is_active_rank = np.ones(n_code,  dtype=np.int32)
  time_init      = np.zeros(n_code, dtype=np.double)
  intra_comm = pycwp.init(comm,
                          n_code,
                          code_name,
                          is_active_rank,
                          time_init)

  comm.Barrier()
  if i_rank == 0:
    print("CWIPI Init OK")


  if verbose:
    f.write("intra_comm : {}\n".format(intra_comm))
    f.flush()

  # Create coupling
  cpl_name = "python_new_api_global_data"
  cpl = []
  for icode in range(n_code):
    if verbose:
      f.write("running {}\n".format(code_name[icode]))
      f.flush()
    cpl.append(pycwp.Coupling(code_name[icode],
                              cpl_name,
                              coupled_code_name[icode],
                              pycwp.INTERFACE_SURFACE,
                              pycwp.COMM_PAR_WITH_PART,
                              spatial_interp_algo,
                              n_part[icode],
                              pycwp.DYNAMIC_MESH_STATIC,
                              pycwp.TIME_EXCH_USER_CONTROLLED))

  comm.Barrier()
  if i_rank == 0:
    print("Create coupling OK")

  if verbose:
    f.write("comm : {}\n".format(comm))
    f.flush()
  for icode in range(n_code):
    if verbose:
      f.write("intra_comm[{}] : {}\n".format(icode, intra_comm[icode]))
      f.flush()
    i_rank_intra = intra_comm[icode].rank
    n_rank_intra = intra_comm[icode].size
    if verbose:
      f.write("code {} : i_rank_intra = {}, n_rank_intra = {}".format(icode, i_rank_intra, n_rank_intra))
      f.flush()

  # Define interface mesh
  mesh = []
  for icode in range(n_code):
    mesh.append(gen_mesh(comm,#intra_comm[icode],
                         n_part[icode],
                         n_subdiv,
                         [0., 0., 0.],
                         1.,
                         part_method))

    for ipart in range(n_part[icode]):
      cpl[icode].mesh_interf_vtx_set(ipart,
                                     mesh[icode]["pvtx_coord"]   [ipart],
                                     mesh[icode]["pvtx_ln_to_gn"][ipart])

      id_block = cpl[icode].mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

      cpl[icode].mesh_interf_f_poly_block_set(ipart,
                                              id_block,
                                              mesh[icode]["pface_vtx_idx"] [ipart],
                                              mesh[icode]["pface_vtx"]     [ipart],
                                              mesh[icode]["pface_ln_to_gn"][ipart])

    cpl[icode].mesh_interf_finalize()


  comm.Barrier()
  if i_rank == 0:
    print("Set mesh OK")


  # Create field
  field_name = "field"
  visu_status = pycwp.STATUS_ON
  stride      = 2
  for icode in range(n_code):

    if code_name[icode] == all_code_name[0]:
      exch_type = pycwp.FIELD_EXCH_SEND
      map_type  = pycwp.FIELD_MAP_SOURCE

      send_val = [
        np.array([(i+1)*g for g in gnum for i in range(stride)], dtype=np.double) \
        for gnum in mesh[icode]["pface_ln_to_gn"]
      ]
      field_val = send_val
    else:
      exch_type = pycwp.FIELD_EXCH_RECV
      map_type  = pycwp.FIELD_MAP_TARGET

      recv_val = [
        np.zeros(stride*len(gnum), dtype=np.double) \
        for gnum in mesh[icode]["pface_ln_to_gn"]
      ]
      field_val = recv_val


    cpl[icode].field_create(field_name,
                            pycwp.DOUBLE,
                            pycwp.FIELD_STORAGE_INTERLACED,
                            stride,
                            pycwp.DOF_LOCATION_CELL_CENTER,
                            exch_type,
                            visu_status)

    for ipart in range(n_part[icode]):
      cpl[icode].field_set(field_name,
                           ipart,
                           map_type,
                           field_val[ipart])


  comm.Barrier()
  if i_rank == 0:
    print("Create fields OK")


  for icode in range(n_code):
    cpl[icode].spatial_interp_property_set("tolerance", "double", "1e-2")

    cpl[icode].spatial_interp_weights_compute()


  comm.Barrier()
  if i_rank == 0:
    print("Interpolation weights computation OK")


  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      cpl[icode].field_issend(code_name[icode])
    else:
      cpl[icode].field_irecv(code_name[icode])


  error = False
  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      cpl[icode].field_wait_issend(code_name[icode])
    else:
      cpl[icode].field_wait_irecv(code_name[icode])

      # check received data
      for ipart in range(n_part[icode]):
        for i, g in enumerate(mesh[icode]["pface_ln_to_gn"][ipart]):
          if verbose:
            f.write("{} received {}\n".format(g, recv_val[stride*i:stride*(i+1)]))
            f.flush()
          for j in range(stride):
            if abs(recv_val[stride*i+j] - (j+1)*g > 1e-9):
              error = True
              print("field: {} received {} (expected {})".format(g, recv_val[stride*i+j], (j+1)*g))

  if error:
    exit(1)

  comm.Barrier()
  if i_rank == 0:
    print("Exchange interpolated field OK")

  # Part data
  if verbose:
    f.write("-- Part Data --\n")
    f.flush()
  part_data_name = "part_data"

  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      exch_type = pycwp.PARTDATA_SEND
    else:
      exch_type = pycwp.PARTDATA_RECV
      # reset recv_val
      for ipart in range(n_part[icode]):
        recv_val[ipart][:] = -1234.


    cpl[icode].part_data_create(part_data_name,
                                exch_type,
                                mesh[icode]["pface_ln_to_gn"])

  request = []
  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      request.append(cpl[icode].part_data_issend(part_data_name,
                                                 stride,
                                                 send_val))
    else:
      request.append(cpl[icode].part_data_irecv(part_data_name,
                                                stride,
                                                recv_val))

  error = False
  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      cpl[icode].part_data_wait_issend(part_data_name,
                                       request[icode])
    else:
      cpl[icode].part_data_wait_irecv(part_data_name,
                                      request[icode])

      # check received data
      for ipart in range(n_part[icode]):
        for i, g in enumerate(mesh[icode]["pface_ln_to_gn"][ipart]):
          if verbose:
            f.write("{} received {}\n".format(g, recv_val[stride*i:stride*(i+1)]))
            f.flush()
          for j in range(stride):
            if abs(recv_val[stride*i+j] - (j+1)*g > 1e-9):
              error = True
              print("part_data: {} received {} (expected {})".format(g, recv_val[stride*i+j], (j+1)*g))

  if error:
    exit(1)

  for icode in range(n_code):
    cpl[icode].part_data_del(part_data_name)

  comm.Barrier()
  if i_rank == 0:
    print("Part data exchange OK")


  # Global data
  if verbose:
    f.write("-- Global Data --\n")
    f.flush()
  global_data_name = "global_data"
  global_stride   = 2
  global_n_entity = 3
  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      send_global_data = np.array([(i+1)*(j+1) for i in range(global_n_entity) for j in range(global_stride)], \
                                  dtype=np.int32)
      cpl[icode].global_data_issend(global_data_name,
                                    global_stride,
                                    global_n_entity,
                                    send_global_data)
    else:
      recv_global_data = np.zeros(n_entity*stride, dtype=np.int32)
      cpl[icode].global_data_irecv(global_data_name,
                                   global_stride,
                                   global_n_entity,
                                   recv_global_data)

  error = False
  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      cpl[icode].global_data_wait_issend(global_data_name)
      if verbose:
        f.write("send_global_data = {}\n".format(send_global_data))
        f.flush()
    else:
      cpl[icode].global_data_wait_irecv(global_data_name)
      if verbose:
        f.write("recv_global_data = {}\n".format(recv_global_data))
        f.flush()
      for i in range(global_n_entity):
        for j in range(global_stride):
          if recv_global_data[global_stride*i+j] != (i+1)*(j+1):
            error = True
            print("global_data: {} received {} (expected {})".format(i, recv_global_data[global_stride*i+j], (i+1)*(j+1))

  if error:
    exit(1)

  comm.Barrier()
  if i_rank == 0:
    print("Global data exchange OK")
    print("End")

  for icode in range(n_code):
    cpl[icode].field_del(cpl_name)
    cpl[icode].mesh_interf_del()

  # FINALIZE
  pycwp.finalize()

  # END
  if verbose:
    f.write("\nEnd.\n")
    f.close()
  comm.Barrier()
  MPI.Finalize()

if __name__ == '__main__':
    runTest()
