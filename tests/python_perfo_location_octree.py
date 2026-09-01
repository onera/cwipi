#!/usr/bin/env python

import argparse
import ctypes
import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM
import sys
import time

from cwipi import cwipi
from pycwp import pycwp
from pycwp.pycwp import gnum_dtype

all_code_name = ["target", "source"]


def compute_elt_centers(connec_idx, connec, vtx_coord):
  """
  Compute barycenters of elements (used to evaluate cell-centered analytical fields)
  """
  n_elt = len(connec_idx)-1
  elt_center = np.zeros(3*n_elt, dtype=np.double)

  for i in range(n_elt):
    for j in connec[connec_idx[i]:connec_idx[i+1]]:
      elt_center[3*i:3*(i+1)] += vtx_coord[3*(j-1):3*j]
    elt_center[3*i:3*(i+1)] /= connec_idx[i+1] - connec_idx[i]

  return elt_center


def setup_mesh(comm, n_vtx_seg, elt_type, xyz_min, length, order, random_factor, i_surface, part_method):
  """
  Generate mesh
  """
  surface = i_surface >= 0

  dcube = PDM.DCubeNodalGenerator(nx     = n_vtx_seg,
                                  ny     = n_vtx_seg,
                                  nz     = n_vtx_seg,
                                  length = length,
                                  zero_x = xyz_min[0]*length,
                                  zero_y = xyz_min[1]*length,
                                  zero_z = xyz_min[2]*length,
                                  t_elmt = elt_type,
                                  order  = order,
                                  comm   = comm)

  dcube.set_random_factor(random_factor)
  dcube.compute()

  dmn = dcube.get_dmesh_nodal()

  n_part = 1
  mpart = PDM.MultiPart(1,
                        np.array([n_part]).astype(np.int32),
                        0,
                        part_method,
                        1,
                        np.ones(1).astype(np.double),
                        comm)

  mpart.dmesh_nodal_set(0, dmn)

  mpart.compute()


  # Deform
  vtx_coord    = mpart.vtx_coord_get(0, 0)
  vtx_ln_to_gn = mpart.ln_to_gn_get(0, 0, PDM._PDM_MESH_ENTITY_VTX)

  y0 = xyz_min[1] + 0.5
  z0 = xyz_min[2] + 0.5
  r2 = ((vtx_coord[1::3]/length - y0)**2 + (vtx_coord[2::3]/length - z0)**2)

  vtx_coord[0::3] += 0.4*length*r2

  if surface:
    # Extract surface mesh
    face_vtx_idx, face_vtx = mpart.connectivity_get(0, 0, PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)

    face_ln_to_gn = mpart.ln_to_gn_get(0, 0, PDM._PDM_MESH_ENTITY_FACE)
    cell_ln_to_gn = mpart.ln_to_gn_get(0, 0, PDM._PDM_MESH_ENTITY_CELL)

    groups = mpart.group_get(0, 0, PDM._PDM_MESH_ENTITY_FACE)

    extrp = PDM.ExtractPart(2,
                            n_part,
                            n_part,
                            PDM._PDM_EXTRACT_PART_KIND_LOCAL,
                            0,
                            1,
                            comm)

    list_faces = groups["group_entity"][groups["group_entity_idx"][i_surface]:groups["group_entity_idx"][i_surface+1]]
    extrp.selected_lnum_set(0, list_faces)

    extrp.part_set(0,
                   0,
                   len(face_ln_to_gn),
                   0,
                   len(vtx_ln_to_gn),
                   None,
                   None,
                   None,
                   None,
                   None,
                   face_vtx_idx,
                   face_vtx,
                   None,
                   face_ln_to_gn,
                   None,
                   vtx_ln_to_gn,
                   vtx_coord)

    extrp.compute()

    extrp_face_vtx_idx, extrp_face_vtx = extrp.connectivity_get(0, PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)

    extrp_vtx_coord = extrp.vtx_coord_get(0)

    return {
      "n_vtx"       : extrp.n_entity_get(0, PDM._PDM_MESH_ENTITY_VTX),
      "n_elt"       : extrp.n_entity_get(0, PDM._PDM_MESH_ENTITY_FACE),
      "coord"       : extrp_vtx_coord,
      "connec_idx"  : extrp_face_vtx_idx,
      "connec"      : extrp_face_vtx,
      "elt_center"  : compute_elt_centers(extrp_face_vtx_idx, extrp_face_vtx, extrp_vtx_coord),
      "elt_ln_to_gn": extrp.ln_to_gn_get(0, PDM._PDM_MESH_ENTITY_FACE),
      "vtx_ln_to_gn": extrp.ln_to_gn_get(0, PDM._PDM_MESH_ENTITY_VTX)
    }

  else:
    # Volume mesh
    pmn = mpart.part_mesh_nodal_get(0)

    section = pmn.get_sections(PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                               0)[0]

    connec = section["np_connec"]
    n_elt  = len(section["np_numabs"])

    strides = {
      PDM._PDM_MESH_NODAL_TETRA4  : 4,
      PDM._PDM_MESH_NODAL_PYRAMID5: 5,
      PDM._PDM_MESH_NODAL_PRISM6  : 6,
      PDM._PDM_MESH_NODAL_HEXA8   : 8
    }
    connec_idx = strides[elt_type] * np.arange(n_elt+1, dtype=np.int32)

    return {
      "n_vtx"       : len(vtx_ln_to_gn),
      "n_elt"       : n_elt,
      "coord"       : vtx_coord,
      "connec_idx"  : connec_idx,
      "connec"      : connec,
      "elt_center"  : compute_elt_centers(connec_idx, connec, vtx_coord),
      "elt_ln_to_gn": section["np_numabs"],
      "vtx_ln_to_gn": vtx_ln_to_gn
    }


def eval_field(coord):
  """
  Analytical field
  """
  return np.array(coord[::3])


def run_coupling():
  """
  Run one-way coupling source -> target
  """
  # Parse command line arguments
  parser = argparse.ArgumentParser()

  parser.add_argument("-n1",
                      "--n_vtx_seg1",
                      help="Number of vertices per side of the target domain. The total number of vertices is this number cubed (or squared in surface mode).",
                      type=int,
                      default=10)
  parser.add_argument("-n2",
                      "--n_vtx_seg2",
                      help="Number of vertices per side of the source domain. The total number of vertices is this number cubed (or squared in surface mode).",
                      type=int,
                      default=10)
  parser.add_argument("-n_ov",
                      "--n_overlap",
                      help="Overlap (0: Surface overlap, N > 0: Partial volume overlap over N cells, N < 0: Full volume overlap)",
                      type=int,
                      default=2)
  parser.add_argument("-part_method",
                      "--part_method",
                      help="Partitioning method (1: ParMETIS*, 2: PT-Scotch*, 3: Hilbert) (*: if available, see ParaDiGM install guide)",
                      type=int,
                      default=PDM._PDM_SPLIT_DUAL_WITH_HILBERT)
  parser.add_argument("-old",
                      "--old",
                      help="Use old API (FVM location algorithm)",
                      action="store_true")
  parser.add_argument("-inactive",
                      "--inactive",
                      help="Enable MPI communicator split to exclude ranks with empty partitions",
                      action="store_true")
  parser.add_argument("-v",
                      "--visu",
                      help="Enable output for visualization",
                      action="store_true")

  args = parser.parse_args()

  tolerance           = 1e-2 if args.n_overlap == 0 else 1e-3
  spatial_interp_algo = pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE
  elt_type            = PDM._PDM_MESH_NODAL_TETRA4

  comm = MPI.COMM_WORLD
  i_rank = comm.rank
  n_rank = comm.size

  # Setup code distribution on processes
  has_code = [i_rank < n_rank//2, False]
  has_code[1] = not has_code[0]

  all_n_vtx_seg = [args.n_vtx_seg1, args.n_vtx_seg2]

  code_name         = []
  coupled_code_name = []
  n_vtx_seg         = []
  mesh_comm         = []
  for i in range(2):
    if has_code[i]:
      code_name        .append(all_code_name[i])
      coupled_code_name.append(all_code_name[(i+1)%2])
      n_vtx_seg        .append(all_n_vtx_seg[i])

    mesh_comm.append(comm.Split(has_code[i], i_rank))
  n_code = len(code_name)

  mesh = []
  n_src = 0
  n_tgt = 0
  for icode in range(n_code):
    xyz_min = np.zeros(3, dtype=np.double)
    i_surface = -1
    if code_name[icode] == all_code_name[0]:
      if args.n_overlap == 0:
        xyz_min[0] = 1
        i_surface  = 2
      else:
        if args.n_overlap > 0:
          xyz_min[0] = 1. - args.n_overlap/(n_vtx_seg[icode] - 1)
    else:
      if args.n_overlap == 0:
        i_surface = 3

    mesh.append(setup_mesh(mesh_comm[icode],
                           n_vtx_seg[icode],
                           elt_type,
                           xyz_min,
                           20,  # length
                           1,   # order
                           0.8, # random_factor
                           i_surface,
                           args.part_method))

    if code_name[icode] == all_code_name[0]:
      n_tgt = mesh[-1]["n_elt"]
    else:
      n_src = mesh[-1]["n_elt"]

  # Initialize CWIPI
  is_active_rank = True
  if args.inactive:
    is_active_rank = (mesh[0]["n_elt"] > 0)

  if args.old:
    intra_comm = cwipi.init(comm, code_name[0])
    intra_comm = [intra_comm]
  else:
    intra_comm = pycwp.init(comm,
                            code_name,
                            is_active_rank)

  i_is_active_rank = int(is_active_rank)
  n_active_rank = comm.allreduce(i_is_active_rank, MPI.SUM)

  ln_elt = np.array([n_tgt, n_src], dtype=np.int64)
  gn_elt = np.empty(2, dtype=np.int64)
  comm.Allreduce(ln_elt, gn_elt, MPI.SUM)

  if i_rank == 0:
    print(f"n_active_rank = {n_active_rank}", flush=True)
    print(f"Global number of targets = {gn_elt[0]}", flush=True)
    print(f"Global number of sources = {gn_elt[1]}", flush=True)

  # Create coupling
  cpl_name = "python_perfo_location_octree"
  cpl = []

  if args.old:
    cpl.append(cwipi.Coupling(cpl_name,
                              cwipi.COUPLING_PARALLEL_WITH_PARTITIONING,
                              coupled_code_name[0],
                              2 if args.n_overlap == 0 else 3,
                              tolerance,
                              cwipi.STATIC_MESH,
                              cwipi.SOLVER_CELL_CENTER,
                              1 if args.visu else -1,
                              "Ensight",
                              "txt"))
  else:
    if is_active_rank:
      for icode in range(n_code):
        cpl.append(pycwp.Coupling(code_name[icode],
                                  cpl_name,
                                  coupled_code_name[icode],
                                  pycwp.INTERFACE_SURFACE if args.n_overlap == 0 else pycwp.INTERFACE_VOLUME,
                                  pycwp.COMM_PAR_WITH_PART,
                                  spatial_interp_algo,
                                  1,
                                  pycwp.DYNAMIC_MESH_STATIC,
                                  pycwp.TIME_EXCH_USER_CONTROLLED))

      if args.visu:
        for icode in range(n_code):
          cpl[icode].visu_set(1,
                              pycwp.VISU_FORMAT_ENSIGHT,
                              "binary")


  # Define interface mesh
  for icode in range(n_code):
    if args.old:
      cpl[icode].define_mesh(mesh[icode]["n_vtx"],
                             mesh[icode]["n_elt"],
                             mesh[icode]["coord"],
                             mesh[icode]["connec_idx"],
                             mesh[icode]["connec"])
    else:
      if is_active_rank:
        cpl[icode].mesh_interf_vtx_set(0,
                                       mesh[icode]["coord"],
                                       mesh[icode]["vtx_ln_to_gn"])

        if args.n_overlap == 0:
          id_block = cpl[icode].mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)
          cpl[icode].mesh_interf_f_poly_block_set(0,
                                                  id_block,
                                                  mesh[icode]["connec_idx"],
                                                  mesh[icode]["connec"],
                                                  mesh[icode]["elt_ln_to_gn"])
        else:
          id_block = cpl[icode].mesh_interf_block_add(elt_type)

          cpl[icode].mesh_interf_block_std_set(0,
                                               id_block,
                                               mesh[icode]["connec"],
                                               mesh[icode]["elt_ln_to_gn"])
        cpl[icode].mesh_interf_finalize()

  # Field
  field = []
  field_name = "field"
  visu_status = pycwp.STATUS_ON
  dof_location = pycwp.DOF_LOCATION_CELL_CENTER
  for icode in range(n_code):

    if code_name[icode] == all_code_name[0]:
      exch_type = pycwp.FIELD_EXCH_RECV
      map_type  = pycwp.FIELD_MAP_TARGET

      recv_val = np.empty(mesh[icode]["n_elt"])
      field_val = recv_val
    else:
      exch_type = pycwp.FIELD_EXCH_SEND
      map_type  = pycwp.FIELD_MAP_SOURCE

      send_val  = eval_field(mesh[icode]["elt_center"])
      field_val = send_val

    if args.old:
      if code_name[icode] == all_code_name[1]:
        mesh[icode]["dummy_pts"] = np.empty(0, dtype=np.double)
        cpl[icode].set_points_to_locate(0, mesh[icode]["dummy_pts"])
    else:
      if is_active_rank:
        field.append(cpl[icode].field_create(field_name,
                                             pycwp.DOUBLE,
                                             pycwp.FIELD_STORAGE_INTERLACED,
                                             1,
                                             dof_location,
                                             exch_type,
                                             visu_status))

        field[icode].data_set(0,
                              map_type,
                              field_val)

  if not args.old and is_active_rank:
    for icode in range(n_code):
      pycwp.time_step_beg(code_name[icode], 0.0)

  # Spatial interpolation
  comm.Barrier()
  tstart = MPI.Wtime()
  if args.old:
    cpl[0].locate()
  else:
    if is_active_rank:
      for icode in range(n_code):
        if spatial_interp_algo != pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT:
          cpl[icode].spatial_interp_property_set("tolerance",
                                                 pycwp.DOUBLE,
                                                 "%f" % tolerance)

        cpl[icode].spatial_interp_weights_compute()

  tend = MPI.Wtime()

  comm.Barrier()

  delta_t = np.array(tend - tstart, dtype=float)
  min_delta_t = np.empty(1, dtype=float)
  max_delta_t = np.empty(1, dtype=float)
  comm.Allreduce([delta_t, MPI.FLOAT], [min_delta_t, MPI.FLOAT], MPI.MIN)
  comm.Allreduce([delta_t, MPI.FLOAT], [max_delta_t, MPI.FLOAT], MPI.MAX)

  if i_rank == 0:
    print(f"Location time: {min_delta_t[0]} {max_delta_t[0]}")

  # Exchange
  tstart = time.time()
  if args.old:
    if code_name[0] == all_code_name[0]:
      result = cpl[0].exchange(field_name,
                               1,
                               1,
                               0.1,
                               field_name + "_s", None,
                               field_name + "_r", recv_val)
    else:
      result = cpl[0].exchange(field_name,
                               1,
                               1,
                               0.1,
                               field_name + "_s", send_val,
                               field_name + "_r", None)
  else:
    if is_active_rank:
      for icode in range(n_code):
        if code_name[icode] == all_code_name[0]:
          field[icode].irecv()
        else:
          field[icode].issend()

      for icode in range(n_code):
        if code_name[icode] == all_code_name[0]:
          field[icode].wait_irecv()
        else:
          field[icode].wait_issend()
  tend = time.time()

  comm.Barrier()

  delta_t = np.array(tend - tstart, dtype=float)
  min_delta_t = np.empty(1, dtype=float)
  max_delta_t = np.empty(1, dtype=float)
  comm.Allreduce([delta_t, MPI.FLOAT], [min_delta_t, MPI.FLOAT], MPI.MIN)
  comm.Allreduce([delta_t, MPI.FLOAT], [max_delta_t, MPI.FLOAT], MPI.MAX)

  if i_rank == 0:
    print(f"Exchange time: {min_delta_t[0]} {max_delta_t[0]}")

  if not args.old and is_active_rank:
    for icode in range(n_code):
      pycwp.time_step_end(code_name[icode])


  # Count uncomputed tgt
  n_computed_tgts   = 0
  n_uncomputed_tgts = 0
  n_involved_srcs   = 0
  for icode in range(n_code):
    if code_name[icode] == all_code_name[0]:
      if args.old:
        n_uncomputed_tgts = cpl[0].get_n_not_located_points()
        n_computed_tgts   = cpl[0].get_n_located_points()
      else:
        if is_active_rank:
          n_uncomputed_tgts = field[icode].n_uncomputed_tgts_get(0)
          n_computed_tgts   = field[icode].n_computed_tgts_get(0)
    else:
      if not args.old and is_active_rank:
        n_involved_srcs = field[icode].n_involved_srcs_get(0)

  # Count global number of unlocated pts (use part-to-block to be more accurate?)
  n_uncomputed = np.array(n_uncomputed_tgts, dtype=np.int64)
  gn_uncomputed = np.empty(1, dtype=np.int64)
  comm.Allreduce([n_uncomputed, MPI.LONG], [gn_uncomputed, MPI.LONG], MPI.SUM)

  n_computed = np.array(n_computed_tgts, dtype=np.int64)
  gn_computed = np.empty(1, dtype=np.int64)
  comm.Allreduce([n_computed, MPI.LONG], [gn_computed, MPI.LONG], MPI.SUM)

  n_involved = np.array(n_involved_srcs, dtype=np.int64)
  gn_involved = np.empty(1, dtype=np.int64)
  comm.Allreduce([n_involved, MPI.LONG], [gn_involved, MPI.LONG], MPI.SUM)


  if i_rank == 0:
    print(f"Global nb of uncomputed tgt: {gn_uncomputed[0]}")
    print(f"Global nb of computed tgt  : {gn_computed[0]}")
    print(f"Global nb of involved src  : {gn_involved[0]}")

  # Finalize
  if args.old:
    cwipi.finalize()
  else:
    del field
    cpl[0].mesh_interf_del()
    del cpl[0]
    pycwp.finalize()

  comm.Barrier()

  if i_rank == 0:
    print("The End :D\n")

  # Finalize MPI
  MPI.Finalize()



if __name__ == '__main__':
  run_coupling()
