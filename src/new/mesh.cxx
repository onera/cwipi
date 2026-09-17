/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2026  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "mesh.hxx"

#include "cwp.h"

#include "block.hxx"
#include "blockCP.hxx"
#include "blockFP.hxx"
#include "blockHO.hxx"
#include "blockStd.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "factory.hpp"

#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_mesh_nodal.h"

#include <mpi.h>

#include <cstdlib>
#include <vector>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \typedef FB
   *
   * \brief Block Factory.
   *
   * A block \ref Factory which makes \ref Block class objects.
   * The type of block objects built depends on the block type \ref CWP_Block_t.
   */
  typedef Factory<Block, CWP_Block_t> FB;


  static void
  _gen_gnum
  (
    const PDM_MPI_Comm                comm,
    const int                         n_part,
    const std::vector<int>           &n_entity,
    const std::vector<double *>      &coords,
          std::vector<CWP_g_num_t *> &global_num
  )
  {
    PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,
                                               n_part,
                                               PDM_FALSE,
                                               1e-3,
                                               comm,
                                               PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_gnum_set_from_coords(gen_gnum,
                               i_part,
                               n_entity[i_part],
                               coords  [i_part],
                               nullptr);
    }

    PDM_gnum_compute(gen_gnum);

    for (int i_part = 0; i_part < n_part; i_part++) {
      global_num[i_part] = const_cast<CWP_g_num_t *>(PDM_gnum_get(gen_gnum,
                                                                  i_part));
    }

    PDM_gnum_free(gen_gnum);
  }


  static void
  _compute_centers_with_connectivity
  (
    const bool    with_poly,
    const bool    with_parent,
    const int     n_cell,
    const int     cell_n_vtx,
    const int    *cell_vtx_idx,
    const int    *cell_vtx,
    const int    *cell_face_idx,
    const int    *cell_face,
    const int     face_n_vtx,
    const int    *face_vtx_idx,
    const int    *face_vtx,
    const double *coords,
    const int    *parent_num,
          double *cell_centers
  )
  {
    double normalization = 0.;
    int    idx_write     = 0;

    if (with_poly) {
      for (int i_cell = 0; i_cell < n_cell; i_cell++) {
        idx_write = with_parent ? parent_num[i_cell] : i_cell;
        normalization = 0.;
        for (int i_face = cell_face_idx[i_cell]; i_face < cell_face_idx[i_cell+1]; i_face++) {
          int id_face = std::abs(cell_face[i_face]) - 1;
          int vtx_0   = (face_n_vtx < 0) ? face_vtx_idx[id_face] : id_face*face_n_vtx;
          int n_vtx   = (face_n_vtx < 0) ? face_vtx_idx[id_face+1] - face_vtx_idx[id_face] : face_n_vtx;
          for (int i_vtx = vtx_0; i_vtx < vtx_0+n_vtx; i_vtx++) {
            int id_vtx = face_vtx[i_vtx] - 1;
            for (int j = 0; j < 3; j++) {
              cell_centers[3*idx_write+j] += coords[3*id_vtx+j];
            }
          }
          normalization += (double) n_vtx;
        }
        for (int j = 0; j < 3; j++) {
          cell_centers[3*idx_write+j] /= normalization;
        }
      }
    }
    else {
      for (int i_cell = 0; i_cell < n_cell; i_cell++) {
        idx_write = with_parent ? parent_num[i_cell] : i_cell;
        int vtx_0 = (cell_n_vtx < 0) ? cell_vtx_idx[i_cell] : i_cell*cell_n_vtx;
        int n_vtx = (cell_n_vtx < 0) ? cell_vtx_idx[i_cell+1] - cell_vtx_idx[i_cell] : cell_n_vtx;
        for (int i_vtx = vtx_0; i_vtx < vtx_0+n_vtx; i_vtx++) {
          int id_vtx = cell_vtx[i_vtx] - 1;
          for (int j = 0; j < 3; j++) {
            cell_centers[3*idx_write+j] += coords[3*id_vtx+j];
          }
        }
        normalization = 1. / (double) n_vtx;
        for (int j = 0; j < 3; j++) {
          cell_centers[3*idx_write+j] *= normalization;
        }
      }
    }
  }


  static void
  _gen_gnum_with_connectivity
  (
    const bool                        with_poly,
    const PDM_MPI_Comm                comm,
    const int                         n_part,
    const std::vector<int>           &n_cell,
    const int                         cell_n_vtx,
    const std::vector<int *>         &cell_vtx_idx,
    const std::vector<int *>         &cell_vtx,
    const std::vector<int *>         &cell_face_idx,
    const std::vector<int *>         &cell_face,
    const int                         face_n_vtx,
    const std::vector<int *>         &face_vtx_idx,
    const std::vector<int *>         &face_vtx,
    const std::vector<double *>      &coords,
          std::vector<CWP_g_num_t *> &global_num
  )
  {
    std::vector<double*> cell_centers;
    cell_centers.resize(n_part, nullptr);

    for (int i_part = 0; i_part < n_part; i_part++) {
      cell_centers[i_part] = PDM_array_const_double(n_cell[i_part] * 3, 0.);

      _compute_centers_with_connectivity(with_poly,
                                         false,
                                         n_cell       [i_part],
                                         cell_n_vtx,
                                         cell_vtx_idx [i_part],
                                         cell_vtx     [i_part],
                                         cell_face_idx[i_part],
                                         cell_face    [i_part],
                                         face_n_vtx,
                                         face_vtx_idx [i_part],
                                         face_vtx     [i_part],
                                         coords       [i_part],
                                         nullptr,
                                         cell_centers [i_part]);

    }

    _gen_gnum(comm,
              n_part,
              n_cell,
              cell_centers,
              global_num);

    for (int i_part = 0; i_part < n_part; i_part++) {
      free(cell_centers[i_part]);
    }
  }


  Mesh::Mesh
  (
    const MPI_Comm           &local_comm,
          int                 n_part,
          CWP_Dynamic_mesh_t  displacement,
          Coupling           *cpl
  )
  : _n_blocks(0),
    _displacement(displacement),
    _cpl(cpl),
    _pmn(nullptr),
    _is_face_edge_method(false),
    _is_face_vtx_method(false),
    _is_cell_face_method(false),
    _is_cell_vtx_method(false),
    _is_vtx_gnum_computed(false),
    _is_elt_gnum_computed(false),
    _is_face_gnum_computed(false)
  {

    _local_comm_pdm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm *>(&local_comm));

    _n_part = n_part;

    /* Vertex data */
    _n_vtx         .resize(_n_part, 0);
    _coords        .resize(_n_part, nullptr);
    _global_num_vtx.resize(_n_part, nullptr);

    /* Element data */
    _n_elts     .resize(_n_part, 0);
    _elt_centers.resize(_n_part, nullptr);

    /* Cell data */
    _n_cells         .resize(_n_part, 0);
    _cell_face_idx   .resize(_n_part, nullptr);
    _cell_face       .resize(_n_part, nullptr);
    _cell_vtx_idx    .resize(_n_part, nullptr);
    _cell_vtx        .resize(_n_part, nullptr);
    _global_num_cells.resize(_n_part, nullptr);

    /* Face data */
    _n_faces         .resize(_n_part, 0);
    _face_edge_idx   .resize(_n_part, nullptr);
    _face_edge       .resize(_n_part, nullptr);
    _face_vtx_idx    .resize(_n_part, nullptr);
    _face_vtx        .resize(_n_part, nullptr);
    _global_num_faces.resize(_n_part, nullptr);

    /* Edge data */
    _n_edges .resize(_n_part, 0);
    _edge_vtx.resize(_n_part, nullptr);

    switch (_cpl->entitiesDimGet()) {
      case CWP_INTERFACE_POINT:
        _geom_kind = PDM_GEOMETRY_KIND_CORNER;
        break;
      case CWP_INTERFACE_LINEAR:
        _geom_kind = PDM_GEOMETRY_KIND_RIDGE;
        break;
      case CWP_INTERFACE_SURFACE:
        _geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
        break;
      case CWP_INTERFACE_VOLUME:
        _geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "Invalid entity dimension %d\n", _cpl->entitiesDimGet());
    }
  }


  Mesh::~Mesh()
  {
    meshDel();
  }


  void
  Mesh::elementCentersCompute
  (
    int i_part
  )
  {
    if (_elt_centers[i_part] == nullptr) {
      _elt_centers[i_part] = (double*) malloc(sizeof(double) * 3 * nEltsGet(i_part));
    }

    for (int i_block = 0; i_block < _n_blocks; i_block++) {
      double *elt_centers_block = _block_DB[i_block]->elementCentersGet(i_part);
      int    *parent_num_block  = _block_DB[i_block]->parentNumGet     (i_part);
      for (int i = 0; i < _block_DB[i_block]->nEltsGet(i_part); i++) {
        for (int j = 0; j < 3; j++) {
          _elt_centers[i_part][3*parent_num_block[i]+j] = elt_centers_block[3*i+j];
        }
      }
    }
  }


  double *
  Mesh::elementCentersGet
  (
    int i_part
  )
  {
    if (_elt_centers[i_part] == nullptr || _displacement != CWP_DYNAMIC_MESH_STATIC) {
      elementCentersCompute(i_part);
    }
    return _elt_centers[i_part];
  }


  CWP_g_num_t *
  Mesh::elementGlobalNumGet
  (
    int i_part
  )
  {
    return PDM_part_mesh_nodal_g_num_get_from_part(_pmn,
                                                   _geom_kind,
                                                   i_part,
                                                   PDM_OWNERSHIP_KEEP);
  }


  void
  Mesh::coordSet
  (
    const int         i_part,
    const int         n_vtx,
          double      coords[],
          CWP_g_num_t global_num[]
  )
  {
    _n_vtx         [i_part] = n_vtx;
    _coords        [i_part] = coords;
    _global_num_vtx[i_part] = global_num;
  }


  void
  Mesh::blockDBUpdate()
  {
    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(_pmn,
                                                                                      _geom_kind);

    int  n_block     = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *section_ids = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

    for (int i_block = 0; i_block < n_block; i_block++) {
      CWP_Block_t block_type = (CWP_Block_t) PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                                        section_ids[i_block]);

      int block_id = blockAdd(block_type);

      for (int i_part = 0; i_part < _n_part; i_part++) {
        int n_elts = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                 section_ids[i_block],
                                                                 i_part);

        if (block_type == CWP_BLOCK_FACE_POLY) {
          int         *connec     = nullptr;
          int         *connec_idx = nullptr;
          CWP_g_num_t *global_num = nullptr;
          PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                       section_ids[i_block],
                                                       i_part,
                                                       &connec_idx,
                                                       &connec,
                                                       PDM_OWNERSHIP_KEEP);

          global_num = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                           section_ids[i_block],
                                                           i_part,
                                                           PDM_OWNERSHIP_KEEP);

          poly2DBlockSet(i_part,
                         block_id,
                         n_elts,
                         connec_idx,
                         connec,
                         global_num);
        }
        else if (block_type == CWP_BLOCK_CELL_POLY) {
          int          n_face          = 0;
          CWP_g_num_t *face_global_num = nullptr;
          int         *face_vtx_idx    = nullptr;
          int         *face_vtx        = nullptr;
          CWP_g_num_t *cell_global_num = nullptr;
          int         *cell_face_idx   = nullptr;
          int         *cell_face       = nullptr;
          int         *fake_num        = nullptr;
          CWP_g_num_t *fake_global_num = nullptr;
          PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
                                                       section_ids[i_block],
                                                       i_part,
                                                       &n_face,
                                                       &face_global_num,
                                                       &face_vtx_idx,
                                                       &face_vtx,
                                                       &cell_global_num,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       &fake_num,
                                                       &fake_global_num,
                                                       PDM_OWNERSHIP_KEEP);

          poly3DBlockSet(i_part,
                         block_id,
                         n_elts,
                         n_face,
                         face_vtx_idx,
                         face_vtx,
                         cell_face_idx,
                         cell_face,
                         cell_global_num);
        }
        else if (_block_DB[i_block]->blockTypeStdIs()) {
          int         *connec          = nullptr;
          CWP_g_num_t *global_num      = nullptr;
          int         *fake_num        = nullptr;
          CWP_g_num_t *fake_global_num = nullptr;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                    section_ids[i_block],
                                                    i_part,
                                                    &connec,
                                                    &global_num,
                                                    &fake_num,
                                                    &fake_global_num,
                                                    PDM_OWNERSHIP_KEEP);

          stdBlockSet(i_part,
                      block_id,
                      n_elts,
                      connec,
                      global_num);
        }

        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                   section_ids[i_block],
                                                                   i_part,
                                                                   PDM_OWNERSHIP_KEEP);

        _block_DB[block_id]->parentNumSet(i_part,
                                          parent_num,
                                          PDM_OWNERSHIP_USER);
      }

      int i_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(_pmn,
                                                                        _geom_kind,
                                                                        section_ids[i_block]);
      _block_DB[block_id]->blockIdPdmSet(i_section);
    }
  }


  void
  Mesh::geomFinalize()
  {
    int mesh_dimension = -1;
    switch (_cpl->entitiesDimGet()) {
      case CWP_INTERFACE_POINT:
        mesh_dimension = 0;
        break;
      case CWP_INTERFACE_LINEAR:
        mesh_dimension = 1;
        break;
      case CWP_INTERFACE_SURFACE:
        mesh_dimension = 2;
        break;
      case CWP_INTERFACE_VOLUME:
        mesh_dimension = 3;
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "Invalid entity dimension %d\n", _cpl->entitiesDimGet());
    }

    _pmn = PDM_part_mesh_nodal_create(mesh_dimension,
                                      _n_part,
                                      _local_comm_pdm);

    int compute_gnum = 0;
    for (int i = 0; i < _n_part; i++) {
      if (_global_num_vtx[i] == nullptr && _n_vtx[i] > 0) {
        compute_gnum = 1;
        break;
      }
    }

    int _compute_gnum = compute_gnum;
    PDM_MPI_Allreduce(&_compute_gnum, &compute_gnum, 1, PDM_MPI_INT, PDM_MPI_MAX,
                      _local_comm_pdm);

    if (compute_gnum == 1) {
      _is_vtx_gnum_computed = true;
      _gen_gnum(_local_comm_pdm,
                _n_part,
                _n_vtx,
                _coords,
                _global_num_vtx);
    }

    for (int i_part = 0; i_part < _n_part; i_part++) {
      PDM_part_mesh_nodal_coord_set(_pmn,
                                    i_part,
                                    _n_vtx [i_part],
                                    _coords[i_part],
                                    PDM_OWNERSHIP_USER);

      PDM_part_mesh_nodal_vtx_gnum_set(_pmn,
                                       i_part,
                                       _global_num_vtx[i_part],
                                       PDM_OWNERSHIP_USER);
    }

    if (_is_face_edge_method || _is_face_vtx_method) {
      assert(mesh_dimension == 2);

      compute_gnum = 0;
      for (int i_part = 0; i_part < _n_part; i_part++) {
        if (_global_num_faces[i_part] == nullptr && _n_faces[i_part] > 0) {
          compute_gnum = 1;
          break;
        }
      }

      _compute_gnum = compute_gnum;
      PDM_MPI_Allreduce(&_compute_gnum, &compute_gnum, 1, PDM_MPI_INT, PDM_MPI_MAX,
                        _local_comm_pdm);

      if (compute_gnum) {
        _is_elt_gnum_computed = true;
        _gen_gnum_with_connectivity(_is_face_edge_method,
                                    _local_comm_pdm,
                                    _n_part,
                                    _n_faces,
                                    -1,
                                    _face_vtx_idx,
                                    _face_vtx,
                                    _face_edge_idx,
                                    _face_edge,
                                    2,
                                    _face_edge_idx, // useless
                                    _edge_vtx,
                                    _coords,
                                    _global_num_faces);
      }

      if (_is_face_edge_method) {
        assert(!_is_face_vtx_method);
        for (int i_part = 0; i_part < _n_part; i_part++) {
          PDM_part_mesh_nodal_face2d_faceedge_add(_pmn,
                                                  i_part,
                                                  _n_faces         [i_part],
                                                  _n_edges         [i_part],
                                                  _edge_vtx        [i_part],
                                                  _face_edge_idx   [i_part],
                                                  _face_edge       [i_part],
                                                  _global_num_faces[i_part],
                                                  PDM_OWNERSHIP_KEEP);
        }
      }
      else if (_is_face_vtx_method) {
        assert(!_is_face_edge_method);
        for (int i_part = 0; i_part < _n_part; i_part++) {
          PDM_part_mesh_nodal_faces_facevtx_add(_pmn,
                                                i_part,
                                                _n_faces         [i_part],
                                                _face_vtx_idx    [i_part],
                                                _face_vtx        [i_part],
                                                _global_num_faces[i_part],
                                                PDM_OWNERSHIP_KEEP);
        }
      }
      blockDBUpdate();
    }

    else if (_is_cell_face_method || _is_cell_vtx_method) {
      assert(mesh_dimension == 3);

      if (_is_cell_face_method) {
        compute_gnum = 0;
        for (int i_part = 0; i_part < _n_part; i_part++) {
          if (_global_num_faces[i_part] == nullptr && _n_faces[i_part] > 0) {
            compute_gnum = 1;
            break;
          }
        }

        _compute_gnum = compute_gnum;
        PDM_MPI_Allreduce(&_compute_gnum, &compute_gnum, 1, PDM_MPI_INT, PDM_MPI_MAX,
                          _local_comm_pdm);

        if (compute_gnum) {
          _is_face_gnum_computed = true;
          _gen_gnum_with_connectivity(false,
                                      _local_comm_pdm,
                                      _n_part,
                                      _n_faces,
                                      -1,
                                      _face_vtx_idx,
                                      _face_vtx,
                                      _face_edge_idx, // useless
                                      _face_edge,     // useless
                                      -1,
                                      _face_edge_idx, // useless
                                      _edge_vtx,      // useless
                                      _coords,
                                      _global_num_faces);
        }
      }

      compute_gnum = 0;
      for (int i_part = 0; i_part < _n_part; i_part++) {
        if (_global_num_cells[i_part] == nullptr && _n_cells[i_part] > 0) {
          compute_gnum = 1;
          break;
        }
      }

      _compute_gnum = compute_gnum;
      PDM_MPI_Allreduce(&_compute_gnum, &compute_gnum, 1, PDM_MPI_INT, PDM_MPI_MAX,
                        _local_comm_pdm);

      if (compute_gnum) {
        _is_elt_gnum_computed = true;
         _gen_gnum_with_connectivity(_is_cell_face_method,
                                     _local_comm_pdm,
                                     _n_part,
                                     _n_cells,
                                     -1,
                                     _cell_vtx_idx,
                                     _cell_vtx,
                                     _cell_face_idx,
                                     _cell_face,
                                     -1,
                                     _face_vtx_idx,
                                     _face_vtx,
                                     _coords,
                                     _global_num_cells);
      }

      if (_is_cell_face_method) {
        assert(!_is_cell_vtx_method);
        for (int i_part = 0; i_part < _n_part; i_part++) {
          PDM_part_mesh_nodal_cell3d_cellface_add(_pmn,
                                                  i_part,
                                                  _n_cells         [i_part],
                                                  _n_faces         [i_part],
                                                  _face_vtx_idx    [i_part],
                                                  _face_vtx        [i_part],
                                                  _global_num_faces[i_part],
                                                  _cell_face_idx   [i_part],
                                                  _cell_face       [i_part],
                                                  _global_num_cells[i_part],
                                                  PDM_OWNERSHIP_KEEP);
        }
      }
      else if (_is_cell_vtx_method) {
        assert(!_is_cell_face_method);
        for (int i_part = 0; i_part < _n_part; i_part++) {
          PDM_part_mesh_nodal_cells_cellvtx_add(_pmn,
                                                i_part,
                                                _n_cells         [i_part],
                                                _cell_vtx_idx    [i_part],
                                                _cell_vtx        [i_part],
                                                _global_num_cells[i_part],
                                                PDM_OWNERSHIP_KEEP);
        }
      }
      blockDBUpdate();
    }

    else {

      assert(!_is_face_edge_method);
      assert(!_is_face_vtx_method );
      assert(!_is_cell_face_method);
      assert(!_is_cell_vtx_method );

      for (int i_part = 0; i_part < _n_part; i_part++) {

        int idx = 0;
        for (int i_block = 0; i_block < _n_blocks; i_block++) {
          int  n_elts     = _block_DB[i_block]->nEltsGet(i_part);
          int *parent_num = (int *) malloc(sizeof(int) * n_elts);
          for (int i = 0; i < n_elts; i++) {
            parent_num[i] = idx++;
          }

          _block_DB[i_block]->parentNumSet(i_part,
                                           parent_num,
                                           PDM_OWNERSHIP_KEEP);
        }

      }

      for (int i_block = 0; i_block < _n_blocks; i_block++) {

        CWP_Block_t block_type = _block_DB[i_block]->blockTypeGet();

        PDM_Mesh_nodal_elt_t pdm_block_type = (PDM_Mesh_nodal_elt_t) block_type;
        int i_section = PDM_part_mesh_nodal_section_add(_pmn,
                                                        pdm_block_type);
        _block_DB[i_block]->blockIdPdmSet(i_section);

        if (block_type == CWP_BLOCK_CELL_POLY) {
          BlockCP *block = dynamic_cast<BlockCP *>(_block_DB[i_block]);

          std::vector<CWP_g_num_t *> global_num;
          global_num.resize(_n_part, nullptr);

           _gen_gnum_with_connectivity(false,
                                       _local_comm_pdm,
                                       _n_part,
                                       block->nFacesGet(),
                                       -1,
                                       block->faceConnecIdxGet(),
                                       block->faceConnecGet(),
                                       _face_edge_idx, // useless
                                       _face_edge,     // useless
                                       -1,
                                       _face_edge_idx, // useless
                                       _edge_vtx,      // useless
                                       _coords,
                                       global_num);

          block->faceGlobalNumSet(global_num);
        }

      }

      compute_gnum = 0;

      for (int i_part = 0; i_part < _n_part; i_part++) {
        for (int i_block = 0; i_block < _n_blocks; i_block++) {
          if (_block_DB[i_block]->globalNumGet(i_part) == nullptr && _block_DB[i_block]->nEltsGet(i_part) > 0) {
            compute_gnum = 1;
            break;
          }
        }
      }

      _compute_gnum = compute_gnum;
      PDM_MPI_Allreduce(&_compute_gnum, &compute_gnum, 1, PDM_MPI_INT, PDM_MPI_MAX,
                        _local_comm_pdm);

      if (compute_gnum) {
        _is_elt_gnum_computed = true;

        std::vector<double*> cell_centers;
        cell_centers.resize(_n_part, nullptr);

        std::vector<CWP_g_num_t*> global_num;
        global_num.resize(_n_part, nullptr);

        for (int i_part = 0; i_part < _n_part; i_part++) {
          cell_centers[i_part] = PDM_array_const_double(_n_elts[i_part] * 3, 0.);

          for (int i_block = 0; i_block < _n_blocks; i_block++) {

            CWP_Block_t block_type = _block_DB[i_block]->blockTypeGet();

            if (block_type == CWP_BLOCK_FACE_POLY) {
              BlockFP *block = dynamic_cast<BlockFP *>(_block_DB[i_block]);
              _compute_centers_with_connectivity(false,
                                                 true,
                                                 block->nEltsGet    (i_part),
                                                 -1,
                                                 block->connecIdxGet(i_part),
                                                 block->connecGet   (i_part),
                                                 nullptr,
                                                 nullptr,
                                                 -1,
                                                 nullptr,
                                                 nullptr,
                                                 _coords     [i_part],
                                                 block->parentNumGet(i_part),
                                                 cell_centers[i_part]);
            }
            else if (block_type == CWP_BLOCK_CELL_POLY) {
              BlockCP *block = dynamic_cast<BlockCP *>(_block_DB[i_block]);
              _compute_centers_with_connectivity(true,
                                                 true,
                                                 block->nEltsGet        (i_part),
                                                 -1,
                                                 nullptr,
                                                 nullptr,
                                                 block->cellConnecIdxGet(i_part),
                                                 block->cellConnecGet   (i_part),
                                                 -1,
                                                 block->faceConnecIdxGet(i_part),
                                                 block->faceConnecGet   (i_part),
                                                 _coords     [i_part],
                                                 block->parentNumGet    (i_part),
                                                 cell_centers[i_part]);
            }
            else if (_block_DB[i_block]->blockTypeStdIs()) {
              BlockStd *block = dynamic_cast<BlockStd *>(_block_DB[i_block]);
              PDM_Mesh_nodal_elt_t elt_type = (PDM_Mesh_nodal_elt_t) block_type;
              int cell_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);
              _compute_centers_with_connectivity(false,
                                                 true,
                                                 block->nEltsGet    (i_part),
                                                 cell_n_vtx,
                                                 nullptr,
                                                 block->connecGet   (i_part),
                                                 nullptr,
                                                 nullptr,
                                                 -1,
                                                 nullptr,
                                                 nullptr,
                                                 _coords     [i_part],
                                                 block->parentNumGet(i_part),
                                                 cell_centers[i_part]);
            }
            else {
              BlockHO *block = dynamic_cast<BlockHO *>(_block_DB[i_block]);
              PDM_Mesh_nodal_elt_t elt_type = (PDM_Mesh_nodal_elt_t) block_type;
              int cell_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, block->orderGet());
              _compute_centers_with_connectivity(false,
                                                 true,
                                                 block->nEltsGet    (i_part),
                                                 cell_n_vtx,
                                                 nullptr,
                                                 block->connecGet   (i_part),
                                                 nullptr,
                                                 nullptr,
                                                 -1,
                                                 nullptr,
                                                 nullptr,
                                                 _coords     [i_part],
                                                 block->parentNumGet(i_part),
                                                 cell_centers[i_part]);
            }
          }
        }

        _gen_gnum(_local_comm_pdm,
                  _n_part,
                  _n_elts,
                  cell_centers,
                  global_num);

        for (int i_part = 0; i_part < _n_part; i_part++) {

          for (int i_block = 0; i_block < _n_blocks; i_block++) {
            int          n_elts           = _block_DB[i_block]->nEltsGet(i_part);
            CWP_g_num_t *global_num_block = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
            if (n_elts > 0) {
              int idx_read = _block_DB[i_block]->parentNumGet(i_part)[0];
              memcpy(global_num_block, &global_num[i_part][idx_read], n_elts * sizeof(CWP_g_num_t));
            }

            _block_DB[i_block]->globalNumSet(i_part,
                                             global_num_block,
                                             PDM_OWNERSHIP_KEEP);
          }

          free(cell_centers[i_part]);
          free(global_num  [i_part]);

        }
      }

      for (int i_part = 0; i_part < _n_part; i_part++) {

        for (int i_block = 0; i_block < _n_blocks; i_block++) {

          CWP_Block_t block_type   = _block_DB[i_block]->blockTypeGet();
          int         pdm_id_block = _block_DB[i_block]->blockIdPdmGet();
          int         n_elts       = _block_DB[i_block]->nEltsGet(i_part);

          if (block_type == CWP_BLOCK_FACE_POLY) {
            BlockFP *block = dynamic_cast<BlockFP *>(_block_DB[i_block]);
            PDM_part_mesh_nodal_section_poly2d_set(_pmn,
                                                   pdm_id_block,
                                                   i_part,
                                                   n_elts,
                                                   block->connecIdxGet(i_part),
                                                   block->connecGet   (i_part),
                                                   block->globalNumGet(i_part),
                                                   block->parentNumGet(i_part),
                                                   PDM_OWNERSHIP_USER);
          }
          else if (block_type == CWP_BLOCK_CELL_POLY) {
            BlockCP *block = dynamic_cast<BlockCP *>(_block_DB[i_block]);
            PDM_part_mesh_nodal_section_poly3d_set(_pmn,
                                                   pdm_id_block,
                                                   i_part,
                                                   n_elts,
                                                   block->nFacesGet       (i_part),
                                                   block->faceConnecIdxGet(i_part),
                                                   block->faceConnecGet   (i_part),
                                                   block->faceGlobalNumGet(i_part),
                                                   block->cellConnecIdxGet(i_part),
                                                   block->cellConnecGet   (i_part),
                                                   block->globalNumGet    (i_part),
                                                   block->parentNumGet    (i_part),
                                                   nullptr,
                                                   PDM_OWNERSHIP_USER);
          }
          else if (_block_DB[i_block]->blockTypeStdIs()) {
            BlockStd *block = dynamic_cast<BlockStd *>(_block_DB[i_block]);
            PDM_part_mesh_nodal_section_std_set(_pmn,
                                                pdm_id_block,
                                                i_part,
                                                n_elts,
                                                block->connecGet   (i_part),
                                                block->globalNumGet(i_part),
                                                block->parentNumGet(i_part),
                                                nullptr,
                                                PDM_OWNERSHIP_USER);
          }
          else {
            BlockHO *block = dynamic_cast<BlockHO *>(_block_DB[i_block]);
            PDM_part_mesh_nodal_section_std_ho_set(_pmn,
                                                   pdm_id_block,
                                                   i_part,
                                                   n_elts,
                                                   block->connecIJKGet(i_part),
                                                   block->globalNumGet(i_part),
                                                   block->parentNumGet(i_part),
                                                   nullptr,
                                                   block->orderGet(),
                                                   nullptr,
                                                   PDM_OWNERSHIP_USER);
          }
        }
      }
    }

    for (int i_block = 0 ; i_block < _n_blocks ; i_block++) {
      _block_DB[i_block]->geomFinalize();
    }

    _n_blocks = PDM_part_mesh_nodal_n_section_in_geom_kind_get(_pmn,
                                                              _geom_kind);

  }


  void
  Mesh::stdBlockSet
  (
    const int         i_part,
    const int         block_id,
    const int         n_elts,
          int         connec[],
          CWP_g_num_t global_num[]
  )
  {
    BlockStd *block = dynamic_cast<BlockStd *>(_block_DB[block_id]);
    block->blockSet(i_part,
                    n_elts,
                    connec,
                    global_num);
    _n_elts[i_part] += n_elts;
  }


  void
  Mesh::stdBlockGet
  (
    const int           i_part,
    const int           block_id,
          int          *n_elts,
          int         **connec,
          CWP_g_num_t **global_num
  )
  {
    BlockStd *block = dynamic_cast<BlockStd *>(_block_DB[block_id]);
    block->blockGet(i_part,
                    n_elts,
                    connec,
                    global_num);
  }


  void
  Mesh::hoBlockSet
  (
    const int          i_part,
    const int          block_id,
    const int          n_elts,
          int          connec[],
          CWP_g_num_t  global_num[],
    const int          order,
    const char        *ho_ordering
  )
  {
    BlockHO *block = dynamic_cast<BlockHO *>(_block_DB[block_id]);
    block->blockSet(i_part,
                    n_elts,
                    connec,
                    global_num,
                    order,
                    ho_ordering);
    _n_elts[i_part] += n_elts;
  }


  void
  Mesh::hoBlockGet
  (
    const int           i_part,
    const int           block_id,
          int          *n_elts,
          int          *order,
          int         **connec,
          CWP_g_num_t **global_num
  )
  {
    BlockHO *block = dynamic_cast<BlockHO *>(_block_DB[block_id]);
    char *ho_ordering;
    block->blockGet(i_part,
                    n_elts,
                    connec,
                    global_num,
                    order,
                    &ho_ordering);
  }


  void
  Mesh::poly2DBlockSet
  (
    const int         i_part,
    const int         block_id,
    const int         n_elts,
          int         connec_idx[],
          int         connec[],
          CWP_g_num_t global_num[]
  )
  {
    BlockFP *block = dynamic_cast<BlockFP *>(_block_DB[block_id]);
    block->blockSet(i_part,
                    n_elts,
                    connec_idx,
                    connec,
                    global_num);
    _n_elts[i_part] += n_elts;
  }


  void
  Mesh::poly2DBlockGet
  (
    const int           i_part,
    const int           block_id,
          int          *n_elts,
          int         **connec_idx,
          int         **connec,
          CWP_g_num_t **global_num
  )
  {
    BlockFP *block = dynamic_cast<BlockFP *>(_block_DB[block_id]);
    block->blockGet(i_part,
                    n_elts,
                    connec_idx,
                    connec,
                    global_num);
  }


  void
  Mesh::poly3DBlockSet
  (
    const int         i_part,
    const int         block_id,
    const int         n_elts,
    const int         n_faces,
          int         connec_faces_idx[],
          int         connec_faces[],
          int         connec_cells_idx[],
          int         connec_cells[],
          CWP_g_num_t global_num[]
  )
  {
    BlockCP *block = dynamic_cast<BlockCP *>(_block_DB[block_id]);
    block->blockSet(i_part,
                    n_elts,
                    n_faces,
                    connec_faces_idx,
                    connec_faces,
                    connec_cells_idx,
                    connec_cells,
                    global_num);
    _n_elts[i_part] += n_elts;
  }


  void
  Mesh::poly3DBlockGet
  (
    const int           i_part,
    const int           block_id,
          int          *n_elts,
          int          *n_faces,
          int         **connec_faces_idx,
          int         **connec_faces,
          int         **connec_cells_idx,
          int         **connec_cells,
          CWP_g_num_t **global_num
  )
  {
    BlockCP *block = dynamic_cast<BlockCP *>(_block_DB[block_id]);
    block->blockGet(i_part,
                    n_elts,
                    n_faces,
                    connec_faces_idx,
                    connec_faces,
                    connec_cells_idx,
                    connec_cells,
                    global_num);
  }


  void
  Mesh::meshDel()
  {
    if (_pmn != nullptr) {
      PDM_part_mesh_nodal_free(_pmn);
    }
    _pmn = nullptr;

    for (int i = 0; i < _n_part; i++) {
      if (_elt_centers[i] != nullptr) {
        free(_elt_centers[i]);
        _elt_centers[i] = nullptr;
      }

      for (int i_block = 0; i_block < _n_blocks; i_block++) {
        _block_DB[i_block]->parentNumFree(i);
        _block_DB[i_block]->globalNumFree(i);

        CWP_Block_t block_type = _block_DB[i_block]->blockTypeGet();
        if (block_type == CWP_BLOCK_CELL_POLY) {
          BlockCP *block = dynamic_cast<BlockCP *>(_block_DB[i_block]);
          block->faceGlobalNumFree(i);
        }
      }

      if (_is_vtx_gnum_computed) {
        if (_global_num_vtx[i] != nullptr) {
          free(_global_num_vtx[i]);
          _global_num_vtx[i] = nullptr;
        }
      }

      if (_is_elt_gnum_computed) {
        if (_is_face_edge_method || _is_face_vtx_method) {
          if (_global_num_faces[i] != nullptr) {
            free(_global_num_faces[i]);
            _global_num_faces[i] = nullptr;
          }
        }
        else if (_is_cell_face_method || _is_cell_vtx_method) {
          if (_global_num_cells[i] != nullptr) {
            free(_global_num_cells[i]);
            _global_num_cells[i] = nullptr;
          }
        }
      }

      if (_is_face_gnum_computed) {
        if (_global_num_faces[i] != nullptr) {
          free(_global_num_faces[i]);
          _global_num_faces[i] = nullptr;
        }
      }
    }

    _n_blocks = 0;
    for (int i = 0; i < (int) _block_DB.size(); i++) {
      delete _block_DB[i];
    }
    _block_DB.clear();
  }


  int
  Mesh::blockAdd
  (
    const CWP_Block_t block_type
  )
  {
    Block *block = FB::getInstance().CreateObject(block_type);
    block->blockAdd(block_type, this);
    int block_id_cwipi = _n_blocks;
    _block_DB.push_back(block);
    block->blockIdCwipiSet(block_id_cwipi);
    _n_blocks = _block_DB.size();
    return block->blockIdCwipiGet();
  }


  void
  Mesh::fromCellFaceSet
  (
    const int         i_part,
    const int         n_cells,
          int         cell_face_idx[],
          int         cell_face[],
          int         n_faces,
          int         face_vtx_idx[],
          int         face_vtx[],
          CWP_g_num_t global_num[]
  )
  {
    _is_cell_face_method = true;
    _n_cells         [i_part] = n_cells;
    _cell_face_idx   [i_part] = cell_face_idx;
    _cell_face       [i_part] = cell_face;
    _global_num_cells[i_part] = global_num;
    _n_faces         [i_part] = n_faces;
    _face_vtx_idx    [i_part] = face_vtx_idx;
    _face_vtx        [i_part] = face_vtx;
  }


  void
  Mesh::fromFaceEdgeSet
  (
    const int         i_part,
    const int         n_faces,
          int         face_edge_idx[],
          int         face_edge[],
    const int         n_edges,
          int         edge_vtx[],
          CWP_g_num_t global_num[]
  )
  {
    _is_face_edge_method = true;
    _n_faces         [i_part] = n_faces;
    _face_edge_idx   [i_part] = face_edge_idx;
    _face_edge       [i_part] = face_edge;
    _global_num_faces[i_part] = global_num;
    _n_edges         [i_part] = n_edges;
    _edge_vtx        [i_part] = edge_vtx;
  }


  void
  Mesh::fromFaceVtxSet
  (
    const int         i_part,
    const int         n_faces,
          int         face_vtx_idx[],
          int         face_vtx[],
          CWP_g_num_t global_num[]
  )
  {
    _is_face_vtx_method = true;
    _n_faces         [i_part] = n_faces;
    _face_vtx_idx    [i_part] = face_vtx_idx;
    _face_vtx        [i_part] = face_vtx;
    _global_num_faces[i_part] = global_num;
  }


  void
  Mesh::fromCellVtxSet
  (
    const int         i_part,
    const int         n_cells,
          int         cell_vtx_idx[],
          int         cell_vtx[],
          CWP_g_num_t global_num[]
  )
  {
    _is_cell_vtx_method = true;
    _n_cells         [i_part] = n_cells;
    _cell_vtx_idx    [i_part] = cell_vtx_idx;
    _cell_vtx        [i_part] = cell_vtx;
    _global_num_cells[i_part] = global_num;
  }


  int
  Mesh::nEltsGet
  (
    int id_part
  ) const
  {
    return PDM_part_mesh_nodal_n_elmts_get(_pmn, _geom_kind, id_part);
  }

}

/**
 * \endcond
 */