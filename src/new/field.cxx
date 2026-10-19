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

#include "field.hxx"

#include "coupling.hxx"
#include "coupling_i.hxx"
#include "field_i.hxx"
#include "mesh.hxx"

#include "pdm_writer.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <string>

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  Field::Field
  (
    string              field_id,
    int                 fieldIDInt,
    CWP_Type_t          dataType,
    Coupling*           cpl,
    CWP_Dof_location_t  fieldType,
    CWP_Field_storage_t storage,
    int                 nComponent,
    CWP_Field_exch_t    exchangeType,
    CWP_Status_t        visuStatus
  )
  : _field_id(field_id),
    _field_id_int(fieldIDInt),
    _data_type(dataType),
    _field_location(fieldType),
    _linked_field_location(CWP_DOF_LOCATION_UNDEF),
    _n_component(nComponent),
    _storage(storage),
    _exchange_type(exchangeType),
    _visu_status(visuStatus),
    _cpl(cpl),
    _mesh(cpl->meshGet()),
    _is_send_yet(false),
    _is_recv_yet(false),
    _writer(nullptr),
    _id_writer_var_send(nullptr),
    _id_writer_var_recv(nullptr),
    _id_writer_var_send_status(-1),
    _id_writer_var_recv_status(-1),
    _interpolation_type(CWP_INTERPOLATION_DEFAULT),
    _interpolation_function(nullptr),
    _interpolation_function_f(nullptr),
    _interpolation_function_p(nullptr),
    _python_object(nullptr),
    _is_computed_tgt_bcast_enabled(false),
    _is_involved_src_bcast_enabled(false)
  {
    _n_part   = _mesh->nPartGet();
    _data_tgt = static_cast<void **>(malloc(sizeof(void *) * _n_part));
    _data_src = static_cast<void **>(malloc(sizeof(void *) * _n_part));

    switch (_data_type) {
      case CWP_DOUBLE:
        _data_type_size = sizeof(double);
        break;
      case CWP_INT:
        _data_type_size = sizeof(int);
        break;
      default:
        PDM_error("CWP_CHAR not usable");
    }

    if (visuStatus == CWP_STATUS_ON) {
      _writer = cpl->writerGet();
    }

    if (_writer != nullptr) {

      _id_writer_var_send = static_cast<int *>(malloc(sizeof(int) * _n_component));
      _id_writer_var_recv = static_cast<int *>(malloc(sizeof(int) * _n_component));
      for (int i = 0; i < _n_component; i++) {
        _id_writer_var_send[i] = -1;
        _id_writer_var_recv[i] = -1;
      }

      // Set field location type
      PDM_writer_var_loc_t pdm_field_type = PDM_WRITER_VAR_ELEMENTS;
      if      (_field_location == CWP_DOF_LOCATION_CELL_CENTER) {
        pdm_field_type = PDM_WRITER_VAR_ELEMENTS;
      }
      else if (_field_location == CWP_DOF_LOCATION_NODE) {
        pdm_field_type = PDM_WRITER_VAR_VERTICES;
      }
      else if (_field_location == CWP_DOF_LOCATION_USER) {
        // Vertices by default. PDM_WRITER_VAR_PARTICLES not implemented
        pdm_field_type = PDM_WRITER_VAR_VERTICES;
      }

      // Create send variables
      if (_exchange_type == CWP_FIELD_EXCH_SEND ||
          _exchange_type == CWP_FIELD_EXCH_SENDRECV) {
        string prefix = "s";

        for (int i_comp = 0; i_comp < _n_component; i_comp++) {

          ostringstream num;
          num << (i_comp + 1);
          string fieldName = prefix + "_" + _field_id + num.str();

          _id_writer_var_send[i_comp] = PDM_writer_var_create(_writer,
                                                              PDM_WRITER_ON,
                                                              PDM_WRITER_VAR_SCALAR,
                                                              pdm_field_type,
                                                              fieldName.c_str());
        }

        // Create variable to tag if field has been exchanged (always time-dependent)
        string fieldName = prefix + "_" + _field_id;
        string fieldStatusName = fieldName + "_status";

        _id_writer_var_send_status = PDM_writer_var_create(_writer,
                                                           PDM_WRITER_ON,
                                                           PDM_WRITER_VAR_SCALAR,
                                                           pdm_field_type,
                                                           fieldStatusName.c_str());
      }

      // Create receive variables
      if (_exchange_type == CWP_FIELD_EXCH_RECV ||
          _exchange_type == CWP_FIELD_EXCH_SENDRECV) {
        string prefix = "r";

        for (int i_comp = 0; i_comp < _n_component; i_comp++) {

          ostringstream num;
          num << (i_comp + 1);
          string fieldName = prefix + "_" + _field_id + num.str();

          _id_writer_var_recv[i_comp] = PDM_writer_var_create(_writer,
                                                              PDM_WRITER_ON,
                                                              PDM_WRITER_VAR_SCALAR,
                                                              pdm_field_type,
                                                              fieldName.c_str());
        }

        // Create variable to tag if field has been computed/exchanged
        string fieldName = prefix + "_" + _field_id;
        string fieldStatusName = fieldName + "_status";

        _id_writer_var_recv_status = PDM_writer_var_create(_writer,
                                                           PDM_WRITER_ON,
                                                           PDM_WRITER_VAR_SCALAR,
                                                           pdm_field_type,
                                                           fieldStatusName.c_str());
      }
    }
  }


  Field::~Field()
  {
    if (_data_tgt != nullptr) {
      free(_data_tgt);
    }
    if (_data_src != nullptr) {
      free(_data_src);
    }

    if (_id_writer_var_send != nullptr) {
      free(_id_writer_var_send);
    }
    if (_id_writer_var_recv != nullptr) {
      free(_id_writer_var_recv);
    }
  }


  int
  Field::nDOFGet
  (
    const int i_part
  )
  {
    if      (_field_location == CWP_DOF_LOCATION_CELL_CENTER) {
      return _mesh->nEltsGet(i_part);
    }
    else if (_field_location == CWP_DOF_LOCATION_NODE) {
      return _mesh->nVtxGet(i_part);
    }
    else if (_field_location == CWP_DOF_LOCATION_USER) {
      return _cpl->userTargetNGet(i_part);
    }
    return 0;
  }


  void
  Field::write
  (
    const CWP_Field_exch_t exch_type
  )
  {
    assert(CWP_FIELD_EXCH_SENDRECV != exch_type);

    if (_writer == nullptr) {
      return;
    }

    if (_cpl->NStepGet() % _cpl->freqWriterGet() != 0) {
      return;
    }

    if (visuStatusGet() == CWP_STATUS_OFF) {
      return;
    }

    // Get data depending on exchange type
    int   *id_writer_var;
    int    id_writer_var_status;
    void **data_var;
    bool   is_exch;

    if (exch_type == CWP_FIELD_EXCH_SEND) {
      id_writer_var        = _id_writer_var_send;
      id_writer_var_status = _id_writer_var_send_status;
      data_var             = _data_src;
      is_exch              = _is_send_yet;
    }
    else {
      id_writer_var        = _id_writer_var_recv;
      id_writer_var_status = _id_writer_var_recv_status;
      data_var             = _data_tgt;
      is_exch              = _is_recv_yet;
    }

    // Get maximum number of DOFs
    int n_dof_max = 0;
    for (int i_part = 0; i_part < _n_part; i_part++) {
      n_dof_max = max(n_dof_max, nDOFGet(i_part));
    }

    double *_comp_data = static_cast<double *>(malloc(sizeof(double) * n_dof_max));
    double *comp_data  = nullptr;

    // Fill buffers and write to file
    for (int i_part = 0; i_part < _n_part; i_part++) {

      int     n_dof     = nDOFGet(i_part);
      double *_data_var = static_cast<double *>(data_var[i_part]);

      // Get data per component
      for (int i_comp = 0; i_comp < _n_component; i_comp++) {

        // Default initialization if field not exchanged or for unlocalized DOFs
        if (!is_exch || exch_type == CWP_FIELD_EXCH_RECV) {
          for (int i = 0; i < n_dof; i++) {
            _comp_data[i] = 123456789;
          }
          comp_data = _comp_data;
        }

        if (is_exch) {
          if (exch_type == CWP_FIELD_EXCH_SEND) {
            // Deinterlace data if necessary
            if (_storage == CWP_FIELD_STORAGE_INTERLACED && _n_component > 1) {
              for (int i = 0; i < n_dof; i++) {
                _comp_data[i] = _data_var[_n_component * i + i_comp];
              }
              comp_data = _comp_data;
            }
            // Pointer to component to write
            else {
              comp_data = _data_var + i_comp * n_dof;
            }
          }
          else {
            // Copy data for localized DOFs
            const int *computed_target   = nullptr;
            int        n_computed_target = 0;
            if (_cpl->has_mesh()) {
              computed_target   = _cpl->computedTargetsGet (_field_id, i_part);
              n_computed_target = _cpl->nComputedTargetsGet(_field_id, i_part);
            }

            if (_storage == CWP_FIELD_STORAGE_INTERLACED) {
              for (int i = 0; i < n_computed_target; i++) {
                int j = computed_target[i] - 1;
                _comp_data[j] = _data_var[_n_component * i + i_comp];
              }
            }
            else {
              for (int i = 0; i < n_computed_target; i++) {
                int j = computed_target[i] - 1;
                _comp_data[j] = _data_var[i_comp * n_computed_target + i];
              }
            }
          }
        }

        PDM_writer_var_set(_writer,
                           id_writer_var[i_comp],
                           _cpl->idGeomWriterGet(_field_location),
                           i_part,
                           comp_data);
      }

      // Get status
      if (!is_exch) {
        for (int i = 0; i < n_dof; i++) {
          _comp_data[i] = -1;
        }
      }
      else {
        if (exch_type == CWP_FIELD_EXCH_SEND) {
          for (int i = 0; i < n_dof; i++) {
            _comp_data[i] = 0;
          }
        }
        else {
          const int *computed_target     = nullptr;
          int        n_computed_target   = 0;
          const int *uncomputed_target   = nullptr;
          int        n_uncomputed_target = 0;
          if (_cpl->has_mesh()) {
            computed_target     = _cpl->computedTargetsGet   (_field_id, i_part);
            n_computed_target   = _cpl->nComputedTargetsGet  (_field_id, i_part);
            uncomputed_target   = _cpl->uncomputedTargetsGet (_field_id, i_part);
            n_uncomputed_target = _cpl->nUncomputedTargetsGet(_field_id, i_part);
          }

          for (int i = 0; i < n_uncomputed_target; i++) {
            _comp_data[uncomputed_target[i] - 1] = 1;
          }
          for (int i = 0; i < n_computed_target; i++) {
            _comp_data[computed_target[i] - 1] = 0;
          }
        }
      }

      PDM_writer_var_set(_writer,
                         id_writer_var_status,
                         _cpl->idGeomWriterGet(_field_location),
                         i_part,
                         _comp_data);
    }

    for (int i_comp = 0; i_comp < _n_component; i_comp++) {
      PDM_writer_var_write    (_writer, id_writer_var[i_comp]);
      PDM_writer_var_data_free(_writer, id_writer_var[i_comp]);
    }
    PDM_writer_var_write    (_writer, id_writer_var_status);
    PDM_writer_var_data_free(_writer, id_writer_var_status);

    free(_comp_data);
  }

}

/**
 * \endcond
 */