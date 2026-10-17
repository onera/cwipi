#ifndef __FIELD_I_H__
#define __FIELD_I_H__
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

namespace cwipi {

  inline std::string
  Field::fieldIDGet() const
  {
    return _field_id;
  }


  inline int
  Field::fieldIDIntGet() const
  {
    return _field_id_int;
  }


  inline CWP_Type_t
  Field::dataTypeGet() const
  {
    return _data_type;
  }


  inline int
  Field::dataTypeSizeGet() const
  {
    return _data_type_size;
  }


  inline CWP_Dof_location_t
  Field::locationGet() const
  {
    return _field_location;
  }


  inline void
  Field::linkedFieldLocationSet
  (
    const CWP_Dof_location_t linkedFieldLocation
  )
  {
    _linked_field_location = linkedFieldLocation;
  }


  inline CWP_Dof_location_t
  Field::linkedFieldLocationGet() const
  {
    return _linked_field_location;
  }


  inline CWP_Field_storage_t
  Field::storageTypeGet() const
  {
    return _storage;
  }


  inline int
  Field::nComponentGet() const
  {
    return _n_component;
  }


  inline CWP_Field_exch_t
  Field::exchangeTypeGet() const
  {
    return _exchange_type;
  }


  inline CWP_Status_t
  Field::visuStatusGet() const
  {
    return _visu_status;
  }


  inline int
  Field::nPartGet() const
  {
    return _n_part;
  }


  inline void
  Field::dataSet
  (
    const int             i_part,
    const CWP_Field_map_t map_type,
          void*           data
  )
  {
    if (map_type == CWP_FIELD_MAP_SOURCE) {
      _data_src[i_part] = data;
    }
    else if (map_type == CWP_FIELD_MAP_TARGET) {
      _data_tgt[i_part] = data;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error : unknown map type.\n");
    }
  }


  inline void*
  Field::dataGet
  (
    const int             i_part,
    const CWP_Field_map_t map_type
  ) const
  {
    if (map_type == CWP_FIELD_MAP_SOURCE) {
      return _data_src[i_part];
    }
    else if (map_type == CWP_FIELD_MAP_TARGET) {
      return _data_tgt[i_part];
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error : unknown map type.\n");
    }
    return nullptr;
  }


  inline void
  Field::interpFunctionSet
  (
    const CWP_Interp_function_t fct
  )
  {
    _interpolation_type     = CWP_INTERPOLATION_USER;
    _interpolation_function = fct;
  }


  inline CWP_Interp_function_t
  Field::interpFunctionGet() const
  {
    return _interpolation_function;
  }


  inline void
  Field::interpFunctionUnset()
  {
    _interpolation_type     = CWP_INTERPOLATION_DEFAULT;
    _interpolation_function = nullptr;
  }


  inline void
  Field::interpFunctionFSet
  (
    const CWP_Interp_function_t fct
  )
  {
    _interpolation_type       = CWP_INTERPOLATION_USER;
    _interpolation_function_f = fct;
  }


  inline CWP_Interp_function_t
  Field::interpFunctionFGet() const
  {
    return _interpolation_function_f;
  }


  inline void
  Field::interpFunctionFUnset()
  {
    _interpolation_type       = CWP_INTERPOLATION_DEFAULT;
    _interpolation_function_f = nullptr;
  }


  inline void
  Field::interpFunctionPSet
  (
    const CWP_Interp_function_p_t fct
  )
  {
    _interpolation_type       = CWP_INTERPOLATION_USER;
    _interpolation_function_p = fct;
  }


  inline CWP_Interp_function_p_t
  Field::interpFunctionPGet() const
  {
    return _interpolation_function_p;
  }


  inline void
  Field::interpFunctionPUnset()
  {
    _interpolation_type       = CWP_INTERPOLATION_DEFAULT;
    _interpolation_function_p = nullptr;
  }


  inline CWP_Interpolation_t
  Field::interpolationTypeGet() const
  {
    return _interpolation_type;
  }


  inline void
  Field::pythonObjectSet
  (
    void *p
  )
  {
    _python_object = p;
  }


  inline void *
  Field::pythonObjectGet() const
  {
    return _python_object;
  }


  inline void
  Field::isSendYetSet
  (
    const bool value
  )
  {
    _is_send_yet = value;
  }


  inline bool
  Field::isSendYetGet() const
  {
    return _is_send_yet;
  }


  inline void
  Field::isRecvYetSet
  (
    const bool value
  )
  {
    _is_recv_yet = value;
  }


  inline bool
  Field::isRecvYetGet() const
  {
    return _is_recv_yet;
  }


  inline int
  Field::idWriterVarSendGet
  (
    const int i_comp
  ) const
  {
    assert(i_comp < _n_component);
    return _id_writer_var_send[i_comp];
  }


  inline int
  Field::idWriterVarSendStatusGet() const
  {
    return _id_writer_var_send_status;
  }


  inline int
  Field::idWriterVarRecvGet
  (
    const int i_comp
  ) const
  {
    assert(i_comp < _n_component);
    return _id_writer_var_recv[i_comp];
  }


  inline int
  Field::idWriterVarRecvStatusGet() const
  {
    return _id_writer_var_recv_status;
  }


  inline void
  Field::computedTgtBcastEnable()
  {
    _is_computed_tgt_bcast_enabled = true;
  }


  inline bool
  Field::computedTgtBcastEnabledIs() const
  {
    return _is_computed_tgt_bcast_enabled;
  }


  inline void
  Field::involvedSrcBcastEnable()
  {
    _is_involved_src_bcast_enabled = true;
  }


  inline bool
  Field::involvedSrcBcastEnabledIs() const
  {
    return _is_involved_src_bcast_enabled;
  }

}

#endif // __FIELD_I_H__