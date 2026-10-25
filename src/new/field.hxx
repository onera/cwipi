#ifndef __FIELD_H__
#define __FIELD_H__
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

#include "cwp_priv.h"

#include "cwp.h"

#include "pdm_writer.h"

#include <string>

/**
 * \cond
 */

namespace cwipi {

  class Coupling;
  class Mesh;

  /**
   * \class Field field.hxx "field.hxx"
   * \brief Abstract field.
   *
   * This class is the field abstract interface
   *
   */
  class Field {

  public:

    /**
     * \brief Constructor.
     *
     * \param [in]  field_id     Field identifier
     * \param [in]  fieldIDInt   Integer field identifier
     * \param [in]  dataType     Field datatype
     * \param [in]  cpl          Coupling
     * \param [in]  fieldType    Value location
     * \param [in]  storage      Storage type
     * \param [in]  nComponent   Number of field components
     * \param [in]  exchangeType Exchange type
     * \param [in]  visuStatus   Visualization status
     *
     */
    Field
    (
      const std::string         field_id,
      const int                 fieldIDInt,
      const CWP_Type_t          dataType,
            Coupling*           cpl,
      const CWP_Dof_location_t  fieldType,
      const CWP_Field_storage_t storage,
      const int                 nComponent,
      const CWP_Field_exch_t    exchangeType,
      const CWP_Status_t        visuStatus
    );


    /**
     * \brief Destructor.
     *
     */
    ~Field();


    /**
     * \brief Get field identifier.
     *
     * \return Field identifier
     *
     */
    inline std::string
    fieldIDGet() const;


    /**
     * \brief Get integer field identifier.
     *
     * \return Integer field identifier
     *
     */
    inline int
    fieldIDIntGet() const;


    /**
     * \brief Get field data type.
     *
     * \return Field data type
     *
     */
    inline CWP_Type_t
    dataTypeGet() const;


    /**
     * \brief Get field data type size.
     *
     * \return Field data type size
     *
     */
    inline int
    dataTypeSizeGet() const;


    /**
     * \brief Get field location.
     *
     * \return Field location
     *
     */
    inline CWP_Dof_location_t
    locationGet() const;


    /**
     * \brief Set linked field location.
     *
     * \param [in]  linkedFieldLocation Linked field location
     *
     */
    inline void
    linkedFieldLocationSet
    (
      const CWP_Dof_location_t linkedFieldLocation
    );


    /**
     * \brief Get linked field location.
     *
     * \return Linked field location
     *
     */
    inline CWP_Dof_location_t
    linkedFieldLocationGet() const;


    /**
     * \brief Get field storage type.
     *
     * \return Field storage type
     *
     */
    inline CWP_Field_storage_t
    storageTypeGet() const;


    /**
     * \brief Get number of field components.
     *
     * \return Number of field components
     *
     */
    inline int
    nComponentGet() const;


    /**
     * \brief Get exchange type.
     *
     * \return Exchange type
     *
     */
    inline CWP_Field_exch_t
    exchangeTypeGet() const;


    /**
     * \brief Get visualization status.
     *
     * \return Visualization status
     *
     */
    inline CWP_Status_t
    visuStatusGet() const;


    /**
     * \brief Get number of field partitions.
     *
     * \return Number of partitions
     *
     */
    inline int
    nPartGet() const;


    /**
     * \brief Get number of field degrees of freedom.
     *
     * \param [in]  i_part Partition identifier
     *
     * \return Number of field degrees of freedom
     *
     */
    int
    nDOFGet
    (
      const int i_part
    );


    /**
     * \brief Set data array.
     *
     * \param [in]  i_part   Partition identifier
     * \param [in]  map_type Data map type
     * \param [in]  data     Data array
     *
     */
    inline void
    dataSet
    (
      const int             i_part,
      const CWP_Field_map_t map_type,
            void*           data
    );


    /**
     * \brief Get data array.
     *
     * \param [in]  i_part   Partition identifier
     * \param [in]  map_type Data map type
     *
     * \return Data array
     *
     */
    inline void*
    dataGet
    (
      const int             i_part,
      const CWP_Field_map_t map_type
    ) const;


    /**
     * \brief Write field.
     *
     */
    void
    write
    (
      const CWP_Field_exch_t exch_type
    );


    /**
     * \brief Set interpolation function.
     *
     * \param [in]  fct Interpolation function
     *
     */
    inline void
    interpFunctionSet
    (
      const CWP_Interp_function_t fct
    );


    /**
     * \brief Get interpolation function.
     *
     * \return Interpolation function
     *
     */
    inline CWP_Interp_function_t
    interpFunctionGet() const;


    /**
     * \brief Unset interpolation function.
     *
     */
    inline void
    interpFunctionUnset();


    /**
     * \brief Set Fortran interpolation function.
     *
     * \param [in]  fct Interpolation function
     *
     */
    inline void
    interpFunctionFSet
    (
      const CWP_Interp_function_t fct
    );


    /**
     * \brief Get Fortran interpolation function.
     *
     * \return Interpolation function
     *
     */
    inline CWP_Interp_function_t
    interpFunctionFGet() const;


    /**
     * \brief Unset Fortran interpolation function.
     *
     */
    inline void
    interpFunctionFUnset();


    /**
     * \brief Set Python interpolation function.
     *
     * \param [in]  fct Interpolation function
     *
     */
    inline void
    interpFunctionPSet
    (
      const CWP_Interp_function_p_t fct
    );


    /**
     * \brief Get Python interpolation function.
     *
     * \return Interpolation function
     *
     */
    inline CWP_Interp_function_p_t
    interpFunctionPGet() const;


    /**
     * \brief Unset Python interpolation function.
     *
     */
    inline void
    interpFunctionPUnset();


    /**
     * \brief Get interpolation type.
     *
     * \return Interpolation type
     *
     */
    inline CWP_Interpolation_t
    interpolationTypeGet() const;


    /**
     * \brief Set python object.
     *
     * \param [in]  p Python object
     *
     */
    inline void
    pythonObjectSet
    (
      void *p
    );


    /**
     * \brief Get python object.
     *
     * \return Python object
     *
     */
    inline void *
    pythonObjectGet() const;


    /**
     * \brief Set send already performed flag.
     *
     * \param [in]  value Flag value
     *
     */
    inline void
    isSendYetSet
    (
      const bool value
    );


    /**
     * \brief Get send already performed flag.
     *
     * \return Send already performed flag
     *
     */
    inline bool
    isSendYetGet() const;


    /**
     * \brief Set receive already performed flag.
     *
     * \param [in]  value Flag value
     *
     */
    inline void
    isRecvYetSet
    (
      const bool value
    );


    /**
     * \brief Get receive already performed flag.
     *
     * \return Receive already performed flag
     *
     */
    inline bool
    isRecvYetGet() const;


    /**
     * \brief Get sent variable writer identifier.
     *
     * \param [in]  i_comp Component identifier
     *
     * \return Sent variable writer identifier
     *
     */
    inline int
    idWriterVarSendGet
    (
      const int i_comp
    ) const;


    /**
     * \brief Get sent variable status writer identifier.
     *
     * \return Sent variable status writer identifier
     *
     */
    inline int
    idWriterVarSendStatusGet() const;


    /**
     * \brief Get received variable writer identifiers.
     *
     * \param [in]  i_comp Component identifier
     *
     * \return Received variable writer identifiers
     *
     */
    inline int
    idWriterVarRecvGet
    (
      const int i_comp
    ) const;


    /**
     * \brief Get received variable status writer identifier.
     *
     * \return Received variable status writer identifier
     *
     */
    inline int
    idWriterVarRecvStatusGet() const;


    /**
     * \brief Enable computed target broadcast.
     *
     */
    inline void
    computedTgtBcastEnable();


    /**
     * \brief Check if computed target broadcast is enabled.
     *
     * \return true if enabled, false otherwise
     *
     */
    inline bool
    computedTgtBcastEnabledIs() const;


    /**
     * \brief Enable involved source broadcast.
     *
     */
    inline void
    involvedSrcBcastEnable();


    /**
     * \brief Check if involved source broadcast is enabled.
     *
     * \return true if enabled, false otherwise
     *
     */
    inline bool
    involvedSrcBcastEnabledIs() const;

  private:

    std::string               _field_id;                      /*!< Field identifier */
    int                       _field_id_int;                  /*!< Integer field identifier */
    CWP_Type_t                _data_type;                     /*!< Field data type */
    int                       _data_type_size;                /*!< Field data type size */
    CWP_Dof_location_t        _field_location;                /*!< Value location */
    CWP_Dof_location_t        _linked_field_location;         /*!< Value location of coupled field */
    int                       _n_component;                   /*!< Number of components */
    CWP_Field_storage_t       _storage;                       /*!< Storage type */
    CWP_Field_exch_t          _exchange_type;                 /*!< Exchange type */
    CWP_Status_t              _visu_status;                   /*!< Visualization status */

    Coupling                 *_cpl;                           /*!< Coupling */
    Mesh                     *_mesh;                          /*!< Mesh */
    int                       _n_part;                        /*!< Number of partitions */
    void                    **_data_tgt;                      /*!< Target data array */
    void                    **_data_src;                      /*!< Source data array */

    bool                      _is_send_yet;                   /*!< Tells if a field has been sent at a given moment */
    bool                      _is_recv_yet;                   /*!< Tells if a field has been received at a given moment */

    PDM_writer_t             *_writer;                        /*!< PDM writer */
    int                      *_id_writer_var_send;            /*!< Sent variable writer identifiers */
    int                      *_id_writer_var_recv;            /*!< Received variable writer identifiers */
    int                       _id_writer_var_send_status;     /*!< Sent variable status writer identifier */
    int                       _id_writer_var_recv_status;     /*!< Received variable status writer identifier */

    typedef enum {

      CWP_FIELD_NOT_EXCHANGED              = -1,              /*!< Not exchanged */
      CWP_FIELD_COMPUTED_AND_EXCHANGED     =  0,              /*!< Computed and exchanged */
      CWP_FIELD_NOT_COMPUTED_AND_EXCHANGED =  1               /*!< Not computed and exchanged */

    } CWP_Field_status_t;

    CWP_Interpolation_t       _interpolation_type;            /*!< Interpolation type */
    CWP_Interp_function_t     _interpolation_function;        /*!< Interpolation function */
    CWP_Interp_function_t     _interpolation_function_f;      /*!< Fortran interpolation function */
    CWP_Interp_function_p_t   _interpolation_function_p;      /*!< Python interpolation function */
    void                     *_python_object;                 /*!< Python object */

    bool                      _is_computed_tgt_bcast_enabled; /*!< Computed target broadcast enabled flag */
    bool                      _is_involved_src_bcast_enabled; /*!< Involved source broadcast enabled flag */
  };

}

/**
 * \endcond
 */

#endif // __FIELD_H__