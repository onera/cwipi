#ifndef __FIELD_H__
#define __FIELD_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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

#include <sstream>
#include <mesh.hxx>
#include <coupling.hxx>
#include "cwp.h"
#include "cwp_priv.h"
#include <map>


/**
 * \cond
 */

namespace cwipi {

  /**
   * \class Field field.hxx "field.hxx"
   * \brief Abstract field
   *
   *  This class is field abstract interface
   *
   */
  class Mesh;
  class Field {

  public:

    /**
     * \brief Constructor
     *
     */

    Field() {}

    Field (std::string            field_id    ,
           int                    fieldIDInt,      
           CWP_Type_t             dataType    ,
           Coupling*              cpl         ,
           CWP_Dof_location_t      fieldType   ,
           CWP_Field_storage_t    storage     ,
           int                    nComponent  ,
           CWP_Field_exch_t       exchangeType,
           CWP_Status_t           visuStatus  ,
           int*                   iteration   ,
           double*                physTime    );

    /**
     * \brief Destructor
     *
     */

    ~Field();

    /**
     * \brief set data array
     *
     * \param [in] data   Data array
     *
     */

    void dataSet (int i_part, const CWP_Field_map_t  map_type, void* data);


    /**
     *
     * \brief Get field storage type
     *
     * \return            Field storage type
     *
     */

    inline CWP_Field_storage_t
    storageTypeGet() const
    {
      return _storage;
    }

    /**
     *
     * \brief Get nunmber of field components
     *
     * \return             Number of field components
     *
     */

    inline int
    nComponentGet() const
    {
      return _nComponent;
    }

    inline std::string
    fieldIDGet() const
    {
      return _fieldID;
    }


    inline int
    fieldIDIntGet() const
    {
      return _fieldIDInt;
    }

    /**
     *
     * \brief Get field nature
     *
     * \return           Field nature
     *
     */

    inline CWP_Dof_location_t
    locationGet() const
    {
      return _fieldLocation;
    }

    CWP_Type_t
    dataTypeGet() const
    {
      return _dataType;
    }

    int
    dataTypeSizeGet() const
    {
      return _dataTypeSize;
    }

    /**
     *
     * \brief Get exchange type
     *
     * \return          Exchange type
     *
     */

    inline CWP_Field_exch_t
    exchangeTypeGet() const
    {
      return _exchangeType;
    }

    /**
     *
     * \brief Get visu status
     *
     * \return          Exchange type
     *
     */

    inline CWP_Status_t
    visuStatusGet()  const
    {
      return _visuStatus;
    }

    /**
     *
     * \brief Get data
     *
     * \return          Data
     *
     */

    void* dataGet(int i_part,  const CWP_Field_map_t  map_type) const
    {
      if (map_type == CWP_FIELD_MAP_SOURCE) {
        return _data_src[i_part];
      }
      else if (map_type == CWP_FIELD_MAP_TARGET) {
        return _data_tgt[i_part];
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "Field::dataSet Error : unknoown data type.\n");
        return nullptr;
      }
    }

    void visuIdSet(int visu_id)
    {
      _visu_id = visu_id;
    }

    void visuIdComputedSet(int visu_id_computed)
    {
      _visu_id_computed = visu_id_computed;
    }

    void interpFromLocationSet(CWP_Interp_from_location_t fct)
    {
      _interpolationType     = CWP_INTERPOLATION_USER ;
      _interpolationFunction = fct                    ;
    }


    CWP_Interp_from_location_t interpolationFunctionGet() {
      return _interpolationFunction;
    }

    CWP_Interpolation_t interpolationTypeGet() {
      return _interpolationType;
    }

    int visuIdGet() const
    {
      return _visu_id;
    }

    int visuIdComputedGet() const
    {
      return _visu_id_computed;
    }



    // std::vector<void*> dataGetAll() const
    // {
    //   return _data;
    // }


    int* iterationGet() const
    {
      return _iteration;
    }

    double* physicalTimeGet() const
    {
      return _physTime;
    }

    void currentStepWasExchangedReset()
    {
      _current_step_was_exchanged = 0;
    }

    int currentStepWasExchangedGet() const
    {
      return _current_step_was_exchanged;
    }


    void ReceptionBufferCreation(int TotLocatedTargets) {
        std::ostringstream strs;

        strs <<"interp"<<_fieldID<<"_"<<_iteration;
        std::string fieldID = strs.str();

        //On alloue l'espace pour la réception si pas déjà fait
        if(_recvBuffer == NULL) {
          _recvBuffer = (void*)malloc(_dataTypeSize*_nComponent*TotLocatedTargets);
        }
    }


  void lastRequestAdd (int i_proc, MPI_Request request) {
    _last_request[i_proc] = request;
  }


  MPI_Request lastRequestGet (int i_proc) {
    return _last_request[i_proc];

  }


  void lastRequestAdd_p2p (int i_proc, std::vector<MPI_Request> request) {
    _last_request_p2p[i_proc] = request;
  }


  std::vector<MPI_Request> lastRequestGet_p2p (int i_proc) {
    return _last_request_p2p[i_proc];
  }



  void* recvBufferGet () {
    return _recvBuffer;
  }

  void* sendBufferGet () {
    return _sendBuffer;
  }

  void sendBufferSet (void* sendBuffer) {
    _sendBuffer = sendBuffer;
  }

  void linkedFieldLocationSet(CWP_Dof_location_t linkedFieldLocation){
    _linkedFieldLocation = linkedFieldLocation;
  }

  CWP_Dof_location_t linkedFieldLocationGet(){
    return _linkedFieldLocation;
  }

  Coupling *couplingGet(){
    return _cpl;
  }

  Mesh *meshGet(){
    return _mesh;
  }

  private:

    CWP_Field_storage_t                      _storage;        /*!< Storage type */
    int                                      _nComponent;     /*!< Number of component */
    CWP_Dof_location_t                       _fieldLocation;  /*!< Value location Interpolation methods for sender and cloud points type for receiver */
    CWP_Dof_location_t                       _linkedFieldLocation; /*!< Value location Interpolation methods for sender and cloud points type for receiver */
    CWP_Field_exch_t                         _exchangeType;   /*!< Exchange type */
    CWP_Status_t                             _visuStatus;     /*!< Visualization status */
    std::vector<void* >                      _data_src;       /*!< Pointer to data array Add a another data poiter for send fields */
    std::vector<void* >                      _data_tgt;       /*!< Pointer to data array Add a another data poiter for recv fields */
    CWP_Type_t                               _dataType;
    std::string                              _fieldID;
    int                                      _fieldIDInt;
    void                                    *_sendBuffer;
    void                                    *_recvBuffer;
    double*                                  _physTime;
    int*                                     _iteration;
    Coupling                                *_cpl;
    Mesh                                    *_mesh;
    int                                      _n_part;
    int                                      _visu_id;
    int                                      _visu_id_computed;
    std::map <int,MPI_Request>               _last_request;
    std::map <int,std::vector<MPI_Request>>  _last_request_p2p;
    int                                      _dataTypeSize;
    CWP_Interp_from_location_t               _interpolationFunction;
    CWP_Interpolation_t                      _interpolationType;

    Field &operator=(const Field &other);       /*!< Assigment operator not available */
    Field (const Field& other);                 /*!< Copy constructor not available */

    int                                      _current_step_was_exchanged;
  };






}

/**
 * \endcond
 */

#endif //__FIELD_H__
