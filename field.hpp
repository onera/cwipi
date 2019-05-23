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
#include <map>


namespace cwipi {
  
  /** 
   * \class Field field.hpp "field.hpp"
   * \brief Abstract field
   *
   *  This class is field abstract interface 
   * 
   */
  class Mesh;
  class Visu;
  template <typename dataType> 
  class Field {
    
  public:
    
    /**
     * \brief Constructor
     *
     */
    Field() {}

    Field (std::string            field_id    ,
           Mesh*                  mesh        ,
           CWP_Field_value_t      fieldType   ,
           CWP_Field_storage_t    storage     ,
           int                    nComponent  ,
           CWP_Field_exch_t       exchangeType,
           CWP_Status_t           visuStatus  ,
           int*                   iteration   ,
           double*                physTime    ): 
           
           _fieldID        (field_id)    ,
           _mesh           (mesh)        ,
           _storage        (storage)     ,
           _nComponent     (nComponent)  ,
           _fieldLocation  (fieldType)   ,
           _exchangeType   (exchangeType),
           _visuStatus     (visuStatus)  ,
           _iteration      (iteration)   ,
           _physTime       (physTime)
    {
       _n_part = _mesh -> getNPart(); 
       _data.resize(_n_part,NULL);
       _sendBuffer = NULL;
       _recvBuffer = NULL;  
       _last_request.resize(300,0);
    }

    /**
     * \brief Destructor
     *
     */

    ~Field(){
       if (_sendBuffer != NULL) free(_sendBuffer);
       if (_recvBuffer != NULL) free(_recvBuffer);           
     }

    /**
     * \brief set data array
     *
     * \param [in] data   Data array
     *
     */

    void dataSet ( int i_part, double data[])
    {
      int size = _mesh->getPartNElts(i_part);
      _data[i_part] = &(data[0]);
  /*    _data[i_part] = new double[size];
      for (int i =0;i<size;i++)
        _data[i_part][i] = data[i]; */
    }

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

    /**
     *
     * \brief Get field nature
     * 
     * \return           Field nature
     * 
     */

    inline CWP_Field_value_t
    typeGet() const
    {
      return _fieldLocation;
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

    dataType* dataGet(int i_part) const
    { 
      return _data[i_part];
    }

    void visuIdSet(int visu_id)
    { 
      _visu_id = visu_id;
    }



    int visuIdGet() const
    { 
      return _visu_id;
    }



    std::vector<dataType*> dataGetAll() const
    {
      return _data;
    }


    int* iterationGet() const
    {
      return _iteration;
    }
    
    double* physicalTimeGet() const
    {
      return _physTime;
    }
    

    void ReceptionBufferCreation(std::vector<int> nLocatedTargets,int TotLocatedTargets) {
        std::ostringstream strs;
        
        strs <<"interp"<<_fieldID<<"_"<<_iteration;
        std::string fieldID = strs.str();
        int* iteration = new int (*_iteration);
        double* physTime = new double (*_physTime); 
                                                  
        //On alloue l'espace pour la réception si pas déjà fait                                               
        if(_recvBuffer == NULL) {
          _recvBuffer = (double*)malloc(sizeof(double)*_nComponent*TotLocatedTargets);
        }

        std::vector<double*> v_interpolatedData;
        v_interpolatedData.resize(_n_part);
        int loc_ind=-1;
          
        for(int i_part=0;i_part<_n_part;i_part++) {
          loc_ind = _nComponent * nLocatedTargets[i_part];
          v_interpolatedData[i_part] = &(_recvBuffer[loc_ind]);
        }

    }


  void lastRequestAdd (int i_proc, MPI_Request request) {
    _last_request[i_proc] = request;
  }

  
  MPI_Request lastRequestGet (int i_proc) {
    return _last_request[i_proc];
  }


  
  double* recvBufferGet () {
    return _recvBuffer;
  }

  double* sendBufferGet () {
    return _sendBuffer;
  }

  void sendBufferSet (double* sendBuffer) {
    _sendBuffer = sendBuffer;
  }
   

  private:

    const CWP_Field_storage_t                _storage;      /*!< Storage type */ 
    const int                                _nComponent;   /*!< Number of component */
    const CWP_Field_value_t                  _fieldLocation;    /*!< Value location */
    const CWP_Field_exch_t                   _exchangeType; /*!< Exchange type */
    const CWP_Status_t                       _visuStatus;   /*!< Visualization status */
    std::vector<dataType* >                  _data;         /*!< Pointer to data array */
    std::string                              _fieldID;
    double                                  *_sendBuffer;
    double                                  *_recvBuffer;
    double*                                  _physTime;
    int*                                     _iteration;
    Mesh                                    *_mesh;
    int                                      _n_part;
    int                                      _visu_id;
    std::vector <MPI_Request>                _last_request;

    Field &operator=(const Field &other);       /*!< Assigment operator not available */
    Field (const Field& other);                 /*!< Copy constructor not available */
  };
  
  
  
  
  

}

#endif //__FIELD_H__
