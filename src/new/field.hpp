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

using namespace std;

namespace cwipi {
  
  /** 
   * \class Field field.hpp "field.hpp"
   * \brief Abstract field
   *
   *  This class is field abstract interface 
   * 
   */
  
  template <typename DataType> 
  class Field {
    
  public:
    
    /**
     * \brief Constructor
     *
     */

    Field
    (
     CWP_Type_t     dataType,
     CWP_Field_storage_t  storage,
     int                    nComponent,
     CWP_Field_value_t   nature,
     CWP_Field_exch_t     exchangeType,
     CWP_Status_t         visuStatus
     ): _dataType(dataType),
        _storage(storage),
        _nComponent(nComponent),
        _nature(nature),
        _exchangeType(exchangeType),
        _visuStatus(visuStatus)
    {
      _data = NULL;
    }

    /**
     * \brief Destructor
     *
     */

    virtual ~Field(){}

    /**
     * \brief set data array
     *
     * \param [in] data   Data array
     *
     */

    inline void 
    mappingSet
    (
     DataType        data[]
    )
    {
      _data = data;
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

    /**
     *
     * \brief Get data type
     * 
     * \return            Field storage type
     * 
     */

    inline CWP_Type_t
    dataTypeGet() const
    {
      return _dataType;
    }

    /**
     *
     * \brief Get field nature
     * 
     * \return           Field nature
     * 
     */

    inline CWP_Field_value_t
    natureGet() const
    {
      return _nature;
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

    inline DataType *
    dataGet() const
    {
      return _data;
    }

  private:

    const CWP_Field_storage_t  _storage;      /*!< Storage type */ 
    const int                    _nComponent;   /*!< Number of component */
    const CWP_Type_t     _dataType;    /*!< Data type */
    const CWP_Field_value_t   _nature;       /*!< Nature */
    const CWP_Field_exch_t     _exchangeType; /*!< Exchange type */
    const CWP_Status_t         _visuStatus;   /*!< Visualization status */
    DataType                     *_data;       /*!< Pointer to data array */

    Field &operator=(const Field &other);       /*!< Assigment operator not available */
    Field (const Field& other);                 /*!< Copy constructor not available */
  };

}

#endif //__FIELD_H__
