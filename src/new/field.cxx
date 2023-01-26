

#include "mesh.hxx"
#include "field.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

#include <coupling_i.hxx>

/**
 * \cond
 */

namespace cwipi {

   Field::Field (std::string            field_id    ,
                 int                     fieldIDInt,      
                 CWP_Type_t             dataType    ,
                 Coupling*              cpl        ,
                 CWP_Dof_location_t      fieldType   ,
                 CWP_Field_storage_t    storage     ,
                 int                    nComponent  ,
                 CWP_Field_exch_t       exchangeType,
                 CWP_Status_t           visuStatus):
   _storage        (storage)     ,
   _nComponent     (nComponent)  ,
   _fieldLocation  (fieldType)   ,
   _linkedFieldLocation (CWP_DOF_LOCATION_UNDEF),
   _exchangeType   (exchangeType),
   _visuStatus     (visuStatus)  ,
   _dataType       (dataType)    ,
   _fieldID        (field_id)    ,
   _fieldIDInt     (fieldIDInt)  ,
   _cpl            (cpl)         ,
   _interpolationFunction   (NULL),
   _interpolationFunction_f (NULL),
   _current_step_was_exchanged(0)

  {
    _mesh = cpl->meshGet();
    _n_part = _mesh->getNPart();
    _data_tgt.resize(_n_part,NULL);
    _data_src.resize(_n_part,NULL);
    _sendBuffer = NULL;
    _recvBuffer = NULL;
    _dataTypeSize = 0;

    switch (_dataType) {
      case CWP_DOUBLE:
        _dataTypeSize = sizeof(double);
        break;
      case CWP_INT:
        _dataTypeSize = sizeof(int);
        break;
      case CWP_CHAR:
        PDM_error(__FILE__, __LINE__, 0, "CWP_CHAR is not usable.\n");

    }

     _interpolationType = CWP_INTERPOLATION_DEFAULT;
  }


  Field::~Field(){

     if (_data_tgt.size() != 0) _data_tgt.clear();
     if (_data_src.size() != 0) _data_src.clear();

     if (_sendBuffer != NULL) free(_sendBuffer);
     _sendBuffer = NULL;
     if (_recvBuffer != NULL) free(_recvBuffer);
     _recvBuffer = NULL;
   }


  void Field::dataSet ( int i_part, const CWP_Field_map_t   map_type, void* data)
  {
    if (map_type == CWP_FIELD_MAP_SOURCE) {
      _data_src[i_part] = data;
    }
    else if (map_type == CWP_FIELD_MAP_TARGET) {
      _data_tgt[i_part] = data;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Field::dataSet Error : unknown data type.\n");
    }

  }

}

/**
 * \endcond
 */
