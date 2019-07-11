

#include <mesh.hxx>
#include <field.hxx>

namespace cwipi {

   Field::Field (std::string            field_id    ,
           CWP_Type_t             dataType    ,
           Coupling*              cpl        ,
           CWP_Field_value_t      fieldType   ,
           CWP_Field_storage_t    storage     ,
           int                    nComponent  ,
           CWP_Field_exch_t       exchangeType,
           CWP_Status_t           visuStatus  ,
           int*                   iteration   ,
           double*                physTime    ): 
           
           _fieldID        (field_id)    ,
           _dataType       (dataType)    ,
           _storage        (storage)     ,
           _nComponent     (nComponent)  ,
           _fieldLocation  (fieldType)   ,
           _exchangeType   (exchangeType),
           _visuStatus     (visuStatus)  ,
           _iteration      (iteration)   ,
           _physTime       (physTime)
    {
    
       _mesh = cpl -> meshGet(); 
       _n_part = _mesh -> getNPart();
       _data.resize(_n_part,NULL);
       _sendBuffer = NULL;
       _recvBuffer = NULL;  
        
       int len = _fieldID.length();
       int id=0;
       for(int i=0;i<len;i++) {
         id+=(int)_fieldID[i];
       }
       _fieldIDInt = id;
       
      _dataTypeSize = 0;
      switch (_dataType) {
        case CWP_DOUBLE:
          _dataTypeSize = sizeof(double);
          break;
        case CWP_INT:
          _dataTypeSize = sizeof(int);
          break;         
        
      }      
    }
    
    
    
    Field::~Field(){
       if (_sendBuffer != NULL) free(_sendBuffer);
       if (_recvBuffer != NULL) free(_recvBuffer);           
     }

  void Field::dataSet ( int i_part, void* data)
  {
    int size = _mesh->getPartNElts(i_part);
    _data[i_part] = data;
  }
    

}
