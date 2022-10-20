

#include "mesh.hxx"
#include "field.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "pdm_writer.h"

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
   _mesh (cpl->meshGet()),
   _writer (cpl->writerGet()),
   _id_writer_var_send(nullptr),
   _id_writer_var_recv(nullptr),
   _id_writer_var_recv_computed(-1),
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
    default:
      PDM_error(__FILE__, __LINE__, 0, "CWP_CHAR is not usable.\n");
  }

  _interpolationType = CWP_INTERPOLATION_DEFAULT;

  if (_writer != nullptr) {

    _id_writer_var_send = (int *) malloc (sizeof (int) * _nComponent);
    _id_writer_var_recv = (int *) malloc (sizeof (int) * _nComponent);

    PDM_writer_var_loc_t pdm_field_type = PDM_WRITER_VAR_ELEMENTS;

    if     ( _fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) { 
      pdm_field_type = PDM_WRITER_VAR_ELEMENTS ;
    }
    else if ( _fieldLocation == CWP_DOF_LOCATION_NODE        ) {
      pdm_field_type = PDM_WRITER_VAR_VERTICES;
    }
    else if ( _fieldLocation == CWP_DOF_LOCATION_USER        ) {
      pdm_field_type = PDM_WRITER_VAR_PARTICLES ;
    }


    if (_exchangeType == CWP_FIELD_EXCH_SEND ||
        _exchangeType == CWP_FIELD_EXCH_SENDRECV) {

  Field::~Field(){

     _data_tgt.clear();
     _data_src.clear();
      std::string prefix = "s";

      for (int i_comp = 0; i_comp < nComponent; i_comp++) {

        std::ostringstream num;
        num << (i_comp + 1);

        std::string fieldName = prefix + "_" + _fieldID + num.str();

        _id_writer_var_send[i_comp] = PDM_writer_var_create(_writer,
                                                            PDM_WRITER_ON,
                                                            PDM_WRITER_VAR_SCALAR,
                                                            pdm_field_type,
                                                            fieldName.c_str());

        printf("WriterFieldCreate - send: '%s' %d\n",fieldName.c_str(), _id_writer_var_send[i_comp]);
        fflush(stdout);
      }
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Field::dataSet Error : unknown data type.\n");
    }

    if (_exchangeType == CWP_FIELD_EXCH_RECV ||
        _exchangeType == CWP_FIELD_EXCH_SENDRECV) {

      std::string prefix = "r";

      for (int i_comp = 0; i_comp < nComponent; i_comp++) {

        std::ostringstream num;
        num << (i_comp + 1);

        std::string fieldName = prefix + "_" + _fieldID + num.str();

        _id_writer_var_recv[i_comp] = PDM_writer_var_create(_writer,
                                                            PDM_WRITER_ON,
                                                            PDM_WRITER_VAR_SCALAR,
                                                            pdm_field_type,
                                                            fieldName.c_str());
        printf("WriterFieldCreate - recv: '%s' %d\n",fieldName.c_str(), _id_writer_var_recv[i_comp]);
        fflush(stdout);
      }

      std::string fieldName = prefix + "_" + _fieldID;
      std::string fieldComputedName = fieldName + "_is_computed";

      PDM_writer_status_t time_dependent = PDM_WRITER_ON;

      CWP_Dynamic_mesh_t topology = cpl->DisplacementGet();

      if (topology == CWP_DYNAMIC_MESH_STATIC) {
        time_dependent = PDM_WRITER_OFF;
      }

      _id_writer_var_recv_computed = PDM_writer_var_create(_writer,
                                                           time_dependent,
                                                           PDM_WRITER_VAR_SCALAR,
                                                           pdm_field_type,
                                                           fieldComputedName.c_str());


      printf("WriterFieldCreate - recv 2: '%s' %d\n",fieldComputedName.c_str(), _id_writer_var_recv_computed);
      fflush(stdout);

    }

  }

}



Field::~Field()
{
  _data_tgt.clear();
  _data_src.clear();

  if (_id_writer_var_send != nullptr) {
    free (_id_writer_var_send);
  }
  if (_id_writer_var_recv != nullptr) {
    free (_id_writer_var_recv);
  }

  if (_sendBuffer != NULL) {
    free(_sendBuffer);
  }
  if (_recvBuffer != NULL) {
    free(_recvBuffer);
  }
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
    PDM_error(__FILE__, __LINE__, 0, "Field::dataSet Error : unknoown data type.\n");
  }

}


}

/**
 * \endcond
 */
