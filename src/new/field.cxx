

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
   _writer (nullptr),
   _id_writer_var_send(nullptr),
   _id_writer_var_recv(nullptr),
   _id_writer_var_recv_computed(-1),
   _interpolationFunction   (NULL),
   _current_step_was_exchanged(0),
   _computed_tgt_bcast_enabled(0),
   _involved_src_bcast_enabled(0)

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


  if (visuStatus == CWP_STATUS_ON) {
    _writer = cpl->writerGet();
  }

  _interpolationType = CWP_INTERPOLATION_DEFAULT;

  if (_writer != nullptr) {

    _id_writer_var_send = (int *) malloc (sizeof (int) * _nComponent);
    _id_writer_var_recv = (int *) malloc (sizeof (int) * _nComponent);

    // Set field location type
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

    // Create send variables
    if (_exchangeType == CWP_FIELD_EXCH_SEND ||
        _exchangeType == CWP_FIELD_EXCH_SENDRECV) {

      std::string prefix = "s";

      for (int i_comp = 0; i_comp < _nComponent; i_comp++) {

        std::ostringstream num;
        num << (i_comp + 1);

        std::string fieldName = prefix + "_" + _fieldID + num.str();

        _id_writer_var_send[i_comp] = PDM_writer_var_create(_writer,
                                                            PDM_WRITER_ON,
                                                            PDM_WRITER_VAR_SCALAR,
                                                            pdm_field_type,
                                                            fieldName.c_str());

        // printf("WriterFieldCreate - send: '%s' %d\n",fieldName.c_str(), _id_writer_var_send[i_comp]);
        // fflush(stdout);
      }
    }

    // Create receive variables
    if (_exchangeType == CWP_FIELD_EXCH_RECV ||
        _exchangeType == CWP_FIELD_EXCH_SENDRECV) {

      std::string prefix = "r";

      for (int i_comp = 0; i_comp < _nComponent; i_comp++) {

        std::ostringstream num;
        num << (i_comp + 1);

        std::string fieldName = prefix + "_" + _fieldID + num.str();

        _id_writer_var_recv[i_comp] = PDM_writer_var_create(_writer,
                                                            PDM_WRITER_ON,
                                                            PDM_WRITER_VAR_SCALAR,
                                                            pdm_field_type,
                                                            fieldName.c_str());
        // printf("WriterFieldCreate - recv: '%s' %d\n",fieldName.c_str(), _id_writer_var_recv[i_comp]);
        // fflush(stdout);
      }

      // Create variable to tag if dof has been located ("computed")
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


      // printf("WriterFieldCreate - recv 2: '%s' %d\n",fieldComputedName.c_str(), _id_writer_var_recv_computed);
      // fflush(stdout);

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


/**
 * \brief Write the field send or receive variables
 *
 * \param [in] exch_type Exchange type
 *
 */

void 
Field::write 
(
CWP_Field_exch_t exch_type
)
{
  assert(CWP_FIELD_EXCH_SENDRECV != exch_type);

  // Get data depending on exchange type
  int *id_writer_var;
  std::vector<void* > &data_var = _data_src;

  if (exch_type == CWP_FIELD_EXCH_SEND) {
    id_writer_var = _id_writer_var_send;
    data_var = _data_src;
  }
  else {
    id_writer_var = _id_writer_var_recv;  
    data_var = _data_tgt;
  }

  // Write
  if (_writer != nullptr) {
    if ((_cpl->NStepGet() % _cpl->freqWriterGet()) == 0) {

      double **comp_data = (double **) malloc (sizeof(double *) * _cpl->nPartGet());

      // (x1, y1, z1, ... , xn, yn, zn) or recv
      // Allocate memory for reordering to write per variable
      if ((_storage == CWP_FIELD_STORAGE_INTERLACED && _nComponent > 1) ||
          (exch_type == CWP_FIELD_EXCH_RECV)) {
        for(int i_part= 0 ; i_part < _cpl->nPartGet(); i_part++){
          int n_elt_size = 0;

          if (_fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            n_elt_size = _mesh->getPartNElts(i_part);
          }
          else if (_fieldLocation == CWP_DOF_LOCATION_NODE) {
            n_elt_size = _mesh->getPartNVertex(i_part);
          }
          else if (_fieldLocation == CWP_DOF_LOCATION_USER) {
            n_elt_size = _cpl->userTargetNGet(i_part);
          }

          comp_data[i_part] = (double *) malloc (sizeof(double) * n_elt_size);

        } // end loop n_part
      }

      // Fill in array and write it
      for (int i_comp = 0; i_comp < _nComponent; i_comp++) {

        for(int i_part= 0 ; i_part < _cpl->nPartGet(); i_part++){

          // size
          int n_elt_size = 0;

          if (_fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            n_elt_size = _mesh->getPartNElts(i_part);
          }
          else if (_fieldLocation == CWP_DOF_LOCATION_NODE) {
            n_elt_size = _mesh->getPartNVertex(i_part);
          }
          else if (_fieldLocation == CWP_DOF_LOCATION_USER) {
            n_elt_size = _cpl->userTargetNGet(i_part);
          }

          // recv computed
          const int  *computed_target    = NULL;
          int         n_computed_target  = 0;
          if (exch_type == CWP_FIELD_EXCH_RECV) {
            if (_cpl->has_mesh()) {
              computed_target   = _cpl->computedTargetsGet (_fieldID, i_part);
              n_computed_target = _cpl->nComputedTargetsGet(_fieldID, i_part);
            }
          }

          // (x1, y1, z1, ... , xn, yn, zn)
          if (_storage == CWP_FIELD_STORAGE_INTERLACED && _nComponent > 1) {
            if ((exch_type == CWP_FIELD_EXCH_RECV) && (n_computed_target != n_elt_size)) {
              for (int j = 0; j < n_elt_size; j++) {
                comp_data[i_part][j] = HUGE_VAL;
              }
              for (int j = 0; j < n_computed_target; j++) {
                int jdx = computed_target[j] -1;
                comp_data[i_part][jdx] = ((double *) data_var[i_part])[j * _nComponent + i_comp];
              }
            } else {
              for (int j = 0; j < n_elt_size; j++) {
                comp_data[i_part][j] = ((double *) data_var[i_part])[j * _nComponent + i_comp];
              }
            }
          }
          // (x1, ... xn, y1, ..., yn, z1, ...zn)
          else {
            if (exch_type == CWP_FIELD_EXCH_RECV) {
              for (int j = 0; j < n_elt_size; j++) {
                comp_data[i_part][j] = HUGE_VAL;
              }
              for (int j = 0; j < n_computed_target; j++) {
                int jdx = computed_target[j] -1;
                comp_data[i_part][jdx] = ((double *) data_var[i_part])[n_computed_target * i_comp + j];
              }
            } else {
              comp_data[i_part] = (double *) &(((double *)data_var[i_part])[n_elt_size * i_comp]);
            }
          }

          PDM_writer_var_set(_writer, id_writer_var[i_comp], _cpl->idGeomWriterGet(), i_part, comp_data[i_part]);
        } // end loop n_part

        PDM_writer_var_write(_writer, id_writer_var[i_comp]);

        PDM_writer_var_data_free(_writer, id_writer_var[i_comp]);

      } // end loop on components

      // write computed tag
      if (exch_type == CWP_FIELD_EXCH_RECV) {

        CWP_Dynamic_mesh_t topology = _cpl->DisplacementGet();

        if (topology != CWP_DYNAMIC_MESH_STATIC || _cpl->NStepGet() == 0) {

          for(int i_part= 0 ; i_part < _cpl->nPartGet(); i_part++){

            // size
            int n_elt_size = 0;

            if (_fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              n_elt_size = _mesh->getPartNElts(i_part);
            }
            else if (_fieldLocation == CWP_DOF_LOCATION_NODE) {
              n_elt_size = _mesh->getPartNVertex(i_part);
            }
            else if (_fieldLocation == CWP_DOF_LOCATION_USER) {
              n_elt_size = _cpl->userTargetNGet(i_part);
            }

            // computed
            const int  *computed_target    = NULL;
            int         n_computed_target  = 0;

            if (_cpl->has_mesh()) {
              computed_target   = _cpl->computedTargetsGet (_fieldID, i_part);
              n_computed_target = _cpl->nComputedTargetsGet(_fieldID, i_part);
            }

            if (n_elt_size == n_computed_target) {
              for (int i = 0; i < n_elt_size; i++) {
                comp_data[i_part][i] = 1.;
              }
            } else {
              for (int i = 0; i < n_elt_size; i++) {
                comp_data[i_part][i] = 0.;
              }
              for (int i = 0; i < n_computed_target; i++) {
                int idx = computed_target[i] -1;
                comp_data[i_part][idx] = 1.;
              }
            }

            PDM_writer_var_set(_writer,
                               _id_writer_var_recv_computed,
                               _cpl->idGeomWriterGet(),
                               i_part,
                               comp_data[i_part]);

          } // end loop n_part

          PDM_writer_var_write(_writer, _id_writer_var_recv_computed);

          PDM_writer_var_data_free(_writer, _id_writer_var_recv_computed);

        }
      }

      if ((_storage == CWP_FIELD_STORAGE_INTERLACED && _nComponent > 1) ||
          (exch_type == CWP_FIELD_EXCH_RECV)) {
        for(int i_part= 0 ; i_part < _cpl->nPartGet(); i_part++){
          free (comp_data[i_part]);
        }
      }
      free (comp_data);

    } // end if frequency
  } // end if writer not null

}


}

/**
 * \endcond
 */
