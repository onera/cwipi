/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-2017  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more detailstr_options.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "visualization.hxx"
#include <string>
#include <iostream>
#include "pdm_writer.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include <unistd.h>
#include <stdio.h>
#include "field.hxx"
#include "cwp.h"
#include "cwp_priv.h"
#include "coupling.hxx"
#include "coupling_i.hxx"


/**
 * \cond
 */

namespace cwipi {

  Visu::Visu(const MPI_Comm &MPIComm,const CWP_Dynamic_mesh_t topology):
                                      _visu_id(nullptr),_visu_mesh_id(-1),
                                      //_freq(-1), 
                                      _output_dir(NULL),
                                      _output_name(NULL),
                                      _divide_polygons(PDM_WRITER_OFF),
                                      _divide_polyhedra(PDM_WRITER_OFF),
                                      _visuCreated(false),
                                      _physical_time(-1),
                                      _topology(topology) {

     _pdmComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&MPIComm));

  }

  Visu::~Visu() {
    // Check if these structures have been filled (ie geometry declared)

    // TODO :  See with Eric but if we call multiple time GeomWrite a leaks appear
    // if (!_partitioning_field_data.empty()) {
    //   for(int i_part= 0 ; i_part < _n_part; i_part++) {
    //     free(_partitioning_field_data[i_part]);
    //   }
    // }
    // if (!_ranking_field_data.empty()) {
    //   for(int i_part= 0 ; i_part < _n_part; i_part++) {
    //     free(_ranking_field_data[i_part]);
    //   }
    // }
    // if (!_blocking_field_data.empty()) {
    //   for(int i_part= 0 ; i_part < _n_part; i_part++) {
    //     free(_blocking_field_data[i_part]);
    //   }
    // }

    // PDM_writer_var_data_free(_visu_id, _id_partitioning_field);
    // PDM_writer_var_data_free(_visu_id, _id_ranking_field);
    // PDM_writer_var_data_free(_visu_id, _id_blocking_field);


    if (_visu_id != nullptr) {
      PDM_writer_step_end(_visu_id);
      PDM_writer_free(_visu_id);
    }
  }


  void Visu::VisuCreate(const int          freq,
                        CWP_Visu_format_t  format,
                        const char        *format_option,
                        char              *output_dir,
                        char              *output_name) {

    CWP_UNUSED(freq);
    CWP_UNUSED(format);

    PDM_writer_fmt_fic_t fmt_fic      = PDM_WRITER_FMT_BIN;
    const char* fmt                   = "Ensight";
    PDM_writer_status_t st_reprise    = PDM_WRITER_OFF;
    const char *options_comp          = "";
    //Proportion of working node for file acess
    int working_node = 1;
    PDM_io_kind_t acess_type =   PDM_IO_KIND_MPI_SIMPLE;

    PDM_writer_topology_t pdm_topology = PDM_WRITER_TOPO_CST;

    if     ( _topology == CWP_DYNAMIC_MESH_STATIC     ) pdm_topology  = PDM_WRITER_TOPO_CST;
    else if( _topology == CWP_DYNAMIC_MESH_DEFORMABLE ) pdm_topology  = PDM_WRITER_TOPO_DEFORMABLE;
    else if( _topology == CWP_DYNAMIC_MESH_VARIABLE   ) pdm_topology  = PDM_WRITER_TOPO_VARIABLE;

    _output_dir  = output_dir;
    _output_name = output_name;

    std::string str_options = format_option;
    std::string delimiter   = ",";

    std::string chars = "\t\n\v\f\r ";

    size_t pos = 0;
    std::string option;
    do
    {
      pos = str_options.find(delimiter);
      option = str_options.substr(0, pos);
      option.erase(0, option.find_first_not_of(chars));
      option.erase(option.find_last_not_of(chars)+1,option.length());
      str_options.erase(0, pos + delimiter.length());

      if (option == "text")                  fmt_fic = PDM_WRITER_FMT_ASCII;
      else if (option == "binary")           fmt_fic = PDM_WRITER_FMT_BIN;
      else if (option == "divide_polygons")  _divide_polygons = PDM_WRITER_ON;
      else if (option == "divide_polyhedra") _divide_polyhedra = PDM_WRITER_ON;
      else if (option != "")                 PDM_error(__FILE__, __LINE__, 0,
                                                        "Not a valid visualization option.\n");

    }
    while (pos != std::string::npos);

    _visu_id = PDM_writer_create(fmt,
                                 fmt_fic,
                                 pdm_topology,
                                 st_reprise,
                                 _output_dir,
                                 _output_name,
                                 _pdmComm,
                                 acess_type,
                                 working_node,
                                 options_comp);

  }


/*****************************************/
  void Visu::GeomCreate(int n_part) {
     _n_part = n_part;
     _visu_mesh_id = PDM_writer_geom_create(_visu_id,
                                            "geom",
                                             n_part);
     _visuCreated = true;
  }

/*****************************************/
  void Visu::GeomCoordSet(int id_part,
                          int n_pts,
                          double *coords,
                          CWP_g_num_t *global_num) {

    PDM_writer_geom_coord_set(_visu_id,
                              _visu_mesh_id,
                              id_part,
                              n_pts,
                              coords,
                              global_num,
                              PDM_OWNERSHIP_USER);
   }

/*****************************************/

  int Visu::GeomBlockAdd(CWP_Block_t blockType) {

     int id_block = PDM_writer_geom_bloc_add(_visu_id,
                               _visu_mesh_id,
                               PdmWriterBlockTypeFromCwpBlockType(blockType),
                               PDM_OWNERSHIP_USER
                              );
    return id_block;
  }

/*****************************************/

  void Visu::GeomWrite(Mesh* mesh) {
    PDM_writer_geom_write(_visu_id,_visu_mesh_id);

    PDM_writer_var_loc_t PDMfieldType = PDM_WRITER_VAR_ELEMENTS   ;

    PDM_writer_var_dim_t PDMfieldComp = PDM_WRITER_VAR_SCALAR;

    int id_partitioning_field = PDM_writer_var_create(_visu_id,
                                                       PDM_WRITER_OFF,
                                                       PDMfieldComp,
                                                       PDMfieldType,
                                                       "partitioning");

    int id_ranking_field = PDM_writer_var_create(_visu_id,
                                                  PDM_WRITER_OFF,
                                                  PDMfieldComp,
                                                  PDMfieldType,
                                                  "ranking");

    int id_blocking_field = PDM_writer_var_create(_visu_id,
                                                  PDM_WRITER_OFF,
                                                  PDMfieldComp,
                                                  PDMfieldType,
                                                  "blocking");

    _partitioning_field_data.resize(mesh->getNPart());
    _ranking_field_data     .resize(mesh->getNPart());
    _blocking_field_data    .resize(mesh->getNPart());

    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);

//    int* blockDB = mesh->blockDBGet();

    int nBlock = mesh->nBlockGet();
    for(int i_part= 0 ; i_part < _n_part; i_part++){
      _partitioning_field_data[i_part] = (double*) malloc( mesh->getPartNElts(i_part) * sizeof(double) );
      _ranking_field_data     [i_part] = (double*) malloc( mesh->getPartNElts(i_part) * sizeof(double) );
      _blocking_field_data    [i_part] = (double*) malloc( mesh->getPartNElts(i_part) * sizeof(double) );
      int idx=0;
      for(int i_block=0;i_block < nBlock; i_block++){
 //       int id_block = blockDB[i_block];
        for(int i_elt=0; i_elt< mesh->getBlockNElts(i_block,i_part); i_elt++){
          _blocking_field_data[i_part][idx++] = (double)i_block;
        }
      }


      for(int i_elt=0; i_elt<mesh->getPartNElts(i_part); i_elt++){
        _partitioning_field_data[i_part][i_elt] = (double)i_part;
        _ranking_field_data[i_part][i_elt] = (double)worldRank;
      }

      PDM_writer_var_set(_visu_id, id_partitioning_field, _visu_mesh_id, i_part, (double *) _partitioning_field_data[i_part]);
      PDM_writer_var_set(_visu_id, id_ranking_field     , _visu_mesh_id, i_part, (double *) _ranking_field_data     [i_part]);
      PDM_writer_var_set(_visu_id, id_blocking_field    , _visu_mesh_id, i_part, (double *) _blocking_field_data    [i_part]);
    }

    PDM_writer_var_write(_visu_id, id_partitioning_field);
    PDM_writer_var_write(_visu_id, id_ranking_field);
    PDM_writer_var_write(_visu_id, id_blocking_field);

    PDM_writer_var_data_free(_visu_id, id_partitioning_field);
    PDM_writer_var_data_free(_visu_id, id_ranking_field);
    PDM_writer_var_data_free(_visu_id, id_blocking_field);

    /* No longer use and recall function causes leaks */
    for(int i_part= 0 ; i_part < _n_part; i_part++){
      free(_partitioning_field_data[i_part]);
      free(_ranking_field_data     [i_part]);
      free(_blocking_field_data    [i_part]);
    }

  }

/*****************************************/

  void Visu::GeomBlockStdSet (int id_block,
                              int id_part,
                              int n_elt,
                              int *connec,
                              CWP_g_num_t *global_num) {
      PDM_writer_geom_bloc_std_set(_visu_id,
                                   _visu_mesh_id,
                                   id_block,
                                   id_part,
                                   n_elt,
                                   connec,
                                   global_num);

  }

/*****************************************/

  void Visu::GeomBlockPoly2D(int id_block,
                             int id_part,
                             int n_elt,
                             int *connec_idx,
                             int *connec,
                             CWP_g_num_t *global_num) {

      PDM_writer_geom_bloc_poly2d_set(_visu_id,
                                      _visu_mesh_id,
                                      id_block,
                                      id_part,
                                      n_elt,
                                      connec_idx,
                                      connec,
                                      global_num);
  }

/**************************************************/

  void Visu::GeomBlockPoly3D(int         id_block,
                             int         id_part,
                             int         n_elts,
                             int         n_faces,
                             int         connec_faces_idx[],
                             int         connec_faces[],
                             int         connec_cells_idx[],
                             int         connec_cells[],
                             CWP_g_num_t global_num[]) {

      PDM_writer_geom_bloc_poly3d_set(_visu_id,
                                      _visu_mesh_id,
                                      id_block,
                                      id_part,
                                      n_elts,
                                      n_faces,
                                      connec_faces_idx,
                                      connec_faces,
                                      connec_cells_idx,
                                      connec_cells,
                                      global_num);
  }

  void Visu::GeomFree() {
    PDM_writer_geom_data_free(_visu_id,_visu_mesh_id);
//    PDM_writer_geom_free(_visu_id,_visu_mesh_id);
  }


  void Visu::WriterFieldCreate(Field* field) {



      CWP_Dof_location_t CWPfielType = field->locationGet();
      int nComponent = field->nComponentGet();
      PDM_writer_var_loc_t PDMfieldType = PDM_WRITER_VAR_ELEMENTS;

      if     ( CWPfielType == CWP_DOF_LOCATION_CELL_CENTER ) PDMfieldType = PDM_WRITER_VAR_ELEMENTS   ;
      else if( CWPfielType == CWP_DOF_LOCATION_NODE        ) PDMfieldType = PDM_WRITER_VAR_VERTICES   ;
      else if( CWPfielType == CWP_DOF_LOCATION_USER        ) PDMfieldType = PDM_WRITER_VAR_PARTICLES ;


      PDM_writer_var_dim_t PDMfieldComp = PDM_WRITER_VAR_SCALAR;
      if( nComponent == 1) PDMfieldComp = PDM_WRITER_VAR_SCALAR;
      else if( nComponent == 3) PDMfieldComp = PDM_WRITER_VAR_VECTOR;
      else if( nComponent !=0 ) {
        PDM_error(__FILE__, __LINE__, 0, "This field have a number of components which cannot be visualized.\n");
      }

      // std::string prefix;
      // if(field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND) {
      //   prefix = "s";
      // } else if(field->exchangeTypeGet() == CWP_FIELD_EXCH_RECV) {
      //   prefix = "r";
      // } else {
      //   // TODO: if send_recv write 2 files s_ + r_
      //   PDM_error(__FILE__, __LINE__, 0, "You have to choose between CWP_FIELD_EXCH_RECV or CWP_FIELD_EXCH_SEND for field writing type.\n");
      // }

      // std::string fieldName = prefix + "_" + field ->fieldIDGet();
      // std::string fieldComputedName = fieldName + "_is_computed";

      // int id_var = PDM_writer_var_create(_visu_id,
      //                                    PDM_WRITER_ON,
      //                                    PDMfieldComp,
      //                                    PDMfieldType,
      //                                    fieldName.c_str());

      // PDM_writer_status_t time_dependent = PDM_WRITER_ON;
      // if (_topology == CWP_DYNAMIC_MESH_STATIC) {
      //   time_dependent = PDM_WRITER_OFF;
      // }

      // int id_var_computed = -1;
      // if (field->exchangeTypeGet() != CWP_FIELD_EXCH_SEND) {
      //   id_var_computed = PDM_writer_var_create(_visu_id,
      //                                           time_dependent,
      //                                           PDM_WRITER_VAR_SCALAR,
      //                                           PDMfieldType,
      //                                           fieldComputedName.c_str());
      // }

      // field->visuIdSet(id_var);
      // field->visuIdComputedSet(id_var_computed);

      printf("WriterFieldCreate - beg: '%s'\n",field ->fieldIDGet().c_str());;
      fflush(stdout);

      int id_var_send          = -1;
      int id_var_recv          = -1;
      int id_var_recv_computed = -1;

      if (field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND ||
          field->exchangeTypeGet() == CWP_FIELD_EXCH_SENDRECV) {
        std::string prefix = "s";
        std::string fieldName = prefix + "_" + field ->fieldIDGet();

        id_var_send = PDM_writer_var_create(_visu_id,
                                            PDM_WRITER_ON,
                                            PDMfieldComp,
                                            PDMfieldType,
                                            fieldName.c_str());
        field->idVarSendSet(id_var_send);

        printf("WriterFieldCreate - send: '%s' %d\n",fieldName.c_str(), id_var_send);
      fflush(stdout);
      }


      if (field->exchangeTypeGet() == CWP_FIELD_EXCH_RECV ||
          field->exchangeTypeGet() == CWP_FIELD_EXCH_SENDRECV) {
        std::string prefix = "r";
        std::string fieldName = prefix + "_" + field ->fieldIDGet();
        std::string fieldComputedName = fieldName + "_is_computed";

        id_var_recv = PDM_writer_var_create(_visu_id,
                                            PDM_WRITER_ON,
                                            PDMfieldComp,
                                            PDMfieldType,
                                            fieldName.c_str());

        PDM_writer_status_t time_dependent = PDM_WRITER_ON;
        if (_topology == CWP_DYNAMIC_MESH_STATIC) {
          time_dependent = PDM_WRITER_OFF;
        }
        printf("WriterFieldCreate - recv: '%s' %d\n",fieldName.c_str(), id_var_recv);
      fflush(stdout);

        id_var_recv_computed = PDM_writer_var_create(_visu_id,
                                                     time_dependent,
                                                     PDM_WRITER_VAR_SCALAR,
                                                     PDMfieldType,
                                                     fieldComputedName.c_str());

        field->idVarRecvSet(id_var_recv);
        field->idVarRecvComputedSet(id_var_recv_computed);

        printf("WriterFieldCreate - recv: '%s' %d\n",fieldName.c_str(), id_var_recv);
      fflush(stdout);
        printf("WriterFieldCreate - recv 2: '%s' %d\n",fieldComputedName.c_str(), id_var_recv_computed);
      fflush(stdout);


      }

  }


/********************************************************/

  void Visu::fieldDataSet(Field* field, CWP_Field_map_t storage_type, int i_part) {


    int id_var = -1;

    if (storage_type == CWP_FIELD_MAP_SOURCE) {
      id_var = field->idVarSendGet();
    } else {
      id_var = field->idVarRecvGet();
    }

    void* data = field->dataGet(i_part, storage_type);
    //TODO: CHange double for multitype

    PDM_writer_var_set(_visu_id, id_var, _visu_mesh_id, i_part,(double*)data);
  }

  void Visu::fieldDataFree(Field* field) {

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND ||
        field->exchangeTypeGet() == CWP_FIELD_EXCH_SENDRECV) {
      PDM_writer_var_data_free(_visu_id, field->idVarSendGet());
    }

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_RECV ||
        field->exchangeTypeGet() == CWP_FIELD_EXCH_SENDRECV) {
      PDM_writer_var_data_free(_visu_id, field->idVarRecvGet());
      PDM_writer_var_data_free(_visu_id, field->idVarRecvComputedGet());
    }

  }


/********************************************************/

  void Visu::WriterField(Field* field, int* n_ref_values, int **ref_values,  const CWP_Field_map_t  map_type) {



    int id_var          = -1;
    int id_var_computed = -1;

    double default_val = 0;//1e15;

    if (map_type == CWP_FIELD_MAP_SOURCE) {
      id_var          = field->idVarSendGet();

      printf("WriterField - send : '%s' %d\n",field ->fieldIDGet().c_str(), id_var);
      fflush(stdout);

    } else {
      id_var          = field->idVarRecvGet();
      printf("WriterField - recv : '%s' %d\n",field ->fieldIDGet().c_str(), id_var);
      fflush(stdout);
      id_var_computed = field->idVarRecvComputedGet();
      printf("WriterField - recv computed : '%s' %d\n",field ->fieldIDGet().c_str(), id_var_computed);
      fflush(stdout);
    }

    PDM_writer_var_data_free(_visu_id, id_var);
    if (map_type == CWP_FIELD_MAP_TARGET) {
      PDM_writer_var_data_free(_visu_id, id_var_computed);
    }

    double **cp_field_data        = NULL;
    double **cp_field_is_computed = NULL;

    if (n_ref_values == NULL) {

      if (map_type == CWP_FIELD_MAP_TARGET) {
        cp_field_is_computed = (double **) malloc (sizeof (void *) * _n_part);
      }

      for (int i = 0; i < _n_part; i++) {
        PDM_writer_var_set(_visu_id, id_var, _visu_mesh_id, i, (double*) field->dataGet(i, map_type));

        int n_elt_part;

        if (field->locationGet() == CWP_DOF_LOCATION_NODE) {
          n_elt_part = field->meshGet()->getPartNVertex(i);
        }

        else if (field->locationGet() == CWP_DOF_LOCATION_CELL_CENTER) {
          n_elt_part = field->meshGet()->getPartNElts(i);
        }

        else if (field->locationGet() == CWP_DOF_LOCATION_USER) {
          n_elt_part = field->couplingGet()->userTargetNGet(i);
        }

        else {
          PDM_error (__FILE__, __LINE__, 0, "Visu::WriterField : Field location is undefined\n");
        }

        if (map_type == CWP_FIELD_MAP_TARGET) {
          cp_field_is_computed[i] = (double *) malloc (sizeof (double) * n_elt_part);
          for (int j = 0; j < n_elt_part; j++) {
            cp_field_is_computed[i][j] = 1.;
          }

          PDM_writer_var_set(_visu_id, id_var_computed, _visu_mesh_id, i, (double *) cp_field_is_computed[i]);
        }
      }
    }

    else {
      cp_field_data        = (double **) malloc (sizeof (void *) * _n_part);
      cp_field_is_computed = (double **) malloc (sizeof (void *) * _n_part);

      for (int i = 0; i < _n_part; i++) {
        int n_elt_part;

        if (field->locationGet() == CWP_DOF_LOCATION_NODE) {
          n_elt_part = field->meshGet()->getPartNVertex(i);
        }

        else if (field->locationGet() == CWP_DOF_LOCATION_CELL_CENTER) {
          n_elt_part = field->meshGet()->getPartNElts(i);
        }

        else if (field->locationGet() == CWP_DOF_LOCATION_USER) {
          n_elt_part = field->couplingGet()->userTargetNGet(i);
        }

        else {
          PDM_error (__FILE__, __LINE__, 0, "Visu::WriterField : Field location is undefined\n");
        }

        cp_field_data[i] = (double *) malloc (sizeof (double) * field->nComponentGet() * n_elt_part);
        cp_field_is_computed[i] = (double *) malloc (sizeof (double) * n_elt_part);
        for (int j = 0; j < n_elt_part; j++) {
          cp_field_data[i][j] = default_val;
          cp_field_is_computed[i][j] = 0.;
        }

        double* data = (double*) field->dataGet(i, map_type);
        for (int j = 0; j < n_ref_values[i]; j++) {
          for (int k = 0; k < field->nComponentGet(); k++) {
            cp_field_data[i][field->nComponentGet() * (ref_values[i][j]-1) + k] = data[field->nComponentGet() * j + k];
          } 
          cp_field_is_computed[i][ref_values[i][j]-1] = 1.;
        }

        PDM_writer_var_set(_visu_id, id_var,          _visu_mesh_id, i, (double *) cp_field_data[i]);
        PDM_writer_var_set(_visu_id, id_var_computed, _visu_mesh_id, i, (double *) cp_field_is_computed[i]);
      }

    }


      printf("WriterField - PDM_writer_var_write : '%s' %d\n",field ->fieldIDGet().c_str(), id_var);
      fflush(stdout);

    PDM_writer_var_write(_visu_id, id_var);
    if (map_type == CWP_FIELD_MAP_TARGET) {
      printf("WriterField - PDM_writer_var_write var computed: '%s'\n",field ->fieldIDGet().c_str());
      fflush(stdout);
      PDM_writer_var_write(_visu_id, id_var_computed);
    }

    if (n_ref_values != NULL) {
      for (int i = 0; i < _n_part; i++) {
        free (cp_field_data[i]);
      }
      free (cp_field_data);
    }

    if (map_type == CWP_FIELD_MAP_TARGET) {
      for (int i = 0; i < _n_part; i++) {
        free (cp_field_is_computed[i]);
      }
      free (cp_field_is_computed);
    }

      printf("WriterField - fin : '%s'\n",field ->fieldIDGet().c_str());
      fflush(stdout);

  }

/********************************************************/

  void Visu::WriterStepBegin(double physical_time,Mesh* mesh) {
    PDM_writer_step_beg(_visu_id,physical_time);
    _physical_time = physical_time;
    if(_topology != CWP_DYNAMIC_MESH_STATIC){
      for(int i_part= 0 ; i_part < _n_part; i_part++) {
        int nVertex = mesh->getPartNVertex(i_part);
        double* coords = mesh->getVertexCoords(i_part);
        CWP_g_num_t* gnum = mesh->getVertexGNum(i_part);

        GeomCoordSet(i_part,
                        nVertex,
                        coords,
                        gnum);
      }//loop i_part

      // int* blockIDs = mesh->blockDBGet(); // !! block_id from PDM
      int  nBlock   = mesh->nBlockGet();

      for(int i_block=0;i_block<nBlock;i_block++){
        int id_block = i_block;//blockIDs[i_block];
        CWP_Block_t type = mesh->blockTypeGet(id_block);
        int idBlockVisu = GeomBlockAdd(type);

        for(int i_part= 0 ; i_part < _n_part; i_part++) {
          int n_elts = mesh->getBlockNElts(id_block,i_part);
          CWP_g_num_t* gnum = mesh->globalNumGet(id_block,i_part);
          if(type != CWP_BLOCK_FACE_POLY && type != CWP_BLOCK_CELL_POLY) {
            int* connec = mesh->getStdConnectivity(id_block,i_part);
            GeomBlockStdSet(idBlockVisu,
                            i_part,
                            n_elts,
                            connec,
                            gnum);
          }
          else if (type == CWP_BLOCK_FACE_POLY) {
            int* connecIdx = mesh->getPoly2DConnectivityIndex(id_block,i_part);
            int* connec = mesh->getPoly2DConnectivity(id_block,i_part);
            GeomBlockPoly2D(idBlockVisu,
                            i_part,
                            n_elts,
                            connecIdx,
                            connec,
                            gnum);
          }
          // else if (type == CWP_BLOCK_CELL_POLY) {
          //   int  n_faces =;
          //   int* connec_faces_idx =;
          //   int* connec_faces =;
          //   int* connec_cells_idx =; 
          //   int* connec_cells = ;
          //   GeomBlockPoly3D(idBlockVisu,
          //                   i_part,
          //                   n_elts,
          //                   n_faces,
          //                   connec_faces_idx,
          //                   connec_faces,
          //                   connec_cells_idx,
          //                   connec_cells,
          //                   gnum);
          // }
        }//loop on i_part
      } //end loop on block
      GeomWrite(mesh);
    }
  }

/********************************************************/

  void Visu::WriterStepEnd() {

     if(_topology != CWP_DYNAMIC_MESH_STATIC)
       PDM_writer_geom_data_reset (_visu_id, _visu_mesh_id);

     PDM_writer_step_end(_visu_id);
  }


/********************************************************/

  PDM_writer_elt_geom_t Visu::PdmWriterBlockTypeFromCwpBlockType(CWP_Block_t CWP_block_type
                                                                ) {
     PDM_writer_elt_geom_t elt_type;
     switch (CWP_block_type) {
       case CWP_BLOCK_NODE: elt_type = PDM_WRITER_POINT;
       break;

       case CWP_BLOCK_EDGE2: elt_type = PDM_WRITER_BAR2;
       break;

       case CWP_BLOCK_FACE_TRIA3: elt_type = PDM_WRITER_TRIA3;
       break;

       case CWP_BLOCK_FACE_QUAD4: elt_type = PDM_WRITER_QUAD4;
       break;

       case CWP_BLOCK_CELL_TETRA4: elt_type = PDM_WRITER_TETRA4;
       break;

       case CWP_BLOCK_FACE_POLY: elt_type = PDM_WRITER_POLY_2D;
       break;

       case CWP_BLOCK_CELL_HEXA8: elt_type = PDM_WRITER_HEXA8;
       break;

       case CWP_BLOCK_CELL_PYRAM5: elt_type = PDM_WRITER_PYRAMID5;
       break;

       case CWP_BLOCK_CELL_PRISM6: elt_type = PDM_WRITER_PRISM6;
       break;

       case CWP_BLOCK_CELL_POLY: elt_type = PDM_WRITER_POLY_3D;
       break;

       default: elt_type = PDM_WRITER_POINT;
                PDM_error(__FILE__, __LINE__, 0, "This argument does not correspond to a PDM_writer_elt_geom_t.\n");

      }
      return elt_type;
  }




}

/**
 * \endcond
 */
