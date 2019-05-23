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
#include "bftc_error.h"
#include <unistd.h> 
#include <stdio.h>
#include "field.hpp"

namespace cwipi {



  Visu::Visu(const MPI_Comm &MPIComm):_visu_id(-1),_visu_mesh_id(-1),_freq(-1),_physical_time(-1),
                                      _visuCreated(false), 
                                      _output_dir(NULL), 
                                      _output_name(NULL),
                                      _divide_polygons(PDM_WRITER_OFF),
                                      _divide_polyhedra(PDM_WRITER_OFF) {
                                      
                                      
                                      
     _pdmComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&MPIComm)); 

  }
  
  Visu::~Visu() {
    PDM_writer_free(_visu_id);
  }


  void Visu::VisuCreate(const int          freq,
                        CWP_Visu_format_t  format,
                        const char        *format_option,
                        char        *output_dir,
                        char        *output_name) {
   
    PDM_writer_fmt_fic_t fmt_fic      = PDM_WRITER_FMT_BIN;
    const char* fmt                   = "Ensight";
    PDM_writer_topologie_t topologie  = PDM_WRITER_TOPO_CONSTANTE;
    PDM_writer_statut_t st_reprise    = PDM_WRITER_OFF;
    const char *options_comp          = "";
    //Proportion of working node for file acess 
    int working_node = 1;
    PDM_io_acces_t acess_type =   PDM_IO_ACCES_MPI_SIMPLE;

    _output_dir  = output_dir; 
    _output_name = output_name;
    
    std::string str_options = format_option;
    std::string delimiter = ",";
    
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
      else if (option != "")                 bftc_error(__FILE__, __LINE__, 0, 
                                                        "Not a valid visualization option.\n");

    } 
    while (pos != std::string::npos);

    _visu_id = PDM_writer_create(fmt,
                                 fmt_fic,   
                                 topologie,
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
                                             _divide_polygons,
                                             _divide_polyhedra,
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
                              global_num);                   
   }

/*****************************************/

  int Visu::GeomBlockAdd(CWP_Block_t blockType) {
  
     int id_block = PDM_writer_geom_bloc_add(_visu_id,
                               _visu_mesh_id,
                               PDM_WRITER_ON,  
                               PdmWriterBlockTypeFromCwpBlockType(blockType)
                              ); 
    return id_block;
  }

/*****************************************/

  void Visu::GeomWrite() {
    PDM_writer_geom_write(_visu_id,_visu_mesh_id);
  
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
  
  void Visu::GeomBlockGNumMeshSet (int id_block,
                                   int id_part,
                                   CWP_g_num_t *global_num) {
                                  
      PDM_writer_geom_bloc_g_num_mesh_set(_visu_id,
                                          _visu_mesh_id,  
                                          id_block,     
                                          id_part,  
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


  void Visu::WriterFieldCreate(Field<double>* field) {

      CWP_Field_value_t CWPfielType = field -> typeGet();
      PDM_writer_var_loc_t PDMfieldType;
      
      if(CWPfielType == CWP_FIELD_VALUE_CELL_POINT) PDMfieldType = PDM_WRITER_VAR_ELEMENTS;
      
      if(CWPfielType == CWP_FIELD_VALUE_NODE) PDMfieldType = PDM_WRITER_VAR_SOMMETS;      

      int id_var = PDM_writer_var_create(_visu_id, 
                                     PDM_WRITER_ON,
                                     PDM_WRITER_VAR_SCALAIRE, 
                                     PDMfieldType, 
                                     (field ->fieldIDGet()).c_str());
      field -> visuIdSet(id_var);
  }


/********************************************************/

  void Visu::fieldDataSet(Field<double>* field,int i_part) {
    
    int id_var = -1;
 
    id_var = field -> visuIdGet();
    double* data = field -> dataGet(i_part);
   
    for(int i=0;i<20;i++) {
   //   printf("data [%i] %f\n",i,(double) data[i]);
    }

    printf("PDM_writer_var_set\n");

    PDM_writer_var_set(_visu_id, id_var, _visu_mesh_id, i_part, data);
   // if(i_part==1)  {  while(1==1){} }
  }

/********************************************************/

  void Visu::WriterField(Field<double>* field) {
    int id_var = -1;

    id_var = field -> visuIdGet();
    
       
    for (int i_part =0;i_part<_n_part;i_part++)
       {   fieldDataSet(field,i_part);
       printf("OOOOO %i\n",i_part);
        
       }
       
  
    PDM_writer_var_write(_visu_id, id_var);                       
   
  }
  
/********************************************************/

  void Visu::WriterStepBegin(double physical_time) {
     PDM_writer_step_beg(_visu_id,physical_time);
     _physical_time = physical_time;  
  }

/********************************************************/

  void Visu::WriterStepEnd() {
     PDM_writer_step_end(_visu_id);  
  }


/********************************************************/
  
  PDM_writer_elt_geom_t Visu::PdmWriterBlockTypeFromCwpBlockType(CWP_Block_t CWP_block_type
                                                                ) {
                                                                   
     switch (CWP_block_type) {

       case CWP_BLOCK_NODE: return PDM_WRITER_POINT;
       break;
       
       case CWP_BLOCK_EDGE2: return PDM_WRITER_BAR2;
       break;
   
       case CWP_BLOCK_FACE_TRIA3: return PDM_WRITER_TRIA3;
       break;

       case CWP_BLOCK_FACE_QUAD4: return PDM_WRITER_QUAD4;
       break;
                       
       case CWP_BLOCK_CELL_TETRA4: return PDM_WRITER_TETRA4;
       break;

       case CWP_BLOCK_FACE_POLY: return PDM_WRITER_POLY_2D;
       break;
       
       case CWP_BLOCK_CELL_HEXA8: return PDM_WRITER_HEXA8;
       break;

       case CWP_BLOCK_CELL_PYRAM5: return PDM_WRITER_PYRAMID5;
       break;
       
       case CWP_BLOCK_CELL_PRISM6: return PDM_WRITER_PRISM6;
       break;
       
       case CWP_BLOCK_CELL_POLY: return PDM_WRITER_POLY_3D;
       break;
       
      }
      
  }




}
