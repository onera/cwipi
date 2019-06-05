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
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstring>
#include <sstream>

#include <bftc_error.h>
#include <bftc_file.h>

#include <fvmc_parall.h>

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "mesh.hxx"
#include "codeProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "cwp.h"

#include "factory.hpp"
#include "field.hpp"

#include "communication.hxx"
#include "visualization.hxx"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) &&  !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

using namespace std;

namespace cwipi {

  /**
   * \typedef FC
   * 
   * \brief Communication Factory
   *
   *  A communication \ref Factory wich makes \ref Communication 
   *  class objects.
   *  The type of communication objects build depends on the
   *  communication type \ref CWP_Comm_t .
   *
   */
   
  typedef Factory<Communication, CWP_Comm_t> FC;
  
   /**
   * \typedef FG
   * 
   * \brief Geometry Factory
   *
   *  A geometry \ref Factory wich makes \ref Geometry class objects.
   *  The type of \ref Geometry objects build depends on the
   *  geometry algorithm type \ref CWP_Geom_t; .
   *
   */
  
  typedef Factory<Geometry, CWP_Geom_t> FG;

  Coupling::Coupling
  (
   const string               &cplId,
   const CWP_Comm_t           cplType,
   const CodeProperties       &localCodeProperties,
   const CodeProperties       &coupledCodeProperties,
   const CWP_Geom_t           geomAlgo,
   const int                  nPart,
   const CWP_Displacement_t   movingStatus,
   const CWP_Freq_t           recvFreqType,
   CouplingDB                 &cplDB          
   )
  :_cplId(cplId),
   _commType(cplType),
   _communication(*(FC::getInstance().CreateObject(cplType))),
   _localCodeProperties(localCodeProperties),
   _coupledCodeProperties(coupledCodeProperties),
   _recvFreqType (recvFreqType),
   _cplDB(cplDB),
   _fields(*(new map < string, Field<double> * >())),  
   _visu(*new Visu(localCodeProperties.intraCommGet())), 
   _mesh(*new Mesh(localCodeProperties.intraCommGet(),NULL,nPart)),
   _geometry(*new std::map <CWP_Field_value_t, Geometry*>()),
   _iteration(new int)
  {

    //In case where the both codes are on the same MPI process.
    if (coupledCodeProperties.localCodeIs()) {
      if (cplDB.couplingIs(coupledCodeProperties, cplId)) {


        //Communication initialization, MPI communicator creation ... 
        _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB); 
        
        //Get distant coupling object
        Coupling &distCpl = cplDB.couplingGet(coupledCodeProperties, cplId);        
        distCpl._communication.init(_communication);
        Visu* visu_cpl = distCpl.visuGet();
        Mesh* mesh_cpl = distCpl.meshGet(); 

        _mesh.setVisu(&_visu); 
        mesh_cpl->setVisu(visu_cpl);  


        std::map <CWP_Field_value_t, Geometry*>* _geometry_cpl = distCpl.geometryGet();

        //Geometry initialization
        _geometry[CWP_FIELD_VALUE_CELL_MEAN] = FG::getInstance().CreateObject(geomAlgo);
        _geometry[CWP_FIELD_VALUE_CELL_POINT] = FG::getInstance().CreateObject(geomAlgo);
        _geometry[CWP_FIELD_VALUE_NODE] = FG::getInstance().CreateObject(geomAlgo);
        _geometry[CWP_FIELD_VALUE_USER] = FG::getInstance().CreateObject(geomAlgo);
        
        (*_geometry_cpl)[CWP_FIELD_VALUE_CELL_MEAN]  = FG::getInstance().CreateObject(geomAlgo);
        (*_geometry_cpl)[CWP_FIELD_VALUE_CELL_POINT] = FG::getInstance().CreateObject(geomAlgo);
        (*_geometry_cpl)[CWP_FIELD_VALUE_NODE]       = FG::getInstance().CreateObject(geomAlgo);
        (*_geometry_cpl)[CWP_FIELD_VALUE_USER]       = FG::getInstance().CreateObject(geomAlgo);

        //Geometry initialization 
        std::map <CWP_Field_value_t, Geometry*>::iterator it = _geometry.begin();
        while (it != _geometry.end()) {    
         (it -> second) -> init(this,it->first);
         it++;
        }           
        
        it = _geometry_cpl->begin();
        while (it != _geometry_cpl->end()) {    
         (it -> second) -> init(&distCpl,it->first);
          it++;
        }

 
      }   
    } // if (coupledCodeProperties.localCodeIs())     
    else {
      //Communication initialization, MPI communicator creation ... 
      _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);

      _mesh.setVisu(&_visu);      
      //Geometry creation
      _geometry[CWP_FIELD_VALUE_CELL_MEAN] = FG::getInstance().CreateObject(geomAlgo);
      _geometry[CWP_FIELD_VALUE_CELL_POINT] = FG::getInstance().CreateObject(geomAlgo);
      _geometry[CWP_FIELD_VALUE_NODE] = FG::getInstance().CreateObject(geomAlgo);
      _geometry[CWP_FIELD_VALUE_USER] = FG::getInstance().CreateObject(geomAlgo);
    
      //Geometry initialization
      std::map <CWP_Field_value_t, Geometry*>::iterator it = _geometry.begin();
      while (it != _geometry.end()) {    
        (it -> second) -> init(this,it->first);
        it++;
      }
    
           
    } // end else
        
    
  }


  Coupling::~Coupling()
  {
  
    if(_visu.isCreated()) {
       _visu.WriterStepEnd();
    }
       
    #if defined(DEBUG) && 0
    cout << "destroying '" << _name << "' coupling : TODO" << endl;
    #endif
  }

 //TODO: Virer ptFortranInterpolationFct
  void 
  Coupling::issend
  (string &sendingFieldID) {

     std:map <std::string, Field<double> *>::iterator it;
     it = _fields.find(sendingFieldID);

     if (it != _fields.end()) {
       Field <double>* sendingField = it -> second;   
      _geometry[sendingField -> typeGet()] -> issend(sendingField);
     }
  }

  void 
  Coupling::irecv
  (string &recevingFieldID) {

     std:map <std::string, Field<double> *>::iterator it;
     it = _fields.find(recevingFieldID);

     if (it != _fields.end()) {
       Field <double>* recevingField = it -> second;   
       _geometry[recevingField -> typeGet()] ->irecv(recevingField);
     }
  }


  int
  Coupling::fieldNComponentGet
  (
   const string &field_id
  )
  {
  
    map<string,Field<double>*>::iterator It = _fields.find(field_id.c_str());  
    if (It == _fields.end()) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' not existing field\n", field_id.c_str());
    }

    return It->second->nComponentGet();
  }

  bool 
  Coupling::fieldIs
  (
   const string &field_id
  )
  {
    map<string,Field<double>*>::iterator It = _fields.find(field_id.c_str());
    return (It != _fields.end());
  }

  void
  Coupling::fieldCreate
  (
   const string               &field_id,
   const CWP_Type_t           data_type,
   const CWP_Field_storage_t  storage,
   const int                  n_component,
   const CWP_Field_value_t    fieldType,
   const CWP_Field_exch_t     exch_type,
   const CWP_Status_t         visu_status
  )
  {

    if (fieldIs(field_id)) {
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' existing field\n", field_id.c_str());
    }
    
    //
    // Create the new field
    
    
    double physTime=0.0;
    *_iteration = 0;
    
    cwipi::Field<double> *newField = new cwipi::Field<double>(field_id,
                                                              &_mesh,
                                                               fieldType,
                                                               storage,
                                                               n_component,
                                                               exch_type,
                                                               visu_status,
                                                               _iteration, //iteration
                                                               &physTime);  //physTime

    pair<string, Field<double>* > newPair(string(field_id), newField);

    _fields.insert(newPair);
    
    if (_visu.isCreated())
      _visu.WriterFieldCreate(newField);
  }
  

   void 
   Coupling::geomCompute (CWP_Field_value_t geometryLocation)
   {  
     _geometry[geometryLocation] -> compute();
     
   }


  int Coupling::nUncomputedTargetsGet
  (
    const CWP_Field_value_t geometryLocation,
    const int  i_part
  )
  {
    _geometry[geometryLocation] -> nUncomputedTargetsGet(i_part);
  
  }


   void 
   Coupling::waitIssend
   (
    string &sendingFieldID
   )
   {
     std:map <std::string, Field<double> *>::iterator it;

     it = _fields.find(sendingFieldID);

     if (it != _fields.end()) {
       Field <double>* sendingField = it -> second;   
      _geometry[sendingField -> typeGet()] -> waitIssend(sendingField);
     }
   }


   void 
   Coupling::waitIrecv
   (
    string &recevingFieldID
   )
   {
     std:map <std::string, Field<double> *>::iterator it;
     it = _fields.find(recevingFieldID);

     if (it != _fields.end()) {
       Field <double>* recevingField = it -> second;   
       _geometry[recevingField -> typeGet()] -> waitIrecv(recevingField);
     }
   }


   
  /**
   *
   * \brief Get field storage type
   * 
   * \param [in]   field_id       Field identifier
   * 
   */

   CWP_Field_storage_t
   Coupling::fieldStorageGet
   (
     const string &field_id
   )
   {
    map<string,Field<double>*>::iterator It = _fields.find(field_id.c_str());  
    if (It == _fields.end()) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' not existing field\n", field_id.c_str());
    }

    return It->second->storageTypeGet();
    
   }


  /**
   *
   * \brief Get field fieldType
   * 
   * \param [in]   field_id       Field identifier
   * 
   */

  CWP_Field_value_t
  Coupling::fieldTypeGet
  (
    const string &field_id
  )
  {
    map<string,Field<double>*>::iterator It = _fields.find(field_id.c_str());  
    if (It == _fields.end()) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' not existing field\n", field_id.c_str());
    }
    return It->second->typeGet();
    
  }


 /**
  *
  * \brief Set data mapping
  * 
  * \param [in]  field_id       Field identifier
  * \param [in]  data           Storage array (Mapping)
  * 
  */
  //TODO:Change to dataSet
  void
  Coupling::fieldDataSet
  (
    const string &field_id,
    int i_part,
    double data[]
  )
  { 
    map<string,Field<double>*>::iterator It = _fields.find(field_id.c_str());  
    if (It == _fields.end())
      {
         bftc_error(__FILE__, __LINE__, 0,
               "'%s' not existing field\n", field_id.c_str());
      }
    else 
      {
        It->second->dataSet(i_part,data);
        if(_visu.isCreated()) {
          printf("_visu.fieldDataSet(It->second,i_part); \n");
          _visu.fieldDataSet(It->second,i_part);
        }      
      }   
  }



  /**
   *
   * \brief Removing a field
   * 
   * \param [in]   field_id       Field identifier
   * 
   */
  void
  Coupling::fieldDel
  (
    const string &field_id
  )
  {
    map<string,Field<double>*>::iterator It = _fields.find(field_id.c_str());  
    if (It == _fields.end())
      {
         bftc_error(__FILE__, __LINE__, 0,
               "'%s' not existing field\n", field_id.c_str());
      }
    else 
      {
        delete It->second;
      }       
   
  }

  
  
  void Coupling::meshVtcsSet
    (
     const int          i_part,
     const int          n_pts,
     double             coords[],
     CWP_g_num_t        global_num[]
    )
    {
      _mesh.nodal_coord_set(i_part,
                            n_pts,
                            coords,
                            global_num); 
    }
    
  int Coupling::meshBlockAdd
    (const CWP_Block_t     block_type){
      _mesh.blockAdd(block_type);
    }
    
    
  void Coupling::meshStdBlockSet
    (
     const int           i_part,
     const int           block_id,
     const int           n_elts,
     int                 connec[],
     CWP_g_num_t        global_num[]
    )
  {
       _mesh.stdBlockSet( i_part,
                          block_id,
                          n_elts,
                          connec, 
                          global_num
                        );
  }
  
  void Coupling::meshHighOrderBlockSet
    (
     const int           i_part,
     const int           block_id,
     const int           n_elts,
     const int           order,
     int                 connec[],
     CWP_g_num_t         global_num[])
    {
    
    }
  
  void Coupling::meshFPolyBlockSet
    (
     const int            i_part,
     const int            block_id,
     const int            n_elts,
     int                  connec_idx[],
     int                  connec[],
     CWP_g_num_t          global_num[]
    )
    {
     _mesh.poly2DBlockSet(i_part,
                          block_id,
                          n_elts,
                          connec_idx,
                          connec, 
                          global_num
                         );
   }
  

  void Coupling::meshCPolyBlockSet
    (
     const int           i_part,
     const int           block_id,
     const int           n_elts,
     const int           n_faces,
     int                 connec_faces_idx[],
     int                 connec_faces[],
     int                 connec_cells_idx[],
     int                 connec_cells[],
     CWP_g_num_t         global_num[]
    )
    {
       _mesh.poly3DBlockSet(i_part,
                            block_id,
                            n_elts,
                            n_faces,
                            connec_faces_idx,
                            connec_faces    ,
                            connec_cells_idx,
                            connec_cells    , 
                            global_num      
                        );               

   }
  
  
  void Coupling::fvmcNodalShared(const int           i_part,
                      fvmc_nodal_t        *fvmc_nodal)
  {
    
  }
  
  void Coupling::meshFromCellFaceSet(const int   i_part,
                        const int   n_cells,
                        int         cell_face_idx[],
                        int         cell_face[],
                        int         n_faces,
                        int         face_vtx_idx[],
                        int         face_vtx[],
                        CWP_g_num_t parent_num[])
  {
     _mesh.fromCellFaceSet(i_part,
                               n_cells,
                               cell_face_idx,
                               cell_face,
                               n_faces,
                               face_vtx_idx,
                               face_vtx,
                               parent_num);
  } 
  
  
  
  
  
  void Coupling::meshFromFacesEdgeSet(const int   i_part,
                         const int   n_faces,
                         int         face_edge_idx[],
                         int         face_edge[],
                         const int   n_edges,
                         int         edge_vtx_idx[],
                         int         edge_vtx[],
                         CWP_g_num_t parent_num[])
  {
     _mesh.fromFacesEdgeSet(i_part,
                            n_faces,
                            face_edge_idx,
                            face_edge,
                            n_edges,
                            edge_vtx_idx,
                            edge_vtx,
                            parent_num); 
  }


  void Coupling::meshDel()
  {_mesh.meshDel();
  }


  void Coupling::visuSet(const int               freq,
                         const CWP_Visu_format_t format,
                         const char             *format_option
                         ) {
      int max_codename_length = 150;
      printf("_localCodeProperties.nameGet();\n");
      string name = _localCodeProperties.nameGet();
      int local_rank = -1;
      printf("_localCodeProperties.intraCommGet(),\n");
     // MPI_Comm_rank(_localCodeProperties.intraCommGet(),&local_rank);
      

      char output_name [name.length()];
        printf("sprintf(output_name,,name.c_str());\n");
      sprintf(output_name,"%s",name.c_str());

      printf("_visu.VisuCreate(fI%sI\n",output_name);
      _visu.VisuCreate(freq,
                       format,
                       format_option,
                       "cwipi",
                       output_name);   
        printf("_visu.GeomCreate(_mesh.getNPar) %i\n",_mesh.getNPart());  
                                   
      _visu.GeomCreate(_mesh.getNPart());
          
  }


} // namespace cwipi

