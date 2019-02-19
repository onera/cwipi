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
#include "geometry.hxx"
#include "communication.hxx"


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
   _mesh(*new Mesh(localCodeProperties.intraCommGet(),nPart))
   //_geometry(*(FG::getInstance().CreateObject(geomAlgo)))
  {
  
    if (coupledCodeProperties.localCodeIs()) {
      if (cplDB.couplingIs(coupledCodeProperties, cplId)) {
        _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);
        Coupling &distCpl = cplDB.couplingGet(coupledCodeProperties, cplId);
        distCpl._communication.init(_communication);
      }
    }
  
    else {
      _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);
    }
    
    
    
  }


  Coupling::~Coupling()
  {
    #if defined(DEBUG) && 0
    cout << "destroying '" << _name << "' coupling : TODO" << endl;
    #endif
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
   const CWP_Field_value_t    nature,
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

    cwipi::Field<double> *newField = new cwipi::Field<double>(data_type,
                                           storage,
                                           n_component,
                                           nature,
                                           exch_type,
                                           visu_status);  

    pair<string, Field<double>* > newPair(string(field_id), newField);

    _fields.insert(newPair);
    bftc_printf("champ '%s' a été créé\n", field_id.c_str());
  }


  /**
   *
   * \brief Get field data type
   * 
   * \param [in]   field_id       Field identifier
   *
   * \return                      Field data type
   * 
   */
   CWP_Type_t
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
    
    return It->second->dataTypeGet();
    
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
   * \brief Get field nature
   * 
   * \param [in]   field_id       Field identifier
   * 
   */

  CWP_Field_value_t
  Coupling::fieldNatureGet
  (
    const string &field_id
  )
  {
    map<string,Field<double>*>::iterator It = _fields.find(field_id.c_str());  
    if (It == _fields.end()) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' not existing field\n", field_id.c_str());
    }
    return It->second->natureGet();
    
  }


 /**
  *
  * \brief Set data mapping
  * 
  * \param [in]  field_id       Field identifier
  * \param [in]  data           Storage array (Mapping)
  * 
  */
  void
  Coupling::fieldMappingSet
  (
    const string &field_id,
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
        It->second->mappingSet(data);
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
    
  void Coupling::meshEndSet()
  {_mesh.endSet();
  } 
    
    
  void Coupling::meshBlockAdd
  (
    const int           i_part,
    const CWP_Block_t   block_type,
    const int           n_elts,
    int                 connec[],
    CWP_g_num_t         global_num[],
    CWP_g_num_t         parent_num[]
  )
  {
     _mesh.blockAdd(i_part,
                    block_type,
                    n_elts,
                    NULL,
                    connec,
                    -1,
                    NULL,
                    NULL,                            
                    NULL,
                    (int*) parent_num);
  }
  
  void Coupling::meshHighOrderBlockAdd
    (
     const int           i_part,
     const CWP_Block_t   block_type,
     const int           n_elts,
     const int           order,
     int                 connec[],
     CWP_g_num_t         global_num[]
    )
    {
    
    }
  
  
  void Coupling::meshFPolyBlockAdd
  (
    const int            i_part,
    const int            n_elts,
    int                  connec_idx[],
    int                  connec[],
    CWP_g_num_t          parent_num[]
  )
  {
     _mesh.blockAdd(i_part,
                    CWP_BLOCK_FACE_POLY,
                    n_elts,
                    connec_idx,
                    connec,
                    -1,
                    NULL,
                    NULL,                            
                    NULL,
                    (int*)parent_num);
  }
  

  void Coupling::meshCPolyBlockAdd
  (
    const int           i_part,
    const int           n_elts,
    int                 cell_face_idx[],
    int                 cell_face[],
    const int           n_faces,
    int                 face_vtx_idx[],
    int                 face_vtx[],
    CWP_g_num_t         parent_num[]
  )
  {
     _mesh.blockAdd(i_part,
                    CWP_BLOCK_CELL_POLY,
                    n_elts,
                    cell_face_idx,
                    cell_face,
                    n_faces,
                    face_vtx_idx,
                    face_vtx,                            
                    NULL,
                    (int*)parent_num);
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


} // namespace cwipi

