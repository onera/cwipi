/*   This file is part of the CWIPI library.

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
#include "field.hxx"

#include "communication.hxx"
#include "visualization.hxx"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

/**
 * \cond
 */

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
   * \brief Spatial Interpolation Factory
   *
   *  A Spatial Interpolation \ref Factory wich makes  \ref SpatialInterp class objects.
   *  The type of \ref SpatialInterp objects build depends on the
   *  spatial interpolation algorithm type \ref CWP_Spatial_interp_t; .
   *
   */

  typedef Factory<SpatialInterp, CWP_Spatial_interp_t> FG;



  Coupling::Coupling
  (
   const string               &cplId,
   const CWP_Comm_t           cplType,
   const CodeProperties       &localCodeProperties,
   const CodeProperties       &coupledCodeProperties,
   const CWP_Interface_t      entities_dim,
   const CWP_Spatial_interp_t           spatialInterpAlgo,
   const int                  nPart,
   const CWP_Dynamic_mesh_t   displacement,
   const CWP_Time_exch_t           recvFreqType,
   CouplingDB                 &cplDB
   )
  :_cplId(cplId),
   _commType(cplType),
   _communication(*(FC::getInstance().CreateObject(cplType))),
   _localCodeProperties(localCodeProperties),
   _coupledCodeProperties(coupledCodeProperties),
   _entities_dim(entities_dim),
   _spatial_interp(*new std::map <CWP_Dof_location_t, SpatialInterp*>()),
   _mesh(*new Mesh(localCodeProperties.connectableCommGet(),NULL,nPart,displacement,this)),
   _recvFreqType (recvFreqType),
   _visu(*new Visu(localCodeProperties.connectableCommGet(),displacement)),
   _fields(*(new map < string, Field * >())),
   _cplDB(cplDB),
   _iteration(new int),
   _displacement(displacement)
  {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

     //In case where the both codes are on the same MPI process.
    if (coupledCodeProperties.localCodeIs()) {
      if (cplDB.couplingIs(coupledCodeProperties, cplId) ) {
        //Communication initialization, MPI communicator creation ...
        _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);

        //Get distant coupling object
        Coupling &distCpl = cplDB.couplingGet(coupledCodeProperties, cplId);
        distCpl._communication.init(_communication);

        Visu* visu_cpl = distCpl.visuGet();
        Mesh* mesh_cpl = distCpl.meshGet();

        _mesh.setVisu(&_visu);
        mesh_cpl->setVisu(visu_cpl);

        std::map <CWP_Dof_location_t, SpatialInterp*>* _spatial_interp_cpl = distCpl.spatialInterpGet();

        //SpatialInterp initialization
        //_spatial_interp[CWP_FIELD_VALUE_CELL_MEAN] = FG::getInstance().CreateObject(spatialInterpAlgo);
        _spatial_interp[CWP_DOF_LOCATION_CELL_CENTER] = FG::getInstance().CreateObject(spatialInterpAlgo);
        _spatial_interp[CWP_DOF_LOCATION_NODE] = FG::getInstance().CreateObject(spatialInterpAlgo);
        _spatial_interp[CWP_DOF_LOCATION_USER] = FG::getInstance().CreateObject(spatialInterpAlgo);

        //(*_spatial_interp_cpl)[CWP_FIELD_VALUE_CELL_MEAN]  = FG::getInstance().CreateObject(spatialInterpAlgo);
        (*_spatial_interp_cpl)[CWP_DOF_LOCATION_CELL_CENTER] = FG::getInstance().CreateObject(spatialInterpAlgo);
        (*_spatial_interp_cpl)[CWP_DOF_LOCATION_NODE]       = FG::getInstance().CreateObject(spatialInterpAlgo);
        (*_spatial_interp_cpl)[CWP_DOF_LOCATION_USER]       = FG::getInstance().CreateObject(spatialInterpAlgo);

        //SpatialInterp initialization
        std::map <CWP_Dof_location_t, SpatialInterp*>::iterator it = _spatial_interp_cpl->begin();

        while (it != _spatial_interp_cpl->end()) {
         (it -> second) -> init(&distCpl,it->first,1);
          it++;
        }

        it = _spatial_interp.begin();
        while (it != _spatial_interp.end()) {
         (it -> second) -> init(this,it->first,0);
         it++;
        }

      }
    } // if (coupledCodeProperties.localCodeIs())
    else {
      //Communication initialization, MPI communicator creation ...
      _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);

      //Tips to provide the correct visu Communicator in case CWP_COMM_PAR_WITHOUT_PART
      int unionRank;
      MPI_Comm_rank(_communication.unionCommGet(),&unionRank);

      MPI_Comm visuComm;
      if(commTypeGet() == CWP_COMM_PAR_WITHOUT_PART) {
        MPI_Group unionGroup;
        MPI_Group intraGroup;
        MPI_Comm_group(_localCodeProperties.connectableCommGet(), &intraGroup);
        MPI_Comm_group(_communication.unionCommGet(), &unionGroup);
        MPI_Group visuGroup;
        int locRootRank = _communication.unionCommLocCodeRootRanksGet();
        int locRootRankIntra;
        MPI_Group_translate_ranks(unionGroup, 1, &locRootRank,
                                  intraGroup , &locRootRankIntra);

        MPI_Group_incl(intraGroup, 1, &locRootRankIntra, &visuGroup);

        MPI_Comm_create(_localCodeProperties.connectableCommGet(), visuGroup, &visuComm);

        if(unionRank == _communication.unionCommLocCodeRootRanksGet()){
           _visu = *new Visu(visuComm,displacement);
        }
      }

       _mesh.setVisu(&_visu);
      //SpatialInterp creation
     // _spatial_interp[CWP_FIELD_VALUE_CELL_MEAN] = FG::getInstance().CreateObject(spatialInterpAlgo);
      _spatial_interp[CWP_DOF_LOCATION_CELL_CENTER] = FG::getInstance().CreateObject(spatialInterpAlgo);
      _spatial_interp[CWP_DOF_LOCATION_NODE] = FG::getInstance().CreateObject(spatialInterpAlgo);
      _spatial_interp[CWP_DOF_LOCATION_USER] = FG::getInstance().CreateObject(spatialInterpAlgo);
      //SpatialInterp initialization
        std::map <CWP_Dof_location_t, SpatialInterp*>::iterator it = _spatial_interp.begin();
        while (it != _spatial_interp.end()) {
          (it -> second) -> init(this,it->first,0);
          it++;
        }

    } // end else

  }


  Coupling::~Coupling()
  {

    if(_visu.isCreated()) {
       _visu.WriterStepEnd();
    }

    std::map <CWP_Dof_location_t, SpatialInterp*>::iterator it = _spatial_interp.begin();
    while (it != _spatial_interp.end()) {
       // delete it -> second;
        it++;
    }


    std::map < string, Field * >::iterator itf = _fields.begin();
    while (itf != _fields.end()) {
      if(_visu.isCreated() && itf -> second -> visuStatusGet() == CWP_STATUS_ON )
        _visu.fieldDataFree(itf -> second);
      delete itf -> second;
      itf++;
    }

    if(_visu.isCreated()) {
      // _visu.SpatialInterpFree();

      delete _iteration;
    }


    #if defined(DEBUG) && 0
    cout << "destroying '" << _name << "' coupling : TODO" << endl;
    #endif
  }


 //TODO: Virer ptFortranInterpolationFct
  void
  Coupling::issend
  (
    string &sendingFieldID
  )
  {
    map <string, Field *>::iterator it;
    it = _fields.find(sendingFieldID);

    if (it != _fields.end()) {
      Field* sendingField = it -> second;
      if(_spatial_interp[sendingField -> associatedCloudPointTypeGet()] -> _both_codes_are_local == 0){
        _spatial_interp[sendingField -> associatedCloudPointTypeGet()] -> issend_p2p(sendingField);
        return;
      }
      else {
        Coupling &distCpl = _cplDB.couplingGet(_coupledCodeProperties, _cplId);
        map <std::string, Field *>::iterator it_recv = distCpl._fields.find(sendingFieldID);
        if (it_recv != distCpl._fields.end()) {
          Field* recevingField = it_recv -> second;
          _spatial_interp[sendingField -> associatedCloudPointTypeGet()] -> both_codes_on_the_same_process_exchange_p2p(sendingField,recevingField);
        }
      }
    }
  }

  void
  Coupling::irecv
  (
    string &recevingFieldID
  ) 
  {
    map <string, Field *>::iterator it = _fields.find(recevingFieldID);
    if (it != _fields.end()) {
      Field* recevingField = it -> second;
      if(_spatial_interp[recevingField -> associatedCloudPointTypeGet()] -> _both_codes_are_local == 0 ){
        _spatial_interp[recevingField -> associatedCloudPointTypeGet()] -> irecv_p2p(recevingField);
      } 
      return;
    }
  }


  int
  Coupling::fieldNComponentGet
  (
   const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
       PDM_error(__FILE__, __LINE__, 0, "'%s' not existing field\n", field_id.c_str());
    }
    return It->second->nComponentGet();
  }

  bool
  Coupling::fieldIs
  (
   const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    return (It != _fields.end());
  }

  void
  Coupling::fieldCreate
  (
   const string               &field_id,
   const CWP_Type_t           data_type,
   const CWP_Field_storage_t  storage,
   const int                  n_component,
   const CWP_Dof_location_t    fieldType,
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
    cwipi::Field *newField = new cwipi::Field(field_id,
                                              data_type,
                                              this,
                                              fieldType,
                                              storage,
                                              n_component,
                                              exch_type,
                                              visu_status,
                                              _iteration, //iteration
                                             &physTime);  //physTime


    pair<string, Field* > newPair(string(field_id), newField);
    string localName = _localCodeProperties.nameGet();
    _fields.insert(newPair);
    if (_visu.isCreated() && newField -> visuStatusGet() == CWP_STATUS_ON) {
      _visu.WriterFieldCreate(newField);
    }
  }


   void
   Coupling::spatialInterpWeightsCompute (CWP_Dof_location_t pointsCloudLocation, CWP_Field_exch_t exchange_type)
   {
     _spatial_interp[pointsCloudLocation] -> spatialInterpWeightsCompute(exchange_type);

   }


  int Coupling::nUncomputedTargetsGet
  (
    const CWP_Dof_location_t pointsCloudLocation,
    const int  i_part
  )
  {
    return _spatial_interp[pointsCloudLocation] -> nUncomputedTargetsGet(i_part);
  }


   void
   Coupling::waitIssend
   (
    string &sendingFieldID
   )
   {

     map <string, Field *>::iterator it;
     it = _fields.find(sendingFieldID);

     if (it != _fields.end()) {
       Field* sendingField = it -> second;
       if(_spatial_interp[sendingField -> associatedCloudPointTypeGet()] -> _both_codes_are_local == 0){
        _spatial_interp[sendingField -> associatedCloudPointTypeGet()] -> waitIssend_p2p(sendingField);
        return;
       }
       else {
        _spatial_interp[sendingField -> associatedCloudPointTypeGet()] -> waitIssend_p2p(sendingField);
        Coupling &distCpl = _cplDB.couplingGet(_coupledCodeProperties, _cplId);

        map <std::string, Field *>::iterator it_recv = distCpl.fieldsGet() -> find(sendingFieldID);
        if (it_recv != distCpl.fieldsGet() -> end() ) {
          Field* recevingField = it_recv -> second;
          distCpl._spatial_interp[recevingField -> associatedCloudPointTypeGet()] -> waitIrecv_p2p(it_recv -> second);
          return;
        }

       }
     }
   }


   void
   Coupling::waitIrecv
   (
    string &recevingFieldID
   )
   {
     map <string, Field *>::iterator it;
     it = _fields.find(recevingFieldID);

     if (it != _fields.end()) {
       Field* recevingField = it -> second;
       if(_spatial_interp[recevingField -> associatedCloudPointTypeGet()] -> _both_codes_are_local == 0)
         _spatial_interp[recevingField -> associatedCloudPointTypeGet()] -> waitIrecv_p2p(recevingField);
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
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
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

  CWP_Dof_location_t Coupling::fieldTypeGet
  (
    const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' not existing field\n", field_id.c_str());
    }
    return It->second->typeGet();

  }


 /**
  *
  * \brief Set Field data
  *
  * \param [in]  field_id       Field identifier
  * \param [in]  data           Storage array (mapping)
  *
  */
  //TODO:Change to dataSet
  void
  Coupling::fieldDataSet
  (
    const string &field_id,
    int i_part,
    void* data
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end())
      {
         bftc_error(__FILE__, __LINE__, 0,
               "'%s' not existing field\n", field_id.c_str());
      }
    else
      {
        It->second->dataSet(i_part,data);
        if(_visu.isCreated() && It -> second -> visuStatusGet() == CWP_STATUS_ON) {
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
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
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
     return _mesh.blockAdd(block_type);
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
/*
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
  */
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


 /* void Coupling::fvmcNodalShared(const int           i_part,
                      fvmc_nodal_t        *fvmc_nodal)
  {

  }


*/


  void Coupling::meshFinalize() 
  {
    _mesh.geomFinalize();
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


  void Coupling::interpFromLocSet (  const string field_id,
                                     CWP_Interp_from_location_t fct
                                  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end())
      {
         bftc_error(__FILE__, __LINE__, 0,
               "'%s' not existing field\n", field_id.c_str());
      }
    else
      {
        It -> second -> interpFromLocationSet(fct);
      }
  }


  void Coupling::userTgtPtsSet (const int i_part,
                                const int n_pts,
                                double    coord[] )
  {
    _spatial_interp[CWP_DOF_LOCATION_USER] -> user_target_points_set(i_part, n_pts, coord);
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
  {
    _mesh.meshDel();
  }


  void Coupling::visuSet(const int               freq,
                         const CWP_Visu_format_t format,
                         const char             *format_option
                         )
  {
    string CodeName = _localCodeProperties.nameGet();
    string cplCodeName = _coupledCodeProperties.nameGet();
    string cplId = IdGet();

    string visuDir = "cwipi";
    char output_name [CodeName.length()];
    char output_dir   [visuDir.length() + 1 + cplId.length() + 1 + CodeName.length() + 1 + cplCodeName.length()];
    sprintf(output_name,"%s",CodeName.c_str());
    sprintf(output_dir,"%s/%s_%s_%s",visuDir.c_str(),cplId.c_str(),CodeName.c_str(),cplCodeName.c_str());
    int rank;
    MPI_Comm_rank(_communication.unionCommGet(),&rank);

    if (commTypeGet() == CWP_COMM_PAR_WITH_PART || 
        (commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && rank == _communication.unionCommLocCodeRootRanksGet())) {
      _visu.VisuCreate(freq,
                     format,
                     format_option,
                     output_dir,
                     "chr");

      _visu.GeomCreate(_mesh.getNPart());
    }
  }

  CWP_g_num_t*
  Coupling::globalNumGet(int id_block,int i_part) 
  {
    return _mesh.globalNumGet(id_block,i_part);
  }

 void Coupling::recvNextTimeSet (double next_time) 
 {

   if(_visu.isCreated() and _visu.physicalTimeGet() > -1.0) {
       _visu.WriterStepEnd();
   }

   _recvNextTime = next_time;

   if(_visu.isCreated()) {
       _visu.WriterStepBegin(_recvNextTime,&_mesh);
   }


 }


} // namespace cwipi

/**
 * \endcond
 */
