/*  This file is part of the CWIPI library.

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

#include <pdm_error.h>

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


  /**
   * \brief Constructor.
   *
   * This function creates a coupling object and defines its properties.
   *
   * \param [in]  cplId                        Coupling identifier
   * \param [in]  commType                     Communication type
   * \param [in]  localCodeProperties          Local code properties
   * \param [in]  coupledCodeProperties        Coupled code properties
   * \param [in]  spatialInterpAlgo                     SpatialInterp algorithm
   * \param [in]  nPart                        Number of interface partitions
   * \param [in]  movingStatus                 Mesh moving status
   * \param [in]  recvFreqType                 Type of receiving frequency
   * \param [in]  cplDB                        Coupling data base where it coupling is stored
   *
   */

  Coupling::Coupling
  (
   const string               &cplId,
   const CWP_Comm_t           cplType,
   const CodeProperties       &localCodeProperties,
   const CodeProperties       &coupledCodeProperties,
   const CWP_Interface_t      entities_dim,
   const CWP_Spatial_interp_t spatialInterpAlgo,
   const int                  nPart,
   const CWP_Dynamic_mesh_t   displacement,
   const CWP_Time_exch_t      recvFreqType,
   CouplingDB                 &cplDB
   )
  :_cplId(cplId),
   _commType(cplType),
   _communication(*(FC::getInstance().CreateObject(cplType))),
   _localCodeProperties(localCodeProperties),
   _coupledCodeProperties(coupledCodeProperties),
   _entities_dim(entities_dim),
   _spatial_interp_send(*new std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>()),
   _spatial_interp_recv(*new std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>()),
   _mesh(*new Mesh(localCodeProperties.connectableCommGet(),NULL,nPart,displacement,this)),
   _recvFreqType (recvFreqType),
   _visu(*new Visu(localCodeProperties.connectableCommGet(),displacement)),
   _fields(*(new map < string, Field * >())),
   _cplDB(cplDB),
   _iteration(new int),
   _displacement(displacement),
   _spatialInterpAlgo(spatialInterpAlgo),
   _userTargetN(nullptr),
   _userTargetGnum(nullptr),
   _localUserTargetGnum(nullptr),
   _userTargetCoord(nullptr),
   _nPart(nPart)
  {

/*    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
*/
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

    } // end else

    //entitiesDimGet();

  }


  /**
   * \brief Destructor.
   *
   */

  Coupling::~Coupling()
  {

    if(_visu.isCreated()) {
       _visu.WriterStepEnd();
    }

    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>::iterator it = _spatial_interp_send.begin();
    while (it != _spatial_interp_send.end()) {
        delete it -> second;
        it++;
    }

    it = _spatial_interp_recv.begin();
    while (it != _spatial_interp_recv.end()) {
       delete it -> second;
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

    if (_userTargetN != nullptr) {
      if (_localUserTargetGnum != nullptr) {
        for (int iPart; iPart < _nPart; iPart++) {
          free (_localUserTargetGnum[iPart]); 
        }
      }
      delete [] _userTargetN;        
      delete [] _userTargetGnum;     
      delete [] _userTargetCoord;    
    }

    #if defined(DEBUG) && 0
    cout << "destroying '" << _name << "' coupling : TODO" << endl;
    #endif
  }


  /*----------------------------------------------------------------------------*
   * Methods about exchange frequency                                           *
   *----------------------------------------------------------------------------*/




  /**
   * \brief Setting the next receiving time.
   *
   * This function set the next receiving time. It must be used when
   * the type of receiving frequency is \ref CWP_TIME_EXCH_ASYNCHRONOUS
   *
   * \param [in]  next_time     Next receiving time
   *
   */

  // void 
  // Coupling::recvNextTimeSet (
  //   double next_time
  // ) 
  // {
  //   PDM_UNUSED (next_time);
  //   PDM_error(__FILE__, __LINE__, 0, "\nrecvNextTimeSet not implemented yet\n");    
  // }

 
  /*----------------------------------------------------------------------------*
   * Methods about spatial interpolation                                        *
   *----------------------------------------------------------------------------*/


  /**
   * \brief Computation spatial interpolation weights
   *
   * This function compute spatial interpolation weights
   *
   */

  void
  Coupling::spatialInterpWeightsCompute ()
  {

    /////////////////////////////////////////////////////////////////////////////
    //                                                                         //
    // Exchange fields properties to obtain spatial intepolation to build      // 
    //                                                                         //
    /////////////////////////////////////////////////////////////////////////////

    printf("Coupling::spatialInterpWeightsCompute - 1\n");
    fflush(stdout);

    // - Store Data to send 

    if (!_coupledCodeProperties.localCodeIs()) {

      printf("Coupling::spatialInterpWeightsCompute - not cpl local code is - 1\n");
      fflush(stdout);

      std::string localFieldName=""; 
      vector<int> localFieldNameIdx;

      int localNbField = (int) _fields.size();

      localFieldNameIdx.reserve(localNbField + 1);
      localFieldNameIdx.push_back(0);
    
      std::vector<CWP_Field_exch_t> localFieldExch;
      localFieldExch.reserve(localNbField);

      std::vector<CWP_Dof_location_t> localFieldLocationV;
      localFieldLocationV.reserve(localNbField);
    
      std::map <std::string, cwipi::Field *>::iterator it = _fields.begin();

      while(it != _fields.end()){
        cwipi::Field* field = it -> second;
    
        localFieldName += it->first;
        localFieldNameIdx.push_back((int) (localFieldNameIdx[localFieldNameIdx.size()-1]+it->first.size()));
        localFieldLocationV.push_back(field->locationGet());
        localFieldExch.push_back(field->exchangeTypeGet());

        it++;
      }

      printf("Coupling::spatialInterpWeightsCompute - not cpl local code is - 2\n");
      fflush(stdout);

      // - Exchange number of fields 

      int nSendData   = 1;
      int nRecvData   = 1;
      int cplNbField = 0;

      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                 nSendData,
                                                                 (void *) &localNbField,
                                                                 -1,
                                                                 NULL,
                                                                 nRecvData,
                                                                 (void *) &cplNbField,
                                                                 -1,
                                                                 NULL);
      // - Allocate memory to receive data 

      vector<int               > cplFieldNameIdx (cplNbField + 1, 0);
      vector<CWP_Field_exch_t  > cplFieldExch (cplNbField);
      vector<CWP_Dof_location_t> cplFieldLocationV (cplNbField);
      string                     cplFieldName;

      // - Transfer memory to receive data 

      printf("Coupling::spatialInterpWeightsCompute - localNbField, cplNbField : <%d> <%d>\n", 
             localNbField, cplNbField);
      fflush(stdout);
//      abort();

      nSendData   = localNbField + 1;
      nRecvData   = cplNbField + 1;
      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                 nSendData,
                                                                 (void *) &(localFieldNameIdx[0]),
                                                                 -1,
                                                                 NULL,
                                                                 nRecvData,
                                                                 (void *) &(cplFieldNameIdx[0]),
                                                                 -1,
                                                                 NULL);
      localFieldNameIdx.clear();
      printf("Coupling::spatialInterpWeightsCompute - cplFieldNameIdx :");
      for (int j = 0; j < cplNbField + 1; j++) {
        printf(" %d", localFieldNameIdx[j]);
      }
      printf("\n");
      fflush(stdout);

      cplFieldName.resize(cplFieldNameIdx[cplNbField]);
      nSendData   = localFieldNameIdx[localNbField];
      nRecvData   = cplFieldNameIdx[cplNbField];

      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                 nSendData,
                                                                 (void *) localFieldName.c_str(),
                                                                 -1,
                                                                 NULL,
                                                                 nRecvData,
                                                                 (void *) cplFieldName.c_str(),
                                                                 -1,
                                                                 NULL);
      printf("Coupling::spatialInterpWeightsCompute - field name : <%s> <%s>\n", 
             localFieldName.c_str(), cplFieldName.c_str());
      fflush(stdout);
      localFieldName.clear();

      nSendData   = localNbField;
      nRecvData   = cplNbField;
      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Field_exch_t),
                                                                 nSendData,
                                                                 (void *) &(localFieldExch[0]),
                                                                 -1,
                                                                 NULL,
                                                                 nRecvData,
                                                                 (void *) &(cplFieldExch[0]),
                                                                 -1,
                                                                 NULL);
      localFieldExch.clear();

      nSendData   = localNbField;
      nRecvData   = cplNbField;
      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                 nSendData,
                                                                 (void *) &(localFieldLocationV[0]),
                                                                 -1,
                                                                 NULL,
                                                                 nRecvData,
                                                                 (void *) &(cplFieldLocationV[0]),
                                                                 -1,
                                                                 NULL);
      localFieldLocationV.clear();

      // Create spatial interpolation objects

      printf("Coupling::spatialInterpWeightsCompute - not cpl local code is - 3\n");
      fflush(stdout);
      
      it = _fields.begin();
      while(it != _fields.end()) {
        cwipi::Field* field = it -> second;
        string localFieldName = it -> first;
        CWP_Field_exch_t   localFieldExch     = field->exchangeTypeGet();
        CWP_Dof_location_t localFieldLocation = field->locationGet();
        
        for (int j = 0; j < cplNbField; j++) {
          string _cplFieldName = cplFieldName.substr(cplFieldNameIdx[j], cplFieldNameIdx[j+1]-cplFieldNameIdx[j]);
          printf("cplFieldName : %s %d %d\n", _cplFieldName.c_str(), cplFieldNameIdx[j], cplFieldNameIdx[j+1]);
          if (_cplFieldName == localFieldName) {
            field->linkedFieldLocationSet(cplFieldLocationV[j]);
            printf("find field\n");
            if (  (localFieldExch == CWP_FIELD_EXCH_SENDRECV && cplFieldExch[j] == CWP_FIELD_EXCH_SENDRECV)
                ||(localFieldExch == CWP_FIELD_EXCH_SEND     && cplFieldExch[j] == CWP_FIELD_EXCH_RECV)
                ||(localFieldExch == CWP_FIELD_EXCH_RECV     && cplFieldExch[j] == CWP_FIELD_EXCH_SEND)) {
    
              if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_SEND) {
                std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]); 
                if (_spatial_interp_send.find(newKey) == _spatial_interp_send.end()) {
                  printf("Coupling::spatialInterpWeightsCompute - add spatial send <%d> <%d>\n", 
                         localFieldLocation, cplFieldLocationV[j]);
                  fflush(stdout);
                  _spatial_interp_send.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                  printf("Coupling::spatialInterpWeightsCompute - init spatial send <%d> <%d>\n", 
                         localFieldLocation, cplFieldLocationV[j]);
                  fflush(stdout);
                  _spatial_interp_send[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_SEND);
                }
              }

              if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_RECV) {
                std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]); 
                if (_spatial_interp_recv.find(newKey) == _spatial_interp_recv.end()) {
                  printf("Coupling::spatialInterpWeightsCompute - add spatial recv <%d> <%d>\n", 
                         localFieldLocation, cplFieldLocationV[j]);
                  fflush(stdout);
                  _spatial_interp_recv.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                  printf("Coupling::spatialInterpWeightsCompute - init spatial recv <%d> <%d>\n", 
                         localFieldLocation, cplFieldLocationV[j]);
                  fflush(stdout);
                  _spatial_interp_recv[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_RECV);
                }
              }
            }
          }
        }
        it++;
      }

      printf("_spatial_interp_send.size() %d\n", (int) _spatial_interp_send.size());
      printf("_spatial_interp_recv.size() %d\n", (int) _spatial_interp_recv.size());

      printf("Coupling::spatialInterpWeightsCompute - not cpl local code is - 4\n");
      fflush(stdout);

      int codeID    = localCodePropertiesGet()->idGet();
      int cplCodeID = coupledCodePropertiesGet()->idGet();

      int sis_s = (int) _spatial_interp_send.size();

      vector<CWP_Dof_location_t> sis_loc;
      sis_loc.reserve(2*sis_s);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
      while(sis_it != _spatial_interp_send.end()) {
        sis_loc.push_back((sis_it->first).first);
        sis_loc.push_back((sis_it->first).second);
        sis_it++;
      }

      int sis_r;
      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                 1,
                                                                 (void *) &sis_s,
                                                                 -1,
                                                                 NULL,
                                                                 1,
                                                                 (void *) &sis_r,
                                                                 -1,
                                                                 NULL);

      int sir_s = (int) _spatial_interp_recv.size();
      assert(sis_r == sir_s);

      vector<CWP_Dof_location_t> sir_loc;
      sir_loc.reserve(2*sir_s);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sir_it = _spatial_interp_recv.begin();
      while(sir_it != _spatial_interp_recv.end()) {
        sir_loc.push_back((sir_it->first).first);
        sir_loc.push_back((sir_it->first).second);
        sir_it++;
      }

      int sir_r;
      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                 1,
                                                                 (void *) &sir_s,
                                                                 -1,
                                                                 NULL,
                                                                 1,
                                                                 (void *) &sir_r,
                                                                 -1,
                                                                 NULL);

      assert(sir_r == sis_s);

      vector<CWP_Dof_location_t> cpl_sis_loc;
      cpl_sis_loc.reserve(2*sir_s);

      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                 2*sis_s,
                                                                 (void *) &(sis_loc[0]),
                                                                 -1,
                                                                 NULL,
                                                                 2*sis_r,
                                                                 (void *) &(cpl_sis_loc[0]),
                                                                 -1,
                                                                 NULL);

      vector<CWP_Dof_location_t> cpl_sir_loc;
      cpl_sir_loc.reserve(2*sis_s);

      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                 2*sir_s,
                                                                 (void *) &(sir_loc[0]),
                                                                 -1,
                                                                 NULL,
                                                                 2*sir_r,
                                                                 (void *) &(cpl_sir_loc[0]),
                                                                 -1,
                                                                 NULL);

      printf("Coupling::spatialInterpWeightsCompute - not cpl local code is - 5\n");
      fflush(stdout);

      if (codeID < cplCodeID) {

        // spatial_interp send 

        sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          sis_it->second->weightsCompute();
          sis_it++;
        }

        // spatial_interp recv 

        for (int i = 0; i < sir_s; i++) {
          _spatial_interp_recv[make_pair(cpl_sis_loc[2*i+1], cpl_sis_loc[2*i])]->weightsCompute();  
        }

      printf("Coupling::spatialInterpWeightsCompute - not cpl local code is - 6\n");
      fflush(stdout);
      }

      else {

        // spatial_interp recv 

        for (int i = 0; i < sir_s; i++) {
          _spatial_interp_recv[make_pair(cpl_sis_loc[2*i+1], cpl_sis_loc[2*i])]->weightsCompute();  
        }

        // spatial_interp send 

        sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          sis_it->second->weightsCompute();
          sis_it++;
        }
      printf("Coupling::spatialInterpWeightsCompute - not cpl local code is - 7\n");
      fflush(stdout);

      }

    }

    // if an instance of each code on the same rank. All work is made in the call from the smallest code id

    else { 
    
      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 1\n");
      fflush(stdout);

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        std::string localFieldName=""; 
        vector<int> localFieldNameIdx;

        int localNbField = (int) _fields.size();

        localFieldNameIdx.reserve(localNbField + 1);
        localFieldNameIdx.push_back(0);
      
        std::vector<CWP_Field_exch_t> localFieldExch;
        localFieldExch.reserve(localNbField);

        std::vector<CWP_Dof_location_t> localFieldLocationV;
        localFieldLocationV.reserve(localNbField);
      
        std::map <std::string, cwipi::Field *>::iterator it = _fields.begin();

        while(it != _fields.end()){
          cwipi::Field* field = it -> second;
      
          localFieldName += it->first;
          localFieldNameIdx.push_back((int) (localFieldNameIdx[localFieldNameIdx.size() - 1] + it->first.size()));
          localFieldLocationV.push_back(field->locationGet());
          localFieldExch.push_back(field->exchangeTypeGet());

          it++;
        }

        std::string cpl_localFieldName=""; 
        vector<int> cpl_localFieldNameIdx;

        int cpl_localNbField = (int) cpl_cpl._fields.size();

        printf("cpl_localNbField %d\n", cpl_localNbField);
        printf("localNbField %d\n", localNbField);

        cpl_localFieldNameIdx.reserve(cpl_localNbField + 1);
        cpl_localFieldNameIdx.push_back(0);

        std::vector<CWP_Field_exch_t> cpl_localFieldExch;
        cpl_localFieldExch.reserve(cpl_localNbField);

        std::vector<CWP_Dof_location_t> cpl_localFieldLocationV;
        cpl_localFieldLocationV.reserve(cpl_localNbField);

        it = cpl_cpl._fields.begin();

        while(it != cpl_cpl._fields.end()){
          cwipi::Field* field = it -> second;

          cpl_localFieldName += it->first;
          cpl_localFieldNameIdx.push_back((int) (cpl_localFieldNameIdx[cpl_localFieldNameIdx.size()-1]+it->first.size()));
          cpl_localFieldLocationV.push_back(field->locationGet());
          cpl_localFieldExch.push_back(field->exchangeTypeGet());

          it++;
        }

        // - Exchange number of fields

        int nSendData       = 1;
        int nRecvData       = 1;
        int cplNbField      = 0;

        int cpl_nSendData       = 1;
        int cpl_nRecvData       = 1;
        int cpl_cplNbField      = 0;

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   nSendData,
                                                                   (void *) &localNbField,
                                                                   cpl_nSendData,
                                                                   (void *) &cpl_localNbField,
                                                                   nRecvData,
                                                                   (void *) &cplNbField,
                                                                   cpl_nRecvData,
                                                                   (void *) &cpl_cplNbField);
        // - Allocate memory to receive data 

        printf("cplNbField %d\n", cplNbField);
        printf("cpl_cplNbField %d\n", cpl_cplNbField);

        vector<int               > cplFieldNameIdx (cplNbField + 1, 0);
        vector<CWP_Field_exch_t  > cplFieldExch (cplNbField);
        vector<CWP_Dof_location_t> cplFieldLocationV (cplNbField);
        string                     cplFieldName;

        vector<int               > cpl_cplFieldNameIdx (cpl_cplNbField + 1, 0);
        vector<CWP_Field_exch_t  > cpl_cplFieldExch (cpl_cplNbField);
        vector<CWP_Dof_location_t> cpl_cplFieldLocationV (cpl_cplNbField);
        string                     cpl_cplFieldName;

        // - Transfer memory to receive data 

        nSendData   = localNbField + 1;
        nRecvData   = cplNbField + 1;

        cpl_nSendData   = cpl_localNbField + 1;
        cpl_nRecvData   = cpl_cplNbField + 1;

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   nSendData,
                                                                   (void *) &(localFieldNameIdx[0]),
                                                                   cpl_nSendData,
                                                                   (void *) &(cpl_localFieldNameIdx[0]),
                                                                   nRecvData,
                                                                   (void *) &(cplFieldNameIdx[0]),
                                                                   cpl_nRecvData,
                                                                   (void *) &(cpl_cplFieldNameIdx[0]));
        localFieldNameIdx.clear();

        cpl_localFieldNameIdx.clear();

        cplFieldName.resize(cplFieldNameIdx[cplNbField]);
        nSendData   = localFieldNameIdx[localNbField];
        nRecvData   = cplFieldNameIdx[cplNbField];

        cpl_cplFieldName.resize(cpl_cplFieldNameIdx[cpl_cplNbField]);
        cpl_nSendData   = cpl_localFieldNameIdx[cpl_localNbField];
        cpl_nRecvData   = cpl_cplFieldNameIdx[cpl_cplNbField];

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   nSendData,
                                                                   (void *) localFieldName.c_str(),
                                                                   cpl_nSendData,
                                                                   (void *) cpl_localFieldName.c_str(),
                                                                   nRecvData,
                                                                   (void *) cplFieldName.c_str(),
                                                                   cpl_nRecvData,
                                                                   (void *) cpl_cplFieldName.c_str());

        localFieldName.clear();
        cpl_localFieldName.clear();

        nSendData   = localNbField;
        nRecvData   = cplNbField;

        cpl_nSendData   = cpl_localNbField;
        cpl_nRecvData   = cpl_cplNbField;

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Field_exch_t),
                                                                   nSendData,
                                                                   (void *) &(localFieldExch[0]),
                                                                   cpl_nSendData,
                                                                   (void *) &(cpl_localFieldExch[0]),
                                                                   nRecvData,
                                                                   (void *) &(cplFieldExch[0]),
                                                                   cpl_nRecvData,
                                                                   (void *) &(cpl_cplFieldExch[0]));

        localFieldExch.clear();
        cpl_localFieldExch.clear();

        nSendData   = localNbField;
        nRecvData   = cplNbField;

        cpl_nSendData   = cpl_localNbField;
        cpl_nRecvData   = cpl_cplNbField;

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                   nSendData,
                                                                   (void *) &(localFieldLocationV[0]),
                                                                   cpl_nSendData,
                                                                   (void *) &(cpl_localFieldLocationV[0]),
                                                                   nRecvData,
                                                                   (void *) &(cplFieldLocationV[0]),
                                                                   cpl_nRecvData,
                                                                   (void *) &(cpl_cplFieldLocationV[0]));

        localFieldLocationV.clear();
        cpl_localFieldLocationV.clear();

        // Create spatial interpolation objects

        printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 2\n");
        fflush(stdout);

        it = _fields.begin();
        while(it != _fields.end()) {
          cwipi::Field* field = it -> second;
          string localFieldName = it -> first;
          CWP_Field_exch_t   localFieldExch     = field->exchangeTypeGet();
          CWP_Dof_location_t localFieldLocation = field->locationGet();

          printf("cplNbField %d\n", cplNbField);

          for (int j = 0; j < cplNbField; j++) {
            printf("step1 cplFieldName cpl_cplFieldName localFieldName <%s> <%s> <%s>\n", cplFieldName.c_str(), cpl_cplFieldName.c_str(), localFieldName.c_str());

            // TODO What is the puropose of this sustr ?
//            string cplFieldName = cplFieldName.substr( cplFieldNameIdx[j], cplFieldNameIdx[j+1]-cplFieldNameIdx[j] );

            printf("step2 cplFieldName cpl_cplFieldName localFieldName <%s> <%s> <%s>\n", cplFieldName.c_str(), cpl_cplFieldName.c_str(), localFieldName.c_str());

            if (cplFieldName == localFieldName) {
              printf("step3 cplFieldName cpl_cplFieldName localFieldName <%s> <%s> <%s>\n", cplFieldName.c_str(), cpl_cplFieldName.c_str(), localFieldName.c_str());

              field->linkedFieldLocationSet(cplFieldLocationV[j]);

              printf("step4 cplFieldName cpl_cplFieldName localFieldName <%s> <%s> <%s>\n", cplFieldName.c_str(), cpl_cplFieldName.c_str(), localFieldName.c_str());

              if (cpl_cpl._fields.find(cplFieldName) == cpl_cpl._fields.end()) {
                PDM_error(__FILE__, __LINE__, 0, "'%s' Field not found\n",
                          cplFieldName.c_str());
              }

              cpl_cpl._fields[cplFieldName]->linkedFieldLocationSet(localFieldLocation);

              printf("step5 cplFieldName cpl_cplFieldName localFieldName <%s> <%s> <%s>\n", cplFieldName.c_str(), cpl_cplFieldName.c_str(), localFieldName.c_str());

              printf("cpl_cpl._fields[\"cplFieldName\"] %s\n", cpl_cpl._fields[cplFieldName]->fieldIDGet().c_str());
              fflush(stdout);

              printf("_spatial_interp_send.size() %d\n", (int) _spatial_interp_send.size());
              printf("_spatial_interp_recv.size() %d\n", (int) _spatial_interp_recv.size());

              printf("cpl_cpl._spatial_interp_send.size() %d\n", (int) cpl_cpl._spatial_interp_send.size());
              printf("cpl_cpl._spatial_interp_recv.size() %d\n", (int) cpl_cpl._spatial_interp_recv.size());

              if (  (localFieldExch == CWP_FIELD_EXCH_SENDRECV && cplFieldExch[j] == CWP_FIELD_EXCH_SENDRECV)
                  ||(localFieldExch == CWP_FIELD_EXCH_SEND     && cplFieldExch[j] == CWP_FIELD_EXCH_RECV)
                  ||(localFieldExch == CWP_FIELD_EXCH_RECV     && cplFieldExch[j] == CWP_FIELD_EXCH_SEND)) {
      
                if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_SEND) {
                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]); 
                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocationV[j], localFieldLocation); 

                  if (_spatial_interp_send.find(newKey) == _spatial_interp_send.end()) {
                    printf("Creating send spatialInterp\n");

                    _spatial_interp_send.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                    cpl_cpl._spatial_interp_recv.insert(make_pair(cpl_newKey, FG::getInstance().CreateObject(cpl_cpl._spatialInterpAlgo)));
                    _spatial_interp_send[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_SEND);
                    cpl_cpl._spatial_interp_recv[cpl_newKey]->init(&cpl_cpl, cplFieldLocationV[j], localFieldLocation, SPATIAL_INTERP_EXCH_RECV);
                  }
                }

                if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_RECV) {
                  printf("Entering recv\n");

                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]);
                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocationV[j], localFieldLocation);

                  if (_spatial_interp_recv.find(newKey) == _spatial_interp_recv.end()) {
                    printf("Creating recv spatialInterp\n");

                    _spatial_interp_recv.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                    cpl_cpl._spatial_interp_send.insert(make_pair(cpl_newKey, FG::getInstance().CreateObject(cpl_cpl._spatialInterpAlgo)));
                    _spatial_interp_recv[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_RECV);
                    cpl_cpl._spatial_interp_send[cpl_newKey]->init(&cpl_cpl, cplFieldLocationV[j], localFieldLocation, SPATIAL_INTERP_EXCH_SEND);
                  }
                }
              }
            }
          }
          it++;
        }  

        printf("cpl_cpl._fields.size() %d\n", (int) cpl_cpl._fields.size());

        printf("_spatial_interp_send.size() %d\n", (int) _spatial_interp_send.size());
        printf("_spatial_interp_recv.size() %d\n", (int) _spatial_interp_recv.size());
        printf("cpl_cpl._spatial_interp_send.size() %d\n", (int) cpl_cpl._spatial_interp_send.size());
        printf("cpl_cpl._spatial_interp_recv.size() %d\n", (int) cpl_cpl._spatial_interp_recv.size());

        printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 3\n");
        fflush(stdout);
        int codeID    = localCodePropertiesGet()->idGet();
        int cplCodeID = coupledCodePropertiesGet()->idGet();

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
        int sis_s = (int) _spatial_interp_send.size();
        printf("sis_s %d\n", sis_s);

        vector<CWP_Dof_location_t> sis_loc;
        sis_loc.reserve(2*sis_s);

        while(sis_it != _spatial_interp_send.end()) {
          sis_loc.push_back((sis_it->first).first);
          sis_loc.push_back((sis_it->first).second);
          printf("1 sis_it++\n");
          fflush(stdout);
          sis_it++;
        }

        int cpl_sis_s = (int) cpl_cpl._spatial_interp_send.size();
        printf("cpl_sis_s %d\n", cpl_sis_s);

        vector<CWP_Dof_location_t> cpl_sis_loc;
        cpl_sis_loc.reserve(2*cpl_sis_s);

        sis_it = cpl_cpl._spatial_interp_send.begin();
        while(sis_it != cpl_cpl._spatial_interp_send.end()) {
          cpl_sis_loc.push_back((sis_it->first).first);
          cpl_sis_loc.push_back((sis_it->first).second);
          printf("2 sis_it++\n");
          fflush(stdout);
          sis_it++;
        }

      printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 4\n");
      fflush(stdout);

        int sis_r;
        int cpl_sis_r;
        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   1,
                                                                   (void *) &sis_s,
                                                                   1,
                                                                   (void *) &cpl_sis_s,
                                                                   1,
                                                                   (void *) &sis_r,
                                                                   1,
                                                                   (void *) &cpl_sis_r);
      printf("sis_r %d\n", sis_r);
      printf("cpl_sis_r %d\n", cpl_sis_r);

      vector<CWP_Dof_location_t> sis_loc_r;
      sis_loc_r.reserve(2*sis_r);
      vector<CWP_Dof_location_t> cpl_sis_loc_r;
      cpl_sis_loc_r.reserve(2*cpl_sis_r);

      printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 41\n");
      fflush(stdout);

      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                 2*sis_s,
                                                                 (void *) &(sis_loc[0]),
                                                                 2*cpl_sis_s,
                                                                 (void *) &(cpl_sis_loc[0]),
                                                                 2*sis_r,
                                                                 (void *) &(sis_loc_r[0]),
                                                                 2*cpl_sis_r,
                                                                 (void *) &(cpl_sis_loc_r[0]));

       printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 42\n");
       fflush(stdout);

      auto sir_it = _spatial_interp_recv.begin();
      int sir_s = (int) _spatial_interp_recv.size();
      printf("sir_s %d\n", sir_s);

      vector<CWP_Dof_location_t> sir_loc;
      sir_loc.reserve(2*sis_r);

      while(sir_it != _spatial_interp_recv.end()) {
        sir_loc.push_back((sir_it->first).first);
        sir_loc.push_back((sir_it->first).second);
        printf("1 sir_it++\n");
        fflush(stdout);
        sir_it++;
      }

      int cpl_sir_s = (int) cpl_cpl._spatial_interp_recv.size();
      printf("cpl_sir_s %d\n", cpl_sir_s);

      vector<CWP_Dof_location_t> cpl_sir_loc;
      cpl_sir_loc.reserve(2*cpl_sir_s);

      sir_it = cpl_cpl._spatial_interp_recv.begin();
      while(sir_it != cpl_cpl._spatial_interp_recv.end()) {
        cpl_sir_loc.push_back((sir_it->first).first);
        cpl_sir_loc.push_back((sir_it->first).second);
        printf("2 sir_it++\n");
        fflush(stdout);
        sir_it++;
      }

      printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 5\n");
      fflush(stdout);

      int sir_r;
      int cpl_sir_r;
      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                 1,
                                                                 (void *) &sir_s,
                                                                 1,
                                                                 (void *) &cpl_sir_s,
                                                                 1,
                                                                 (void *) &sir_r,
                                                                 1,
                                                                 (void *) &cpl_sir_r);
      printf("sir_s %d\n", sir_r);
      printf("cpl_sir_s %d\n", cpl_sir_r);

      vector<CWP_Dof_location_t> sir_loc_r;
      sir_loc_r.reserve(2*sir_r);
      vector<CWP_Dof_location_t> cpl_sir_loc_r;
      cpl_sir_loc_r.reserve(2*cpl_sir_r);

      printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 51\n");
      fflush(stdout);

      _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                 2*sir_s,
                                                                 (void *) &(sir_loc[0]),
                                                                 2*cpl_sir_s,
                                                                 (void *) &(cpl_sir_loc[0]),
                                                                 2*sir_r,
                                                                 (void *) &(sir_loc_r[0]),
                                                                 2*cpl_sir_r,
                                                                 (void *) &(cpl_sir_loc_r[0]));

      assert(sir_r     == sis_s);
      assert(sis_r     == sir_s);
      assert(cpl_sir_r == cpl_sis_s);
      assert(cpl_sis_r == cpl_sir_s);

      printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 6\n");
      fflush(stdout);

      printf("codeID cplCodeID %d %d\n", codeID, cplCodeID);

      if (codeID < cplCodeID) {
        // spatial_interp send
        printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 61\n");
        fflush(stdout);

        sis_it = _spatial_interp_send.begin();
        int i = 0;
        while(sis_it != _spatial_interp_send.end()) {
          sis_it->second->weightsCompute();
          cpl_cpl._spatial_interp_recv[make_pair(cpl_sis_loc_r[2 * i + 1], cpl_sis_loc_r[2 * i])]->weightsCompute();
          sis_it++;
          i++;
        }

        printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 611\n");
        fflush(stdout);

        sis_it = cpl_cpl._spatial_interp_send.begin();
        i = 0;
        while(sis_it != cpl_cpl._spatial_interp_send.end()) {
          _spatial_interp_recv[make_pair(sis_loc_r[2 * i + 1], sis_loc_r[2 * i])]->weightsCompute();
          sis_it->second->weightsCompute();
          sis_it++;
          i++;
        }
      }
      else {
        // spatial_interp recv
        printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 62\n");
        fflush(stdout);

        int i = 0;
        sis_it = cpl_cpl._spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          _spatial_interp_recv[make_pair(sis_loc_r[2*i+1], sis_loc_r[2*i])]->weightsCompute();
          sis_it->second->weightsCompute();
          sis_it++;
          i++;
        }

        // spatial_interp send
        sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          sis_it->second->weightsCompute();
          cpl_cpl._spatial_interp_recv[make_pair(cpl_sis_loc_r[2*i+1], cpl_sis_loc_r[2*i])]->weightsCompute();
          sis_it++;
          i++;
        }
      }
    }
      else {
        printf("Coupling::spatialInterpWeightsCompute - cpl local code is - 7\n");
      }
    }

      printf("Coupling::spatialInterpWeightsCompute - fin\n");
      fflush(stdout);
  }


  /*----------------------------------------------------------------------------*
   * Methods about visualization                                                *
   *----------------------------------------------------------------------------*/

  /**
   * \brief Enable visualization output
   *
   * This function enable visualization output.
   *
   * \param [in]  freq             Output frequency
   * \param [in]  format           Output format to visualize exchanged fieldsDouble
   *                               on the coupled mesh. Choice between :
   *                               - "EnSight Gold"
   *                               - "MED_ficher"
   *                               - "CGNS"
   *                               .
   * \param [in]  format_option   Output options "opt1, opt2, ..." :
   *                         - text               output text files
   *                         - binary             output binary files (default)
   *                         - big_endian         force binary files
   *                                              to big-endian
   *                         - discard_polygons   do not output polygons
   *                                              or related values
   *                         - discard_polyhedra  do not output polyhedra
   *                                              or related values
   *                         - divide_polygons    tesselate polygons
   *                                              with triangles
   *                         - divide_polyhedra   tesselate polyhedra
   *                                              with tetrahedra and pyramids
   *                                              (adding a vertex near
   *                                               each polyhedron's center)
   *                         .
   *
   */

  void 
  Coupling::visuSet (
    const int               freq,
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
                     (char*) "chr");

      _visu.GeomCreate(_mesh.getNPart());
    }
  }

  /*----------------------------------------------------------------------------*
   * Methods  about mesh                                                        *
   *----------------------------------------------------------------------------*/




  /*----------------------------------------------------------------------------*
   * Methods about field                                                        *
   *----------------------------------------------------------------------------*/

  /**
   *
   * \brief Create a new field
   *
   * \param [in]  field_id       Field id
   * \param [in]  data_type      Data type
   * \param [in]  storage        Storage type
   * \param [in]  n_component    Number of componenent
   * \param [in]  nature         Nature
   * \param [in]  exch_type      Exchange type
   * \param [in]  visu_status    Visualization status
   *
   */

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
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' existing field\n", field_id.c_str());
    }

    //
    // Create the new field

    double physTime=0.0;
    *_iteration = 0;

    printf("sizeof fields : %d\n", (int) _fields.size());

    cwipi::Field *newField = new cwipi::Field(field_id,
                                              _fields.size(),
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


  /**
   * \brief Return if a field identifier exists
   *
   * \param [in]  field_id         Field identifier
   *
   * \return status
   */


  bool
  Coupling::fieldIs
  (
   const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    return (It != _fields.end());
  }


 /**
  *
  * \brief Set Field data
  *
  * \param [in]  field_id       Field identifier
  * \param [in]  data           Storage array (mapping)
  *
  */

  void
  Coupling::fieldDataSet
  (
    const string &field_id,
    int i_part,
    const CWP_Field_map_t   map_type,
    void* data
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It->second->dataSet(i_part, map_type, data);
      if (_visu.isCreated() && It -> second -> visuStatusGet() == CWP_STATUS_ON) {
        _visu.fieldDataSet(It->second,map_type, i_part);
      }
    }
  }

  /*----------------------------------------------------------------------------*
   * Methods about exchange                                                     *
   *----------------------------------------------------------------------------*/


  /**
   * \brief data exchange <b>(Not implemented yet)</b>
   *
   * Exchange depending on exchange frequency
   *
   */

  void
  Coupling::exchange ()
  {
    PDM_error(__FILE__, __LINE__, 0, "\nexchange not implemented yet\n");
  }


  /**
   * \brief Exchange data field with the coupled code with blocking
   *        communications. <b>(Not implemented yet)</b>
   *
   * This function exchanges interpolated fieldsDouble between coupled codes.
   *
   * \warning  The size of tgt_field_id size is n_computed_tgt.
   *           If \f$ n\_uncomputed\_tgt \ne n\_tgt\_pts \f$,
   *           user himself must set values for uncomputed target points.
   *
   * \param [in]  src_field_id              Source field (NULL -> no sending)
   * \param [in]  tgt_field_id              Target field (NULL -> no receiving)
   * \param [in]  ptFortranInterpolationFct Fortran user interpolation (or NULL)
   * \param [out] n_uncomputed_tgt          Number of uncomputed target
   *
   * \return                                Exchange status
   *
   */

  void
  Coupling::sendrecv (
    const string &field_id
  )
  {
    PDM_UNUSED (field_id);
    PDM_error(__FILE__, __LINE__, 0, "\nsendrecv not implemented yet\n");
  }


  /**
   *
   * \brief Sending of data field to the coupled code with nonblocking
   *        communications.
   *
   * This function sends interpolated field to the coupled code.
   *
   * \param [in]  src_id                    Source field
   *
   */

  void
  Coupling::issend (
    const string &sendingFieldID
  )
  {

    if (!_coupledCodeProperties.localCodeIs()) {

      map <string, Field *>::iterator it;
      it = _fields.find(sendingFieldID);

      if (it != _fields.end()) {
        Field* sendingField = it->second;

        CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_send.find(newKey); 

        if (it2 == _spatial_interp_send.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->issend(sendingField);
        }
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
      }
    }

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(sendingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(sendingFieldID);

        Field* recvField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          recvField = it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* sendingField = it->second;

          CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);       

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;                   
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;                   

          it2 = _spatial_interp_send.find(newKey); 
          cpl_it2 = cpl_cpl._spatial_interp_recv.find(newKey); 

          if (cpl_it2 == cpl_cpl._spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->irecv(recvField);
          }

          if (it2 == _spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->issend(sendingField);
          }
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }
  }


  /**
   *
   * \brief Waiting of the end of exchange related to request.
   *
   * This function waits the end of exchange related to request
   * from \ref CWP_Issend
   *
   * \param [in] src_id                    Source field
   *
   */

  void
  Coupling::waitIssend (
    const string &sendingFieldID
  )
  {
    if (!_coupledCodeProperties.localCodeIs()) {
      map <string, Field *>::iterator it;
      it = _fields.find(sendingFieldID);

      if (it != _fields.end()) {
        Field* sendingField = it->second;

        CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_send.find(newKey); 

        if (it2 == _spatial_interp_send.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->waitIssend(sendingField);
        }
      }
    }
    else {
      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {
        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(sendingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(sendingFieldID);

        Field* recvField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          recvField = it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* sendingField = it->second;

          CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);       

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;                   
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;                   

          it2 = _spatial_interp_send.find(newKey); 
          cpl_it2 = cpl_cpl._spatial_interp_recv.find(newKey); 

          if (cpl_it2 == cpl_cpl._spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->waitIrecv(recvField);
          }

          if (it2 == _spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->waitIssend(sendingField);
          }
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }
  }


  /**
   *
   * \brief Receiving of Data field from the coupled code with nonblocking
   *        communications.
   *
   * This function receives interpolated field from the coupled code
   *
   * \param [in]  receving_field_id       Target field ID
   *
   *
   */

  void
  Coupling::irecv
  (
    const string &receivingFieldID
  ) 
  {
    if (!_coupledCodeProperties.localCodeIs()) {

      map <string, Field *>::iterator it;
      it = _fields.find(receivingFieldID);

      if (it != _fields.end()) {
        Field* receivingField = it->second;

        CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_recv.find(newKey); 

        if (it2 == _spatial_interp_recv.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->irecv(receivingField);
        }
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
      }
    }

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(receivingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(receivingFieldID);

        Field* sendField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          sendField = it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* receivingField = it->second;

          CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);       

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;                   
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;                   

          it2 = _spatial_interp_recv.find(newKey); 
          cpl_it2 = cpl_cpl._spatial_interp_send.find(newKey); 

          if (it2 == _spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->irecv(receivingField);
          }

          if (cpl_it2 == cpl_cpl._spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->issend(sendField);
          }

        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }

  }


  /**
   *
   * \brief Waiting of the end of exchange related to request.
   *
   * This function waits the end of exchange related to request
   * from \ref CWP_Irecv
   *
   * \param [in]  receving_field_id       Target field ID
   *
   */

  void
  Coupling::waitIrecv (
    const string &receivingFieldID
  )
  {
    if (!_coupledCodeProperties.localCodeIs()) {

      map <string, Field *>::iterator it;
      it = _fields.find(receivingFieldID);

      if (it != _fields.end()) {
        Field* receivingField = it->second;

        CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_recv.find(newKey); 

        if (it2 == _spatial_interp_recv.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->waitIrecv(receivingField);
        }
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
      }
    }

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(receivingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(receivingFieldID);

        Field* sendField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          sendField = it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* receivingField = it->second;

          CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation); 
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);       

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;                   
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;                   

          it2 = _spatial_interp_recv.find(newKey); 
          cpl_it2 = cpl_cpl._spatial_interp_send.find(newKey); 

          if (it2 == _spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->waitIrecv(receivingField);
          }

          if (cpl_it2 == cpl_cpl._spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->waitIssend(sendField);
          }

        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }
  }


  /*----------------------------------------------------------------------------*
   * methods about user interpolation                                           *
   *----------------------------------------------------------------------------*/


  /*----------------------------------------------------------------------------*
   * Private methods                                                            *
   *----------------------------------------------------------------------------*/

  /**
   *
   * \brief Compute user target global number (if not given by user)
   *
   */

  void 
  Coupling::userTargetGnumCompute() 
  {
    if (_userTargetN != nullptr) {
      if (_userTargetGnum == nullptr) {

        PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm (_localCodeProperties.intraCommGet());

        PDM_gen_gnum_t *pgg  = PDM_gnum_create (3, _nPart, PDM_FALSE, 1e-3, comm,
                                                   PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);
        for (int iPart = 0; iPart <_nPart; iPart++){
          PDM_gnum_set_from_coords (pgg, iPart, _userTargetN[iPart], _userTargetCoord[iPart], NULL);
        }

        PDM_gnum_compute (pgg);
  
        _localUserTargetGnum = new  CWP_g_num_t * [_nPart];

        for (int iPart = 0; iPart < _nPart; iPart++){
          _localUserTargetGnum[iPart] = const_cast <CWP_g_num_t*> (PDM_gnum_get (pgg, iPart));
        }

        PDM_gnum_free (pgg);

        _userTargetGnum = const_cast <const CWP_g_num_t**> (_localUserTargetGnum);
      }
    }
  }


// A supprimer

  CWP_g_num_t*
  Coupling::globalNumGet(int id_block,int i_part) 
  {
    return _mesh.globalNumGet(id_block,i_part);
  }



} // namespace cwipi

/**
 * \endcond
 */
