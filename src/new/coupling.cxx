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
#include "globalData.hxx"

#include "communication.hxx"
// #include "visualization.hxx"
#include "pdm_writer.h"
#include "pdm_logging.h"


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
         CodeProperties       &localCodeProperties,
         CodeProperties       &coupledCodeProperties,
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
   _mesh(*new Mesh(localCodeProperties.connectableCommGet(),nPart,displacement,this)),
   _recvFreqType (recvFreqType),
   _id_geom_writer(-1),
   _id_field_partitioning_writer(-1),
   _id_field_ranking_writer(-1),
   _freq_writer(-1),
   _writer(nullptr),
   _fields(*(new map < string, Field * >())),
   _globalData(*(new map < string, GlobalData * >())),
   _cplDB(cplDB),
   _displacement(displacement),
   _spatialInterpAlgo(spatialInterpAlgo),
   _nPart(nPart),
   _cplNPart(-1),
   _userTargetN(nullptr),
   _userTargetGnum(nullptr),
   _localUserTargetGnum(nullptr),
   _userTargetCoord(nullptr),
   _spatial_interp_send(*new std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>()),
   _spatial_interp_recv(*new std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>()),
   _spatial_interp_properties_double(*new std::map<std::string, double>),
   _spatial_interp_properties_int(*new std::map<std::string, int>),
   // _visu(*new Visu(localCodeProperties.connectableCommGet(),displacement)),
   _is_mesh_finalized(0),
   _is_first_field_created(0),
   _n_step(0)

   {

/*    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
*/
     //In case where the both codes are on the same MPI process.
    if (coupledCodeProperties.localCodeIs()) {
      if (cplDB.couplingIs(coupledCodeProperties, cplId) ) {
        const int localRootRank = localCodeProperties.rootRankGet();
        const int cplRootRank = coupledCodeProperties.rootRankGet();

        const MPI_Comm& globalComm = localCodeProperties.globalCommGet();

        int globalRank;
        MPI_Comm_rank(globalComm, &globalRank);

        Coupling &distCpl = cplDB.couplingGet(coupledCodeProperties, cplId);

        if ((cplRootRank == globalRank) && (cplRootRank != localRootRank)) {
          distCpl._communication.init(_coupledCodeProperties, _localCodeProperties, cplId, cplDB);
          _communication.init(distCpl._communication);
        }
        else {
          _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);
          distCpl._communication.init(_communication);
        }

        // Visu* visu_cpl = distCpl.visuGet();
        // Mesh* mesh_cpl = distCpl.meshGet();

        // _mesh.setVisu(&_visu);
        // mesh_cpl->setVisu(visu_cpl);

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

        // if(unionRank == _communication.unionCommLocCodeRootRanksGet()){
        //   _visu = *new Visu(visuComm,displacement);
        // }
      }

       // _mesh.setVisu(&_visu);

    } // end else

    //entitiesDimGet();

  }


  /**
   * \brief Destructor.
   *
   */

  Coupling::~Coupling()
  {

    delete &_spatial_interp_properties_double;
    delete &_spatial_interp_properties_int;

    // if(_visu.isCreated()) {
    //    _visu.WriterStepEnd();
    // }

    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>::iterator it = _spatial_interp_send.begin();
    while (it != _spatial_interp_send.end()) {
        delete it->second;
        it++;
    }

    delete &_spatial_interp_send;

    it = _spatial_interp_recv.begin();
    while (it != _spatial_interp_recv.end()) {
       delete it->second;
        it++;
    }

    delete &_spatial_interp_recv;


    std::map < string, Field * >::iterator itf = _fields.begin();
    while (itf != _fields.end()) {
      // if(_visu.isCreated() && itf->second->visuStatusGet() == CWP_STATUS_ON )
      //   _visu.fieldDataFree(itf->second);
      if (itf->second != NULL) delete itf->second;
      itf++;
    }

    delete &_fields;

    std::map < string, GlobalData * >::iterator itgd = _globalData.begin();
    while (itgd != _globalData.end()) {
      delete itgd->second;
      itgd++;
    }

    delete &_globalData;

    // if(_visu.isCreated()) {
    //   // _visu.SpatialInterpFree();
    // }

    // delete &_visu;

    delete &_communication;

    delete &_mesh;

    if (_userTargetN != nullptr) {
      if (_localUserTargetGnum != nullptr) {
        for (int iPart = 0; iPart < _nPart; iPart++) {
          free (_localUserTargetGnum[iPart]);
        }
      }
      free ( _userTargetN);
      free ( _userTargetGnum);
      free ( _userTargetCoord);
    }

    #if defined(DEBUG) && 0
    cout << "destroying '" << _name << "' coupling : TODO" << endl;
    #endif
  }

  /*----------------------------------------------------------------------------*
   * Methods about global data                                                  *
   *----------------------------------------------------------------------------*/

  /**
   * \brief Send a data array.
   *
   * \param [in] global_data_id
   * \param [in] s_send_entity
   * \param [in] send_stride
   * \param [in] n_send_entity
   * \param [in] send_data
   *
   */

  void
  Coupling::globalDataIsend
  (
   const string    &global_data_id,
   size_t          s_send_entity,
   int             send_stride,
   int             n_send_entity,
   void           *send_data
   )
  {
    cout << "globalDataIsend - in\n" << endl;

    // Create an instance of GlobalData
    map<string,GlobalData*>::iterator it = _globalData.find(global_data_id.c_str());
    if (it == _globalData.end()) {
      cwipi::GlobalData *newGlobalData = new cwipi::GlobalData(global_data_id,
                                                               s_send_entity,
                                                               send_stride,
                                                               n_send_entity,
                                                               send_data);

      pair<string, GlobalData* > newPair(global_data_id, newGlobalData);
      _globalData.insert(newPair);
    } // end if does not exist
    it = _globalData.find(global_data_id.c_str());

    assert(it != _globalData.end());

    size_t        &s_entity       = it->second->s_send_entity_get();
    int           &stride         = it->second->send_stride_get();
    int           &n_entity       = it->second->n_send_entity_get();
    void *        data           = it->second->send_data_get();

    MPI_Request s_entity_request;
    MPI_Request stride_request;
    MPI_Request n_entity_request;
    MPI_Request data_request;

    _communication.isendGlobalDataBetweenCodesThroughUnionCom(global_data_id,
                                                              &s_entity_request,
                                                              &stride_request,
                                                              &n_entity_request,
                                                              &data_request,
                                                              s_entity,
                                                              stride,
                                                              n_entity,
                                                              data);

    it->second->s_entity_request_set(s_entity_request);
    it->second->stride_request_set(stride_request);
    it->second->n_entity_request_set(n_entity_request);
    it->second->data_request_set(data_request);

    cout << "globalDataIsend - out\n" << endl;
  }

  /**
   * \brief Receive a data array.
   *
   * \param [in] global_data_id
   * \param [in] s_recv_entity
   * \param [in] recv_stride
   * \param [in] n_recv_entity
   * \param [in] recv_data
   *
   */

  void
  Coupling::globalDataIrecv
  (
   const string    &global_data_id,
   size_t         *s_recv_entity,
   int            *recv_stride,
   int            *n_recv_entity,
   void          **recv_data
  )
  {
    cout << "globalDataIrecv - in\n" << endl;

    assert(s_recv_entity != NULL);

    // Create an instance of GlobalData
    map<string,GlobalData*>::iterator it = _globalData.find(global_data_id.c_str());
    if (it == _globalData.end()) {
      cwipi::GlobalData *newGlobalData = new cwipi::GlobalData(global_data_id,
                                                               s_recv_entity,
                                                               recv_stride,
                                                               n_recv_entity,
                                                               recv_data);

      pair<string, GlobalData* > newPair(global_data_id, newGlobalData);
      _globalData.insert(newPair);
    } // end if does not exist
    it = _globalData.find(global_data_id.c_str());

    assert(it != _globalData.end());

    size_t * s_entity  = it->second->s_recv_entity_get();
    int *    stride    = it->second->recv_stride_get();
    int *    n_entity  = it->second->n_recv_entity_get();

    assert(s_entity != NULL);

    MPI_Request s_entity_request;
    MPI_Request stride_request;
    MPI_Request n_entity_request;

    _communication.irecvGlobalDataBetweenCodesThroughUnionCom(global_data_id,
                                                              &s_entity_request,
                                                              &stride_request,
                                                              &n_entity_request,
                                                              s_entity,
                                                              stride,
                                                              n_entity);

    it->second->s_entity_request_set(s_entity_request);
    it->second->stride_request_set(stride_request);
    it->second->n_entity_request_set(n_entity_request);

    cout << "globalDataIrecv - out\n" << endl;
  }

  /**
   * \brief Wait of send a data array.
   *
   * \param [in] global_data_id
   *
   */

  void
  Coupling::globalDataWaitIsend
  (
   const string    &global_data_id
  )
  {
    cout << "globalDataWaitIsend - in\n" << endl;

    // Get local
    map<string,GlobalData*>::iterator it = _globalData.find(global_data_id.c_str());
    MPI_Request s_entity_request = it->second->s_entity_request_get();
    MPI_Request stride_request   = it->second->stride_request_get();
    MPI_Request n_entity_request = it->second->n_entity_request_get();
    MPI_Request data_request     = it->second->data_request_get();
    size_t        s_entity         = it->second->s_send_entity_get();
    int           stride           = it->second->send_stride_get();
    int           n_entity         = it->second->n_send_entity_get();
    void *        data             = it->second->send_data_get();

    // Get coupled
    size_t *      cpl_s_entity = NULL;
    int *         cpl_stride   = NULL;
    int *         cpl_n_entity = NULL;
    void **       cpl_data     = NULL;
    if (_coupledCodeProperties.localCodeIs()) {
      cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);
      map<string,GlobalData*>::iterator cpl_it = cpl_cpl._globalData.find(global_data_id.c_str());
      cpl_s_entity       = cpl_it->second->s_recv_entity_get();
      cpl_stride         = cpl_it->second->recv_stride_get();
      cpl_n_entity       = cpl_it->second->n_recv_entity_get();
      cpl_data           = cpl_it->second->recv_data_get();
    }

    _communication.waitIsendGlobalDataBetweenCodesThroughUnionCom(global_data_id,
                                                                  &s_entity_request,
                                                                  &stride_request,
                                                                  &n_entity_request,
                                                                  &data_request,
                                                                  s_entity,
                                                                  stride,
                                                                  n_entity,
                                                                  data,
                                                                  cpl_s_entity,
                                                                  cpl_stride,
                                                                  cpl_n_entity,
                                                                  cpl_data);

    cout << "globalDataWaitIsend - out\n" << endl;
  }

  /**
   * \brief Wait of receive a data array.
   *
   * \param [in] global_data_id
   *
   */

  void
  Coupling::globalDataWaitIrecv
  (
   const string    &global_data_id
  )
  {
    cout << "globalDataWaitIrecv - in\n" << endl;

    // Get local
    map<string,GlobalData*>::iterator it = _globalData.find(global_data_id.c_str());
    MPI_Request  s_entity_request = it->second->s_entity_request_get();
    MPI_Request  stride_request   = it->second->stride_request_get();
    MPI_Request  n_entity_request = it->second->n_entity_request_get();
    size_t *      s_entity       = it->second->s_recv_entity_get();
    int *         stride         = it->second->recv_stride_get();
    int *         n_entity       = it->second->n_recv_entity_get();
    void **       data           = it->second->recv_data_get();

    // Get coupled
    size_t        cpl_s_entity = 0;
    int           cpl_stride   = 0;
    int           cpl_n_entity = 0;
    void *        cpl_data     = NULL;
    if (_coupledCodeProperties.localCodeIs()) {
      cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);
      map<string,GlobalData*>::iterator cpl_it = cpl_cpl._globalData.find(global_data_id.c_str());
      cpl_s_entity       = cpl_it->second->s_send_entity_get();
      cpl_stride         = cpl_it->second->send_stride_get();
      cpl_n_entity       = cpl_it->second->n_send_entity_get();
      cpl_data           = cpl_it->second->send_data_get();
    }

    _communication.waitIrecvGlobalDataBetweenCodesThroughUnionCom(global_data_id,
                                                                  &s_entity_request,
                                                                  &stride_request,
                                                                  &n_entity_request,
                                                                  cpl_s_entity,
                                                                  cpl_stride,
                                                                  cpl_n_entity,
                                                                  cpl_data,
                                                                  s_entity,
                                                                  stride,
                                                                  n_entity,
                                                                  data);

    cout << "globalDataWaitIrecv - out\n" << endl;
  }

  /*----------------------------------------------------------------------------*
   * Methods about communicators                                                *
   *----------------------------------------------------------------------------*/

 /**
  * \brief Get coupling communicator and coupling ranks.
  *
  * \param [out] cpl_comm             Coupling communicator
  * \param [out] cpl_ranks            Coupling ranks
  *
  * \return Size of \ref cpl_ranks vector
  *
  */

  int
  Coupling::commGet (
    MPI_Comm  *cpl_comm,
    int      **cpl_ranks
  )
  {
    *cpl_comm  = _communication.cplCommGet();

    std::vector<int>* vect_cpl_ranks = _communication.cplCommCplRanksGet();
    *cpl_ranks = vect_cpl_ranks->data();

    return vect_cpl_ranks->size();
  }


  /*----------------------------------------------------------------------------*
   * Methods about exchange frequency                                           *
   *----------------------------------------------------------------------------*/




  /*----------------------------------------------------------------------------*
   * Methods about spatial interpolation                                        *
   *----------------------------------------------------------------------------*/


 /**
  * \brief Set a spatial interpolation property of type double.
  *
  * \param [in]   name       Name of the property
  * \param [in]   value      Value of the property
  *
  */
 void
  Coupling::spatialInterpPropertyDoubleSet (
    std::string name,
    double      value
  )
{
  std::map<std::string, double>::iterator it;
  it = _spatial_interp_properties_double.find(name);

  if (it == _spatial_interp_properties_double.end()) {
    pair<string, double> newPair(name, value);
    _spatial_interp_properties_double.insert(newPair);
  } else {
    it->second = value;
  }
}


/**
  * \brief Set a spatial interpolation property of type int.
  *
  * \param [in]   name       Name of the property
  * \param [in]   value      Value of the property
  *
  */
 void
  Coupling::spatialInterpPropertyIntSet (
    std::string name,
    int         value
  )
  {
    std::map<std::string, int>::iterator it;
    it = _spatial_interp_properties_int.find(name);

    if (it == _spatial_interp_properties_int.end()) {
      pair<string, int> newPair(name, value);
      _spatial_interp_properties_int.insert(newPair);
    } else {
      it->second = value;
    }
  }


  /**
   *
   * \brief Export mesh to Ensight format
   *
   */

  void 
  Coupling::exportMesh(Coupling &cpl)
  {

    printf("Coupling::exportMesh\n");


    if (cpl._writer != NULL) {

      printf("exportMesh : %s\n", cpl._localCodeProperties.nameGet().c_str());
      fflush(stdout);

      /* First, create geometry and variables if necessary */
      if (cpl._n_step == 0) {

        cpl._id_geom_writer =PDM_writer_geom_create_from_mesh_nodal (cpl._writer,
                                                                     "geom",
                                                                     cpl._mesh.getPdmNodalIndex());

        PDM_writer_var_loc_t PDMfieldType = PDM_WRITER_VAR_ELEMENTS;

        PDM_writer_var_dim_t PDMfieldComp = PDM_WRITER_VAR_SCALAR;

        PDM_writer_status_t  st_dep_tps = PDM_WRITER_ON;

        if (cpl._displacement == CWP_DYNAMIC_MESH_STATIC) {
          st_dep_tps = PDM_WRITER_OFF;
        }

        cpl._id_field_partitioning_writer = PDM_writer_var_create(cpl._writer,
                                                                  st_dep_tps,
                                                                  PDMfieldComp,
                                                                  PDMfieldType,
                                                                  "partitioning");

        cpl._id_field_ranking_writer = PDM_writer_var_create(cpl._writer,
                                                             st_dep_tps,
                                                             PDMfieldComp,
                                                             PDMfieldType,
                                                             "ranking");
      }

      /* Then begin a new time step if none is currently open */
      if (!PDM_writer_is_open_step (cpl._writer) && ((cpl._n_step % cpl._freq_writer) == 0)) {

        double current_time;

        cpl._localCodeProperties.ctrlParamGet("time", &current_time);

        PDM_writer_step_beg (cpl._writer, current_time);

      }


      /* Finally write geometry and variables */
      if (cpl._n_step == 0) {

        PDM_writer_geom_write(cpl._writer, cpl._id_geom_writer);

        std::vector <double *> partitioning_field_data(cpl._mesh.getNPart());
        std::vector <double *> ranking_field_data(cpl._mesh.getNPart());

        int worldRank;
        MPI_Comm_rank (cpl._localCodeProperties.connectableCommGet(), &worldRank);

        int n_part = cpl._mesh.getNPart();
        int g_n_part = 0;

        MPI_Scan (&n_part, &g_n_part, 1, MPI_INT, MPI_SUM, cpl._localCodeProperties.connectableCommGet());

        g_n_part += -n_part;        

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          partitioning_field_data[i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );
          ranking_field_data     [i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );

          for(int i_elt = 0; i_elt < cpl._mesh.getPartNElts(i_part); i_elt++){
            partitioning_field_data[i_part][i_elt] = (double) (i_part + g_n_part);
            ranking_field_data[i_part][i_elt] = (double) worldRank;
          }

          PDM_writer_var_set(cpl._writer, cpl._id_field_partitioning_writer, cpl._id_geom_writer, i_part, (double *) partitioning_field_data[i_part]);
          PDM_writer_var_set(cpl._writer, cpl._id_field_ranking_writer     , cpl._id_geom_writer, i_part, (double *) ranking_field_data     [i_part]);
        }

        PDM_writer_var_write(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_write(cpl._writer, cpl._id_field_ranking_writer);

        PDM_writer_var_data_free(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_data_free(cpl._writer, cpl._id_field_ranking_writer);

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          free (partitioning_field_data[i_part]);
          free (ranking_field_data     [i_part]);
          partitioning_field_data[i_part] = nullptr;
          ranking_field_data     [i_part] = nullptr;
        }
      }

      else if (((cpl._n_step % cpl._freq_writer) == 0) && (cpl._displacement != CWP_DYNAMIC_MESH_STATIC)) {

        if (cpl._displacement == CWP_DYNAMIC_MESH_VARIABLE) {
        
          PDM_writer_geom_set_from_mesh_nodal (cpl._writer, 
                                               cpl._id_geom_writer,
                                               cpl._mesh.getPdmNodalIndex());
        }

        PDM_writer_geom_write(cpl._writer, cpl._id_geom_writer);

        std::vector <double *> partitioning_field_data(cpl._mesh.getNPart());
        std::vector <double *> ranking_field_data(cpl._mesh.getNPart());

        int worldRank;
        MPI_Comm_rank (cpl._localCodeProperties.connectableCommGet(), &worldRank);

        int n_part = cpl._mesh.getNPart();
        int g_n_part = 0;

        MPI_Scan (&n_part, &g_n_part, 1, MPI_INT, MPI_SUM, cpl._localCodeProperties.connectableCommGet());

        g_n_part += -n_part;        

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          partitioning_field_data[i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );
          ranking_field_data     [i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );

          for(int i_elt = 0; i_elt < cpl._mesh.getPartNElts(i_part); i_elt++){
            partitioning_field_data[i_part][i_elt] = (double) (i_part + g_n_part);
            ranking_field_data[i_part][i_elt] = (double) worldRank;
          }

          PDM_writer_var_set(cpl._writer, cpl._id_field_partitioning_writer, cpl._id_geom_writer, i_part, (double *) partitioning_field_data[i_part]);
          PDM_writer_var_set(cpl._writer, cpl._id_field_ranking_writer     , cpl._id_geom_writer, i_part, (double *) ranking_field_data     [i_part]);
        }

        PDM_writer_var_write(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_write(cpl._writer, cpl._id_field_ranking_writer);

        PDM_writer_var_data_free(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_data_free(cpl._writer, cpl._id_field_ranking_writer);

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          free (partitioning_field_data[i_part]);
          free (ranking_field_data     [i_part]);
          partitioning_field_data[i_part] = nullptr;
          ranking_field_data     [i_part] = nullptr;
        }
      }
    }

    else {
      printf(" sortie NULL\n");
    }

    fflush(stdout);
  }


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
    // Export mesh and associted fields                                        //
    //                                                                         //
    /////////////////////////////////////////////////////////////////////////////



    if (!_coupledCodeProperties.localCodeIs()) {

      exportMesh(*this);

    }

    else {
      int codeID    = localCodePropertiesGet()->idGet();
      int cplCodeID = coupledCodePropertiesGet()->idGet();

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        if (codeID < cplCodeID) {

          exportMesh(*this);
          exportMesh(cpl_cpl);   

        }
      }
    }

    /////////////////////////////////////////////////////////////////////////////
    //                                                                         //
    // Exchange fields properties to obtain spatial intepolation to build      //
    //                                                                         //
    /////////////////////////////////////////////////////////////////////////////

    // - Store Data to send

    if (!_coupledCodeProperties.localCodeIs()) {

      int codeID    = localCodePropertiesGet()->idGet();
      int cplCodeID = coupledCodePropertiesGet()->idGet();

      if (_n_step == 0) {

        std::string localFieldsName="";
        vector<int> localFieldsNameIdx;

        int localNbField = (int) _fields.size();

        localFieldsNameIdx.reserve(localNbField + 1);
        localFieldsNameIdx.push_back(0);

        std::vector<CWP_Field_exch_t> localFieldsExch;
        localFieldsExch.reserve(localNbField);

        std::vector<CWP_Dof_location_t> localFieldLocationV;
        localFieldLocationV.reserve(localNbField);

        std::map <std::string, cwipi::Field *>::iterator it = _fields.begin();

        while(it != _fields.end()){
          cwipi::Field* field = it->second;

          localFieldsName += it->first;
          localFieldsNameIdx.push_back((int) (localFieldsNameIdx[localFieldsNameIdx.size()-1]+it->first.size()));
          localFieldLocationV.push_back(field->locationGet());
          localFieldsExch.push_back(field->exchangeTypeGet());

          it++;
        }

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   1,
                                                                   (void *) &_nPart,
                                                                   -1,
                                                                   NULL,
                                                                   1,
                                                                   (void *) &_cplNPart,
                                                                   -1,
                                                                   NULL);
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

        nSendData   = localNbField + 1;
        nRecvData   = cplNbField + 1;
        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   nSendData,
                                                                   (void *) &(localFieldsNameIdx[0]),
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) &(cplFieldNameIdx[0]),
                                                                   -1,
                                                                   NULL);

        cplFieldName.resize(cplFieldNameIdx[cplNbField]);
        nSendData   = localFieldsNameIdx[localNbField];
        nRecvData   = cplFieldNameIdx[cplNbField];

        localFieldsNameIdx.clear();


        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(char),
                                                                   nSendData,
                                                                   (void *) localFieldsName.c_str(),
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) cplFieldName.c_str(),
                                                                   -1,
                                                                   NULL);

        localFieldsName.clear();

        nSendData   = localNbField;
        nRecvData   = cplNbField;
        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Field_exch_t),
                                                                   nSendData,
                                                                   (void *) &(localFieldsExch[0]),
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) &(cplFieldExch[0]),
                                                                   -1,
                                                                   NULL);
        localFieldsExch.clear();

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

        it = _fields.begin();
        while(it != _fields.end()) {
          cwipi::Field* field = it->second;
          string localFieldName = it->first;
          CWP_Field_exch_t   localFieldExch     = field->exchangeTypeGet();
          CWP_Dof_location_t localFieldLocation = field->locationGet();

          for (int j = 0; j < cplNbField; j++) {
            string _cplFieldName = cplFieldName.substr(cplFieldNameIdx[j], cplFieldNameIdx[j+1]-cplFieldNameIdx[j]);

            if (_cplFieldName == localFieldName) {
              field->linkedFieldLocationSet(cplFieldLocationV[j]);

              if (  (localFieldExch == CWP_FIELD_EXCH_SENDRECV && cplFieldExch[j] == CWP_FIELD_EXCH_SENDRECV)
                  ||(localFieldExch == CWP_FIELD_EXCH_SEND     && cplFieldExch[j] == CWP_FIELD_EXCH_RECV)
                  ||(localFieldExch == CWP_FIELD_EXCH_RECV     && cplFieldExch[j] == CWP_FIELD_EXCH_SEND)) {

                if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_SEND) {
                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]);
                  if (_spatial_interp_send.find(newKey) == _spatial_interp_send.end()) {

                    _spatial_interp_send.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));

                    _spatial_interp_send[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_SEND);
                  }
                }

                if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_RECV) {
                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]);
                  if (_spatial_interp_recv.find(newKey) == _spatial_interp_recv.end()) {

                    _spatial_interp_recv.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));

                    _spatial_interp_recv[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_RECV);

                  }
                }
              }
            }
          }
          it++;
        }

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

        _sis_loc_r.resize(2*sir_s);

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                   2*sis_s,
                                                                   (void *) &(sis_loc[0]),
                                                                   -1,
                                                                   NULL,
                                                                   2*sis_r,
                                                                   (void *) &(_sis_loc_r[0]),
                                                                   -1,
                                                                   NULL);

        vector<CWP_Dof_location_t> sir_loc_r;
        sir_loc_r.resize(2*sis_s);

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                   2*sir_s,
                                                                   (void *) &(sir_loc[0]),
                                                                   -1,
                                                                   NULL,
                                                                   2*sir_r,
                                                                   (void *) &(sir_loc_r[0]),
                                                                   -1,
                                                                   NULL);
      }

      if (codeID < cplCodeID) {

        // spatial_interp send

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          sis_it->second->weightsCompute();
          sis_it++;
        }

        // spatial_interp recv

        for (int i = 0; i < (int) (_sis_loc_r.size()/2); i++) {
          _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->weightsCompute();
        }

      }

      else {

        // spatial_interp recv

        for (int i = 0; i < (int) (_sis_loc_r.size()/2); i++) {
          _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->weightsCompute();
        }

        // spatial_interp send

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          sis_it->second->weightsCompute();
          sis_it++;
        }

      }

    }

    // if an instance of each code on the same rank. All work is made in the call from the smallest code id

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        int codeID    = localCodePropertiesGet()->idGet();
        int cplCodeID = coupledCodePropertiesGet()->idGet();

        if (_n_step == 0) {

          std::string localFieldsName="";
          vector<int> localFieldsNameIdx;

          int localNbField = (int) _fields.size();

          localFieldsNameIdx.reserve(localNbField + 1);
          localFieldsNameIdx.push_back(0);

          std::vector<CWP_Field_exch_t> localFieldsExch;
          localFieldsExch.reserve(localNbField);

          std::vector<CWP_Dof_location_t> localFieldLocationV;
          localFieldLocationV.reserve(localNbField);

          std::map <std::string, cwipi::Field *>::iterator it = _fields.begin();

          while(it != _fields.end()){
            cwipi::Field* field = it->second;

            localFieldsName += it->first;
            localFieldsNameIdx.push_back((int) (localFieldsNameIdx[localFieldsNameIdx.size()-1]+it->first.size()));
            localFieldLocationV.push_back(field->locationGet());
            localFieldsExch.push_back(field->exchangeTypeGet());

            it++;
          }

          std::string cpl_localFieldsName="";
          vector<int> cpl_localFieldsNameIdx;

          int cpl_localNbField = (int) cpl_cpl._fields.size();

          cpl_localFieldsNameIdx.reserve(cpl_localNbField + 1);
          cpl_localFieldsNameIdx.push_back(0);

          std::vector<CWP_Field_exch_t> cpl_localFieldsExch;
          cpl_localFieldsExch.reserve(cpl_localNbField);

          std::vector<CWP_Dof_location_t> cpl_localFieldLocationV;
          cpl_localFieldLocationV.reserve(cpl_localNbField);

          it = cpl_cpl._fields.begin();

          while(it != cpl_cpl._fields.end()){
            cwipi::Field* field = it->second;

            cpl_localFieldsName += it->first;
            cpl_localFieldsNameIdx.push_back((int) (cpl_localFieldsNameIdx[cpl_localFieldsNameIdx.size()-1]+it->first.size()));
            cpl_localFieldLocationV.push_back(field->locationGet());
            cpl_localFieldsExch.push_back(field->exchangeTypeGet());

            it++;
          }

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     1,
                                                                     (void *) &_nPart,
                                                                     1,
                                                                     (void *) &cpl_cpl._nPart,
                                                                     1,
                                                                     (void *) &_cplNPart,
                                                                     1,
                                                                     (void *) &cpl_cpl._cplNPart);

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

          vector<int               > cplFieldsNameIdx (cplNbField + 1, 0);
          vector<CWP_Field_exch_t  > cplFieldsExch (cplNbField);
          vector<CWP_Dof_location_t> cplFieldsLocationV (cplNbField);
          string                     cplFieldsName;

          vector<int               > cpl_cplFieldsNameIdx (cpl_cplNbField + 1, 0);
          vector<CWP_Field_exch_t  > cpl_cplFieldsExch (cpl_cplNbField);
          vector<CWP_Dof_location_t> cpl_cplFieldsLocationV (cpl_cplNbField);
          string                     cpl_cplFieldsName;

          // - Transfer memory to receive data

          nSendData   = localNbField + 1;
          nRecvData   = cplNbField + 1;

          cpl_nSendData   = cpl_localNbField + 1;
          cpl_nRecvData   = cpl_cplNbField + 1;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     nSendData,
                                                                     (void *) &(localFieldsNameIdx[0]),
                                                                     cpl_nSendData,
                                                                     (void *) &(cpl_localFieldsNameIdx[0]),
                                                                     nRecvData,
                                                                     (void *) &(cplFieldsNameIdx[0]),
                                                                     cpl_nRecvData,
                                                                     (void *) &(cpl_cplFieldsNameIdx[0]));

          cplFieldsName.resize(cplFieldsNameIdx[cplNbField]);
          nSendData   = localFieldsNameIdx[localNbField];
          nRecvData   = cplFieldsNameIdx[cplNbField];

          cpl_cplFieldsName.resize(cpl_cplFieldsNameIdx[cpl_cplNbField]);
          cpl_nSendData   = cpl_localFieldsNameIdx[cpl_localNbField];
          cpl_nRecvData   = cpl_cplFieldsNameIdx[cpl_cplNbField];

          localFieldsNameIdx.clear();
          cpl_localFieldsNameIdx.clear();


          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(char),
                                                                     nSendData,
                                                                     (void *) localFieldsName.c_str(),
                                                                     cpl_nSendData,
                                                                     (void *) cpl_localFieldsName.c_str(),
                                                                     nRecvData,
                                                                     (void *) cplFieldsName.c_str(),
                                                                     cpl_nRecvData,
                                                                     (void *) cpl_cplFieldsName.c_str());

          localFieldsName.clear();
          cpl_localFieldsName.clear();

          nSendData   = localNbField;
          nRecvData   = cplNbField;

          cpl_nSendData   = cpl_localNbField;
          cpl_nRecvData   = cpl_cplNbField;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Field_exch_t),
                                                                     nSendData,
                                                                     (void *) &(localFieldsExch[0]),
                                                                     cpl_nSendData,
                                                                     (void *) &(cpl_localFieldsExch[0]),
                                                                     nRecvData,
                                                                     (void *) &(cplFieldsExch[0]),
                                                                     cpl_nRecvData,
                                                                     (void *) &(cpl_cplFieldsExch[0]));

          localFieldsExch.clear();
          cpl_localFieldsExch.clear();

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
                                                                     (void *) &(cplFieldsLocationV[0]),
                                                                     cpl_nRecvData,
                                                                     (void *) &(cpl_cplFieldsLocationV[0]));

          localFieldLocationV.clear();
          cpl_localFieldLocationV.clear();

          // Create spatial interpolation objects

          it = _fields.begin();
          while(it != _fields.end()) {
            cwipi::Field* field = it->second;
            string localFieldName = it->first;
            CWP_Field_exch_t   localFieldExch     = field->exchangeTypeGet();
            CWP_Dof_location_t localFieldLocation = field->locationGet();

            for (int j = 0; j < cplNbField; j++) {

              string cplFieldName = cplFieldsName.substr( cplFieldsNameIdx[j], cplFieldsNameIdx[j+1]-cplFieldsNameIdx[j] );

              if (cplFieldName == localFieldName) {
                field->linkedFieldLocationSet(cplFieldsLocationV[j]);
                if (cpl_cpl._fields.find(cplFieldName) == cpl_cpl._fields.end()) {
                  PDM_error(__FILE__, __LINE__, 0, "'%s' Field not found\n",
                            cplFieldName.c_str());
                }

                cpl_cpl._fields[cplFieldName]->linkedFieldLocationSet(localFieldLocation);

                if (  (localFieldExch == CWP_FIELD_EXCH_SENDRECV && cplFieldsExch[j] == CWP_FIELD_EXCH_SENDRECV)
                    ||(localFieldExch == CWP_FIELD_EXCH_SEND     && cplFieldsExch[j] == CWP_FIELD_EXCH_RECV)
                    ||(localFieldExch == CWP_FIELD_EXCH_RECV     && cplFieldsExch[j] == CWP_FIELD_EXCH_SEND)) {

                  if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_SEND) {
                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldsLocationV[j]);
                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldsLocationV[j], localFieldLocation);

                    if (_spatial_interp_send.find(newKey) == _spatial_interp_send.end()) {

                      _spatial_interp_send.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                      cpl_cpl._spatial_interp_recv.insert(make_pair(cpl_newKey, FG::getInstance().CreateObject(cpl_cpl._spatialInterpAlgo)));
              
                      _spatial_interp_send[newKey]->init(this, localFieldLocation, cplFieldsLocationV[j], SPATIAL_INTERP_EXCH_SEND);
                      cpl_cpl._spatial_interp_recv[cpl_newKey]->init(&cpl_cpl, cplFieldsLocationV[j], localFieldLocation, SPATIAL_INTERP_EXCH_RECV);        

                    }
                  }

                  if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_RECV) {

                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldsLocationV[j]);
                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldsLocationV[j], localFieldLocation);

                    if (_spatial_interp_recv.find(newKey) == _spatial_interp_recv.end()) {

                      _spatial_interp_recv.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                      cpl_cpl._spatial_interp_send.insert(make_pair(cpl_newKey, FG::getInstance().CreateObject(cpl_cpl._spatialInterpAlgo)));

                      _spatial_interp_recv[newKey]->init(this, localFieldLocation, cplFieldsLocationV[j], SPATIAL_INTERP_EXCH_RECV);
                      cpl_cpl._spatial_interp_send[cpl_newKey]->init(&cpl_cpl, cplFieldsLocationV[j], localFieldLocation, SPATIAL_INTERP_EXCH_SEND);

                    }
                  }
                }
              }
            }
            it++;
          }

          int sis_s = (int) _spatial_interp_send.size();

          vector<CWP_Dof_location_t> sis_loc;
          sis_loc.reserve(2*sis_s);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
          while(sis_it != _spatial_interp_send.end()) {
            sis_loc.push_back((sis_it->first).first);
            sis_loc.push_back((sis_it->first).second);
            sis_it++;
          }

          int cpl_sis_s = (int) cpl_cpl._spatial_interp_send.size();

          vector<CWP_Dof_location_t> cpl_sis_loc;
          cpl_sis_loc.reserve(2*cpl_sis_s);

          sis_it = cpl_cpl._spatial_interp_send.begin();
          while(sis_it != cpl_cpl._spatial_interp_send.end()) {
            cpl_sis_loc.push_back((sis_it->first).first);
            cpl_sis_loc.push_back((sis_it->first).second);
            sis_it++;
          }

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

          int sir_s = (int) _spatial_interp_recv.size();
    
          vector<CWP_Dof_location_t> sir_loc;
          sir_loc.reserve(2*sis_r);
    
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sir_it = _spatial_interp_recv.begin();
          while(sir_it != _spatial_interp_recv.end()) {
            sir_loc.push_back((sir_it->first).first);
            sir_loc.push_back((sir_it->first).second);
            sir_it++;
          }
    
          int cpl_sir_s = (int) cpl_cpl._spatial_interp_recv.size();
    
          vector<CWP_Dof_location_t> cpl_sir_loc;
          cpl_sir_loc.reserve(2*cpl_sir_s);
    
          sir_it = cpl_cpl._spatial_interp_recv.begin();
          while(sir_it != cpl_cpl._spatial_interp_recv.end()) {
            cpl_sir_loc.push_back((sir_it->first).first);
            cpl_sir_loc.push_back((sir_it->first).second);
            sir_it++;
          }
    
          int sir_r;
          int cpl_sir_r;

          // sir_s = 5;
          // cpl_sir_s = 7;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     1,
                                                                     (void *) &sir_s,
                                                                     1,
                                                                     (void *) &cpl_sir_s,
                                                                     1,
                                                                     (void *) &sir_r,
                                                                     1,
                                                                     (void *) &cpl_sir_r);


          _sis_loc_r.resize(2*sir_s);

          _cpl_sis_loc_r.resize(2*cpl_sir_s);
    
          assert(sir_r     == sis_s);
          assert(sis_r     == sir_s);
          assert(cpl_sir_r == cpl_sis_s);
          assert(cpl_sis_r == cpl_sir_s);

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                     2*sis_s,
                                                                     (void *) &(sis_loc[0]),
                                                                     2*cpl_sis_s,
                                                                     (void *) &(cpl_sis_loc[0]),
                                                                     2*sis_r,
                                                                     (void *) &(_sis_loc_r[0]),
                                                                     2*cpl_sis_r,
                                                                     (void *) &(_cpl_sis_loc_r[0]));
    
          vector<CWP_Dof_location_t> sir_loc_r;
          sir_loc_r.resize(2*sis_s);

          vector<CWP_Dof_location_t> cpl_sir_loc_r;
          cpl_sir_loc_r.resize(2*cpl_sis_s);
    
          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                     2*sir_s,
                                                                     (void *) &(sir_loc[0]),
                                                                     2*cpl_sir_s,
                                                                     (void *) &(cpl_sir_loc[0]),
                                                                     2*sir_r,
                                                                     (void *) &(sir_loc_r[0]),
                                                                     2*cpl_sir_r,
                                                                     (void *) &(cpl_sir_loc_r[0]));
        }

        if (codeID < cplCodeID) {
  
          int i = 0;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
          while(sis_it != _spatial_interp_send.end()) {
            sis_it->second->weightsCompute();
            cpl_cpl._spatial_interp_recv[make_pair(_cpl_sis_loc_r[2 * i + 1], _cpl_sis_loc_r[2 * i])]->weightsCompute();
            sis_it++;
            i++;
          }
  
          i = 0;
          sis_it = cpl_cpl._spatial_interp_send.begin();
          while(sis_it != cpl_cpl._spatial_interp_send.end()) {
            _spatial_interp_recv[make_pair(_sis_loc_r[2 * i + 1], _sis_loc_r[2 * i])]->weightsCompute();
            sis_it->second->weightsCompute();
            sis_it++;
            i++;
          }
        }

        else {
  
          int i = 0;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = cpl_cpl._spatial_interp_send.begin();
          while(sis_it != cpl_cpl._spatial_interp_send.end()) {
            _spatial_interp_recv[make_pair(_sis_loc_r[2 * i + 1], _sis_loc_r[2 * i])]->weightsCompute();
            sis_it->second->weightsCompute();
            sis_it++;
            i++;
          }
  
          i = 0;
          sis_it = _spatial_interp_send.begin();
          while(sis_it != _spatial_interp_send.end()) {
            sis_it->second->weightsCompute();
            cpl_cpl._spatial_interp_recv[make_pair(_cpl_sis_loc_r[2 * i + 1], _cpl_sis_loc_r[2 * i])]->weightsCompute();
            sis_it++;
            i++;
          }
        }
      }
    }

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
    CWP_UNUSED(format);

    if (_is_first_field_created) {
      PDM_error(__FILE__, __LINE__, 0,
                "Error : CWP_Visu_set has to be called before CWP_Field_create.\n");
      abort();
    }

    if (_is_mesh_finalized) {
      PDM_error(__FILE__, __LINE__, 0,
                "Error : CWP_Visu_set has to be called before CWP_Mesh_interf_finalize.\n");
      abort();
    }

    string CodeName = _localCodeProperties.nameGet();
    string cplCodeName = _coupledCodeProperties.nameGet();
    string cplId = IdGet();

    string visuDir = "cwipi";
    string output_dir = visuDir+"/"+cplId+"_"+CodeName+"_"+cplCodeName;
    string output_dir_tmp = visuDir+"_writer/"+cplId+"_"+CodeName+"_"+cplCodeName;

    int rank;

    MPI_Comm_rank(_communication.unionCommGet(),&rank);

    if ((commTypeGet() == CWP_COMM_PAR_WITH_PART ||
        (commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && rank == _communication.unionCommLocCodeRootRanksGet())) && freq > 0) {

      _freq_writer = freq;

      PDM_MPI_Comm pdmComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_localCodeProperties.connectableCommGet()));

      PDM_writer_topology_t pdm_topology = PDM_WRITER_TOPO_CST;

      if     ( _displacement == CWP_DYNAMIC_MESH_STATIC     ) {
        pdm_topology  = PDM_WRITER_TOPO_CST;
      }
      
      else if( _displacement == CWP_DYNAMIC_MESH_DEFORMABLE ) {
        pdm_topology  = PDM_WRITER_TOPO_DEFORMABLE;
      }

      else { 
        pdm_topology  = PDM_WRITER_TOPO_VARIABLE;
      }

      PDM_writer_fmt_fic_t fmt_fic      = PDM_WRITER_FMT_BIN;
      const char* fmt                   = "Ensight";
      PDM_writer_status_t st_reprise    = PDM_WRITER_OFF;
      const char *options_comp          = "";
      //Proportion of working node for file acess
      int working_node = 1;
  
      PDM_io_kind_t acess_type =   PDM_IO_KIND_MPI_SIMPLE;

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

        if (option == "text") {
          fmt_fic = PDM_WRITER_FMT_ASCII;
        }
        else if (option == "binary") {
          fmt_fic = PDM_WRITER_FMT_BIN;
        }
        else if (option != "") {
          PDM_error(__FILE__, __LINE__, 0,
                    "Not a valid visualization option.\n");
        } 

      }
      while (pos != std::string::npos);

      _writer = PDM_writer_create(fmt,
                                  fmt_fic,
                                  pdm_topology,
                                  st_reprise,
                                  (char *) output_dir_tmp.c_str(),
                                  (char *) string("chr").c_str(),
                                  pdmComm,
                                  acess_type,
                                  working_node,
                                  options_comp);
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

    _is_first_field_created = 1;

    if (fieldIs(field_id)) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' existing field\n", field_id.c_str());
    }

    //
    // Create the new field

    cwipi::Field *newField = new cwipi::Field(field_id,
                                              _fields.size(),
                                              data_type,
                                              this,
                                              fieldType,
                                              storage,
                                              n_component,
                                              exch_type,
                                              visu_status);  

    pair<string, Field* > newPair(string(field_id), newField);
    string localName = _localCodeProperties.nameGet();
    _fields.insert(newPair);
    // if (_visu.isCreated() && newField->visuStatusGet() == CWP_STATUS_ON) {
    //   _visu.WriterFieldCreate(newField);
    // }
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
      // if (_visu.isCreated() && It->second->visuStatusGet() == CWP_STATUS_ON) {
      //   _visu.fieldDataSet(It->second,map_type, i_part);
      // }
    }
  }

  /**
  *
  * \brief Get Field data
  *
  * \param [in]   field_id       Field identifier
  * \param [out]  data           Storage array (mapping)
  *
  */

  void
  Coupling::fieldDataGet
  (
    const string &field_id,
    int i_part,
    const CWP_Field_map_t   map_type,
    void** data
  )
  {

    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      *data = It->second->dataGet(i_part, map_type);
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
   * \param [in]  src_field_id              Source field (NULL->no sending)
   * \param [in]  tgt_field_id              Target field (NULL->no receiving)
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

        PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm ((void *) &(_localCodeProperties.intraCommGet()));

        PDM_gen_gnum_t *pgg  = PDM_gnum_create (3, _nPart, PDM_FALSE, 1e-3, comm,
                                                   PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);
        for (int iPart = 0; iPart <_nPart; iPart++){
          PDM_gnum_set_from_coords (pgg, iPart, _userTargetN[iPart], _userTargetCoord[iPart], NULL);
        }

        PDM_gnum_compute (pgg);

        _localUserTargetGnum = (CWP_g_num_t **) malloc (sizeof(CWP_g_num_t *) * _nPart);

        for (int iPart = 0; iPart < _nPart; iPart++){
          _localUserTargetGnum[iPart] = const_cast <CWP_g_num_t*> (PDM_gnum_get (pgg, iPart));
        }

        PDM_gnum_free (pgg);

        _userTargetGnum = const_cast <const CWP_g_num_t**> (_localUserTargetGnum);
      }
    }
  }

  /**
   * \brief Update time.
   *
   * \param [in]  current_time     Current time
   *
   */

  void
  Coupling::timeUpdate (
    double current_time
  )
  {

    _n_step++;

    if (_writer != NULL) {
      PDM_writer_step_end(_writer);
    }

    // if(_visu.isCreated() and _visu.physicalTimeGet() > -1.0) {
    //    _visu.WriterStepEnd();
    // }

    std::map < string, Field * >::iterator itf = _fields.begin();
    while (itf != _fields.end()) {
      itf->second->currentStepWasExchangedReset();
      itf++;
    }

    if (_writer != NULL && (_n_step % _freq_writer == 0)) {
      PDM_writer_step_beg(_writer, current_time);
    }

    // if(_visu.isCreated()) {
    //    _visu.WriterStepBegin(current_time, &_mesh);
    // }

  }

// A supprimer

  CWP_g_num_t*
  Coupling::globalNumGet(int id_block,int i_part)
  {
    return _mesh.globalNumGet(id_block,i_part);
  }

  // int
  // Coupling::isUpToDateGet ()
  // {
  //   return _is_up_to_date;
  // }

  // void
  // Coupling::isUpToDateSet ()
  // {
  //   _is_up_to_date = 1;
  // }



} // namespace cwipi

/**
 * \endcond
 */
