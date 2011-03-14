/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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

#include <algorithm>
#include <cmath>
#include <cassert>

#include <bft_printf.h>

#include "locationToLocalMesh.hxx"
#include "locationToDistantMesh.hxx"
#include "applicationProperties.hxx"

namespace cwipi
{

LocationToLocalMesh::LocationToLocalMesh(
                                         const cwipi_solver_type_t  &solverType,
                                         const double &tolerance,
                                         const MPI_Comm& couplingComm,
                                         const int &coupledApplicationNRankCouplingComm,
                                         const int &coupledApplicationBeginningRankCouplingComm,
                                         const bool isCoupledRank,
                                         const int entitiesDim,
                                         const ApplicationProperties& localApplicationProperties,
                                         LocationToDistantMesh &locationToDistantMesh)

: _solverType(solverType),  _tolerance(tolerance),
 _locationToDistantMesh(locationToDistantMesh), _couplingComm(couplingComm),
 _coupledApplicationNRankCouplingComm(coupledApplicationNRankCouplingComm),
 _coupledApplicationBeginningRankCouplingComm(coupledApplicationBeginningRankCouplingComm),
 _isCoupledRank(isCoupledRank), _entitiesDim(entitiesDim), _localApplicationProperties(localApplicationProperties)

{
  _fvmLocator = NULL;
  _barycentricCoordinatesIndex = NULL;
  _barycentricCoordinates = NULL;
  _nDistantPoint = 0;
  _location = NULL;
  _toLocate = true;
  _maxElementContainingNVertex = -1;
  _nVertex = NULL;
  _supportMesh = NULL;
}

LocationToLocalMesh::~LocationToLocalMesh()
{

  if (_barycentricCoordinatesIndex != NULL) {
    delete _barycentricCoordinatesIndex;
    delete _barycentricCoordinates;
  }

  if (_fvmLocator != NULL)
    fvm::fvm_locator_destroy(_fvmLocator);

  if (_nVertex != NULL)
    delete _nVertex;

}

void LocationToLocalMesh::locate()
{
  //
  // Create locator

  // TODO: Exchange MPI of _toLocate param between master rank !!!
  // TODO: Attention la fonction clear ne devrait-elle pas etre dans le if suivant !

  _locationToDistantMesh.clear();

  if ( _isCoupledRank && (_toLocate || _locationToDistantMesh.getToLocateStatus())) {

    if (_supportMesh == NULL)
      bft::bft_error(__FILE__, __LINE__, 0, "undefined support mesh\n");

    if (_fvmLocator == NULL)
      _fvmLocator = fvm::fvm_locator_create(_tolerance,
                                       _couplingComm,
                                       _coupledApplicationNRankCouplingComm,
                                       _coupledApplicationBeginningRankCouplingComm);

    // TODO: Revoir les coordonnees des points a localiser (cas centres sommets + centres faces + autres points)
    // TODO: Ajouter un locator pour les sommets pour les centres faces,...

    double* coords = NULL;
    if (_locationToDistantMesh._coordsPointsToLocate != NULL)
      coords = _locationToDistantMesh._coordsPointsToLocate;

    else if(_solverType == CWIPI_SOLVER_CELL_CENTER) {
      _locationToDistantMesh._nPointsToLocate = _supportMesh->getNElts();
      coords = const_cast <double*> (&(_supportMesh->getCellCenterCoords()[0]));
    }

    else if(_solverType == CWIPI_SOLVER_CELL_VERTEX) {
      _locationToDistantMesh._nPointsToLocate = _supportMesh->getNVertex();
      coords = const_cast <double*> (_supportMesh->getVertexCoords());
    }

    fvm::fvm_locator_set_nodal(_fvmLocator,
                               &_supportMesh->getFvmNodal(),
                               0,
                               3,
                               _locationToDistantMesh._nPointsToLocate,
                               NULL,
                               coords);

    _toLocate = false;
    const int nLocatedPoint = fvm::fvm_locator_get_n_interior(_fvmLocator);
    const int nNotLocatedPoint = _locationToDistantMesh._nPointsToLocate - nLocatedPoint;
    const int* exteriorList = fvm::fvm_locator_get_exterior_list(_fvmLocator);
    const int* interiorList = fvm::fvm_locator_get_interior_list(_fvmLocator);
    const int* locationList = fvm::fvm_locator_get_dist_locations(_fvmLocator);
    const int nExterior = fvm::fvm_locator_get_n_exterior(_fvmLocator);
    assert(nNotLocatedPoint == nExterior);

    _locationToDistantMesh._unlocatedPoint = const_cast<int *> (exteriorList);
    _locationToDistantMesh._locatedPoint = const_cast<int *> (interiorList);
    _locationToDistantMesh._nUnlocatedPoint = nNotLocatedPoint;
    _locationToDistantMesh._nLocatedPoint = nLocatedPoint;
    _locationToDistantMesh._toLocate = false;

    _location = const_cast<int *> (locationList);
    _nDistantPoint = fvm::fvm_locator_get_n_dist_points(_fvmLocator);

    if (_barycentricCoordinatesIndex != NULL) {
      delete _barycentricCoordinatesIndex;
      delete _barycentricCoordinates;
    }

    //
    // TODO: Prevoir une fabrique pour supprimer les tests if sur _entitiesDim
    //       Le calcul des coordonnees barycentriques se fera dans cette fabrique

    if (_barycentricCoordinatesIndex == NULL) {
      if (_entitiesDim == 1) {
        const int nDistantPoint      = fvm::fvm_locator_get_n_dist_points(_fvmLocator);
        const int *distantLocation   = fvm::fvm_locator_get_dist_locations(_fvmLocator);
        const double *distantCoords   = fvm::fvm_locator_get_dist_coords(_fvmLocator);

        const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
        const int *eltsConnec = _supportMesh->getEltConnectivity();
        const double *localCoords    = _supportMesh->getVertexCoords();

        if ( nDistantPoint > 0 ) {
          _barycentricCoordinatesIndex = new std::vector <int> (nDistantPoint + 1);
          _barycentricCoordinates = new std::vector <double> (2 * nDistantPoint);
          std::vector <int> &  _refBarycentricCoordinatesIndex = *_barycentricCoordinatesIndex;
          std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;

          _refBarycentricCoordinatesIndex[0] = 0;

          for (int ipoint = 0; ipoint < nDistantPoint ; ipoint++) {
            int iel = distantLocation[ipoint] - 1;
            _refBarycentricCoordinatesIndex[ipoint+1] = _refBarycentricCoordinatesIndex[ipoint] + 2;
            int index = eltsConnecPointer[iel];
            int nVertex = eltsConnecPointer[iel+1] - eltsConnecPointer[iel];
            assert(nVertex == 2);
            int pt1 = eltsConnec[index] - 1;
            int pt2 = eltsConnec[index+1] - 1;
            double coef1 = sqrt((localCoords[3*pt1]-distantCoords[3*ipoint])*(localCoords[3*pt1]-distantCoords[3*ipoint])+
                                (localCoords[3*pt1+1]-distantCoords[3*ipoint+1])*(localCoords[3*pt1+1]-distantCoords[3*ipoint+1])+
                                (localCoords[3*pt1+2]-distantCoords[3*ipoint+2])*(localCoords[3*pt1+2]-distantCoords[3*ipoint+2]));
            double coef2 = sqrt((localCoords[3*pt2]-distantCoords[3*ipoint])*(localCoords[3*pt2]-distantCoords[3*ipoint])+
                                (localCoords[3*pt2+1]-distantCoords[3*ipoint+1])*(localCoords[3*pt2+1]-distantCoords[3*ipoint+1])+
                                (localCoords[3*pt2+2]-distantCoords[3*ipoint+2])*(localCoords[3*pt2+2]-distantCoords[3*ipoint+2]));
            _refBarycentricCoordinates[_refBarycentricCoordinatesIndex[ipoint]] = coef2/(coef1+coef2);
            _refBarycentricCoordinates[_refBarycentricCoordinatesIndex[ipoint]+1] = coef1/(coef1+coef2);
          }
        }
      }

      else if (_entitiesDim == 2) {
        compute2DMeanValues();
      }

      else if (_entitiesDim == 3) {
        //TODO: calcul des coord barycentriques 3D
        int nPoints;
        //         coo_baryc(_fvmLocator,
        //                   _supportMesh->getNVertex(),
        //                   _supportMesh->getVertexCoords(),
        //                   _supportMesh->getNElts(),
        //                   _supportMesh->getEltConnectivityIndex(),
        //                   _supportMesh->getEltConnectivity(),
        //                   &nPoints,
        //                   &_barycentricCoordinatesIndex,
        //                   &_barycentricCoordinates);
      }
    }

    //
    // Exchange info status between root ranks

    int currentRank;
    MPI_Comm_rank(_couplingComm, &currentRank);
    cwipi_located_point_info_t distantInfo;
    cwipi_located_point_info_t localInfo = _locationToDistantMesh._locationInfo;
    MPI_Status MPIStatus;

    const bool isRootRank = (currentRank == 0 ||
                            currentRank == _coupledApplicationBeginningRankCouplingComm +
                            _coupledApplicationNRankCouplingComm);

    int rootRank;
    if (isRootRank) {

      MPI_Sendrecv((int*) &localInfo,   1, MPI_INT,
                   _coupledApplicationBeginningRankCouplingComm, 0,
                   (int*) &distantInfo, 1, MPI_INT,
                   _coupledApplicationBeginningRankCouplingComm, 0,
                   _couplingComm, &MPIStatus);
    }

    //
    // application rootRank send the value of distantInfo to the other applications in the couplingComm

    int sizeComm;
    MPI_Comm_size(_couplingComm, &sizeComm);
    MPI_Status status;

    if ( _coupledApplicationBeginningRankCouplingComm == 0)
      rootRank = _coupledApplicationNRankCouplingComm;
    else
      rootRank = 0;

    if (! isRootRank)
      MPI_Recv(&distantInfo, 1, MPI_INT, rootRank, 0, _couplingComm, &status);
    else {
      int rank1;
      int rank2;

      if (rootRank == 0) {
        rank1 = 1;
        rank2 = _coupledApplicationBeginningRankCouplingComm;
      }
      else {
        rank1 = _coupledApplicationNRankCouplingComm+1;
        rank2 = sizeComm;
      }

      for (int i = rank1; i < rank2; i++)
        MPI_Send(&distantInfo, 1, MPI_INT, i, 0, _couplingComm);

    }

    //
    // Exchange info about distant mesh

    if (distantInfo == CWIPI_DISTANT_MESH_INFO ||
        localInfo == CWIPI_DISTANT_MESH_INFO) {

      //
      // Location

      if (localInfo == CWIPI_DISTANT_MESH_INFO)
        _locationToDistantMesh._elementContaining = new int[_locationToDistantMesh._nLocatedPoint];

      int *pLocation = NULL;
      if (distantInfo == CWIPI_DISTANT_MESH_INFO)
        pLocation = _location;

      fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                     (void *) pLocation,
                                     (void *) _locationToDistantMesh._elementContaining,
                                     NULL,
                                     sizeof(int),
                                     1,
                                     0);

      //
      // elementContainingMPIrankContaining

      if (localInfo == CWIPI_DISTANT_MESH_INFO)
        _locationToDistantMesh._elementContainingMPIrankContaining = new int[_locationToDistantMesh._nLocatedPoint];

      int *pMPIrank = NULL;
      std::vector <int> *MPIrank = NULL;

      if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
        int rankInGlobalComm;
        MPI_Comm_rank(_localApplicationProperties.getGlobalComm(), &rankInGlobalComm);
        MPIrank = new std::vector <int> (_nDistantPoint, rankInGlobalComm);
        pMPIrank = &((*MPIrank)[0]);
      }

      fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                     (void *) pMPIrank,
                                     (void *) _locationToDistantMesh._elementContainingMPIrankContaining,
                                     NULL,
                                     sizeof(int),
                                     1,
                                     0);

      if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
        delete MPIrank;
      }

      //
      // ElementContainingNVertex

      if (localInfo == CWIPI_DISTANT_MESH_INFO)
        _locationToDistantMesh._elementContainingNVertex = new int[_locationToDistantMesh._nLocatedPoint + 1];

      int *p_nVertex = NULL;
      _maxElementContainingNVertex = 0;

      if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
        _nVertex = new std::vector <int> (_nDistantPoint, 0);
        std::vector <int> & _nVertexRef = *_nVertex;

        for (int i = 0; i < _nDistantPoint; i++)
          _nVertexRef[i] = _supportMesh->getEltConnectivityIndex()[_location[i]] -
                                           _supportMesh->getEltConnectivityIndex()[_location[i]-1];

        p_nVertex = &(*_nVertex)[0];

        int localMaxElementContainingNVertex = *std::max_element(_nVertexRef.begin(), _nVertexRef.end());
        MPI_Allreduce (&localMaxElementContainingNVertex,
                       &_maxElementContainingNVertex,
                       1, MPI_INT, MPI_MAX,
                       _localApplicationProperties.getLocalComm());
      }

      int distantMaxElementContainingNVertex = 0;
      if (isRootRank) {
        MPI_Sendrecv((int*) &_maxElementContainingNVertex,   1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     (int*) &distantMaxElementContainingNVertex, 1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     _couplingComm, &MPIStatus);
      }

      MPI_Bcast(&distantMaxElementContainingNVertex, 1, MPI_INT, 0, _localApplicationProperties.getLocalComm());

      _maxElementContainingNVertex = std::max(_maxElementContainingNVertex, distantMaxElementContainingNVertex);

      fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                     (void *) p_nVertex,
                                     (void *) _locationToDistantMesh._elementContainingNVertex,
                                     NULL,
                                     sizeof(int),
                                     1,
                                     0);

      if (localInfo == CWIPI_DISTANT_MESH_INFO) {
        int previous1 = _locationToDistantMesh._elementContainingNVertex[0];
        _locationToDistantMesh._elementContainingNVertex[0] = 0;
        for (int i = 1; i < _locationToDistantMesh._nLocatedPoint + 1; i++) {
          int previous2 = _locationToDistantMesh._elementContainingNVertex[i];
          _locationToDistantMesh._elementContainingNVertex[i] = previous1 + _locationToDistantMesh._elementContainingNVertex[i-1];
          previous1 = previous2;
        }
      }

      //
      // ElementContainingVertex

      int *tmpLocal = NULL;
      int *tmpDistant = NULL;

      if (localInfo == CWIPI_DISTANT_MESH_INFO) {
        tmpDistant = new int [_locationToDistantMesh._nLocatedPoint];
        _locationToDistantMesh._elementContainingVertex = new int [_locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]];
      }

      if (distantInfo == CWIPI_DISTANT_MESH_INFO)
        tmpLocal = new int [_nDistantPoint];

      for (int i = 0; i < _maxElementContainingNVertex; i++) {
        if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
          std::vector <int> & _nVertexRef = *_nVertex;
          for (int j = 0; j < _nDistantPoint; j++) {
            if (i < _nVertexRef[j]) {
              int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
              tmpLocal[j] = _supportMesh->getEltConnectivity()[index];
            }
            else
              tmpLocal[j] = -1;
          }
        }

        fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                       (void *) tmpLocal,
                                       (void *) tmpDistant,
                                       NULL,
                                       sizeof(int),
                                       1,
                                       0);

        if (localInfo == CWIPI_DISTANT_MESH_INFO) {
          for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
            int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
            if (i < nVertices) {
              assert(tmpDistant[j] != -1);
              _locationToDistantMesh._elementContainingVertex[_locationToDistantMesh._elementContainingNVertex[j] + i] = tmpDistant[j];
            }
          }
        }
      }

      if (localInfo == CWIPI_DISTANT_MESH_INFO)
        delete [] tmpDistant;

      if (distantInfo == CWIPI_DISTANT_MESH_INFO)
        delete [] tmpLocal;

      //
      // ElementContainingBarycentricCoordinates

      double *tmpLocal1 = NULL;
      double *tmpDistant1 = NULL;

      if (localInfo == CWIPI_DISTANT_MESH_INFO) {
        tmpDistant1 = new double [_locationToDistantMesh._nLocatedPoint];
        _locationToDistantMesh._elementContainingBarycentricCoordinates = new double [_locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]];
      }

      if (distantInfo == CWIPI_DISTANT_MESH_INFO)
        tmpLocal1 = new double [_nDistantPoint];

      for (int i = 0; i < _maxElementContainingNVertex; i++) {
        std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;

        if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
          std::vector <int> & _nVertexRef = *_nVertex;
          for (int j = 0; j < _nDistantPoint; j++) {
            if (i < _nVertexRef[j]) {
              int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
              tmpLocal1[j] = _refBarycentricCoordinates[index];
            }
          }
        }

        fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                       (void *) tmpLocal1,
                                       (void *) tmpDistant1,
                                       NULL,
                                       sizeof(double),
                                       1,
                                       0);

        if (localInfo == CWIPI_DISTANT_MESH_INFO) {
          for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
            int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
            if (i < nVertices) {
              _locationToDistantMesh._elementContainingBarycentricCoordinates[_locationToDistantMesh._elementContainingNVertex[j] + i] = tmpDistant1[j];
            }
          }
        }
      }

      if (localInfo == CWIPI_DISTANT_MESH_INFO)
        delete [] tmpDistant1;

      if (distantInfo == CWIPI_DISTANT_MESH_INFO)
        delete [] tmpLocal1;

      // TODO: Optimisation a réaliser dans fvm en mettant un stride[dim]
      //       pour l'instant on calcule le max des nombre de sommets
      //       ou bien creer son propre graphe de communication MPI

      //
      // ElementContainingVertexCoords
      tmpLocal1 = NULL;
      tmpDistant1 = NULL;

      int stride = 3;
      if (localInfo == CWIPI_DISTANT_MESH_INFO) {
        tmpDistant1 = new double [stride * _locationToDistantMesh._nLocatedPoint];
        _locationToDistantMesh._elementContainingVertexCoords = new double [stride * _locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]];
      }

      if (distantInfo == CWIPI_DISTANT_MESH_INFO)
        tmpLocal1 = new double [stride * _nDistantPoint];

      for (int i = 0; i < _maxElementContainingNVertex; i++) {
        std::vector <int> & _nVertexRef = *_nVertex;
        if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
          for (int j = 0; j < _nDistantPoint; j++) {
            if (i < _nVertexRef[j]) {
              int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
              int icel = _supportMesh->getEltConnectivity()[index];
              for (int k = 0; k < stride; k++)
                tmpLocal1[stride *j + k] = _supportMesh->getVertexCoords()[stride * (icel - 1) + k];
            }
          }
        }

        fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                       (void *) tmpLocal1,
                                       (void *) tmpDistant1,
                                       NULL,
                                       sizeof(double),
                                       stride,
                                       0);

        if (localInfo == CWIPI_DISTANT_MESH_INFO) {
          for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
            int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
            if (i < nVertices) {
              for (int k = 0; k < stride; k++) {
                _locationToDistantMesh._elementContainingVertexCoords[stride * (_locationToDistantMesh._elementContainingNVertex[j] + i) + k] =
                  tmpDistant1[stride * j + k];
              }
            }
          }
        }
      }

      if (localInfo == CWIPI_DISTANT_MESH_INFO)
        delete [] tmpDistant1;

      if (distantInfo == CWIPI_DISTANT_MESH_INFO)
        delete [] tmpLocal1;
    }
  }

  // TODO: Attention la fonction synchronise ne devrait-elle pas etre dans le if precedent !

  _locationToDistantMesh.synchronize();

}

///
/// \brief Exchange field on vertices of cells that contain each located points
///

void LocationToLocalMesh::exchangeCellVertexFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride)
{

  // TODO: Pas optimal : une fois FVM amélioré ou remplacé : transmettre le champ  a une liste de sommets et plus sur la conectivité
  // TODO: exchangeCellVertexFieldOfElementContaining : doute pour un fonctionnement en parallele sans partitionnement
  double *tmpDistant1 = NULL;
  double *tmpLocal1 = NULL;

  if (receivingField != NULL) {
    tmpDistant1 = new double [stride * _locationToDistantMesh._nLocatedPoint];
    _locationToDistantMesh._elementContainingVertexCoords = new double [stride * _locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]];
  }

  if (sendingField != NULL)
    tmpLocal1 = new double [stride * _nDistantPoint];

  for (int i = 0; i < _maxElementContainingNVertex; i++) {
    std::vector <int> & _nVertexRef = *_nVertex;
    if (sendingField != NULL) {
      for (int j = 0; j < _nDistantPoint; j++) {
        if (i < _nVertexRef[j]) {
          int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
          int ivertex = _supportMesh->getEltConnectivity()[index];
          for (int k = 0; k < stride; k++)
            tmpLocal1[stride *j + k] = sendingField[stride*(ivertex - 1) + k];
        }
      }
    }

    fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                   (void *) tmpLocal1,
                                   (void *) tmpDistant1,
                                   NULL,
                                   sizeof(double),
                                   stride,
                                   0);

    if (_locationToDistantMesh._elementContainingNVertex != NULL) {
      for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
        int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
        if (i < nVertices) {
          for (int k = 0; k < stride; k++)
            receivingField[stride * (_locationToDistantMesh._elementContainingNVertex[j] + i) + k] = tmpDistant1[stride * j + k];
        }
      }
    }
  }

  delete [] tmpLocal1;
  delete [] tmpDistant1;
}

///
/// \brief Exchange field on cells that contain each located points
///

void LocationToLocalMesh::exchangeCellCenterFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride)
{

  // TODO: exchangeCellCenterFieldOfElementContaining : doute pour un fonctionnement en parallele sans partitionnement

  fvm::fvm_locator_exchange_point_var(_fvmLocator,
                                 (void *) sendingField,
                                 (void *) receivingField,
                                 NULL,
                                 sizeof(double),
                                 stride,
                                 0);
}

///
/// \brief Compute Mean Values
///
///

void  LocationToLocalMesh::midplaneProjection
(
 const int     nbr_som_fac,
 double *const coo_som_fac,
 double *const coo_point_dist
)

{

  int    icoo ;
  int    isom ;
  int    itri ;

  double   cost ;
  double   sint ;

  double   coo_tmp ;

  double   vect1[3] ;
  double   vect2[3] ;

  double  prod_vect[3] ;

  double  barycentre_fac[3] ;
  double  normale_fac[3] ;

  double *coo_som_fac_tmp = NULL;
  double  coo_point_dist_tmp[3] ;

  const double eps = 1e-15;


  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Calcul des coordonnées du barycentre B du polygone P */
  /*======================================================*/

  for (icoo = 0 ; icoo < 3 ; icoo++) {

    barycentre_fac[icoo] = 0. ;

    for (isom = 0 ; isom < nbr_som_fac ; isom++)
      barycentre_fac[icoo] += coo_som_fac[3*isom+icoo] ;

    barycentre_fac[icoo] /= nbr_som_fac ;

  }

  for (icoo = 0 ; icoo < 3 ; icoo++)
    normale_fac[icoo] = 0. ;

  /* Calcul de la normale */
  /*======================*/

  for (itri = 0 ; itri < nbr_som_fac ; itri++) {

    for (icoo = 0 ; icoo < 3 ; icoo++) {

      vect1[icoo] = coo_som_fac[3*itri+icoo] - barycentre_fac[icoo] ;

      if (itri < nbr_som_fac - 1)
        vect2[icoo] = coo_som_fac[3*(itri+1)+icoo] - barycentre_fac[icoo] ;
      else
        vect2[icoo] = coo_som_fac[icoo] - barycentre_fac[icoo] ;

    }

    normale_fac[0] += vect1[1] * vect2[2] - vect2[1] * vect1[2] ;
    normale_fac[1] += vect2[0] * vect1[2] - vect1[0] * vect2[2] ;
    normale_fac[2] += vect1[0] * vect2[1] - vect2[0] * vect1[1] ;

  }

  /* Projection dans un plan parallèle à la face */
  /*=============================================*/

  /* On ramène l'origine au centre de gravité de la fac */

  for (isom = 0 ; isom < nbr_som_fac ; isom++)
    for (icoo = 0 ; icoo < 3 ; icoo++)
      coo_som_fac[3*isom+icoo] -= barycentre_fac[icoo] ;

  for (icoo = 0 ; icoo < 3 ; icoo++)
    coo_point_dist[icoo] -= barycentre_fac[icoo] ;

  if (LocationToLocalMesh::abs(normale_fac[0]) > eps || LocationToLocalMesh::abs(normale_fac[1]) > eps) {

    /* Première rotation d'axe (Oz) et d'angle (Ox, proj normale sur Oxy) */

    coo_som_fac_tmp = new double [3 * nbr_som_fac];

    vect1[0] = 1. ;
    vect1[1] = 0. ;
    vect1[2] = 0. ;

    vect2[0] = normale_fac[0] ;
    vect2[1] = normale_fac[1] ;
    vect2[2] = 0. ;

    computeVectorProduct(prod_vect, vect1, vect2) ;

    cost = computeCrossProduct(vect1, vect2) / computeNorm(vect2) ;

    if (prod_vect[2] > 0.)
      sint =  computeNorm(prod_vect) / computeNorm(vect2) ;
    else
      sint = -computeNorm(prod_vect) / computeNorm(vect2) ;

    for (isom = 0 ; isom < nbr_som_fac ; isom++) {

      coo_som_fac_tmp[3*isom] =
         cost*coo_som_fac[3*isom] + sint*coo_som_fac[3*isom+1] ;
      coo_som_fac_tmp[3*isom+1] =
        -sint*coo_som_fac[3*isom] + cost*coo_som_fac[3*isom+1] ;
      coo_som_fac_tmp[3*isom+2] = coo_som_fac[3*isom+2] ;

    }

    coo_point_dist_tmp[0] = cost*coo_point_dist[0] + sint*coo_point_dist[1] ;
    coo_point_dist_tmp[1] = -sint*coo_point_dist[0] + cost*coo_point_dist[1] ;
    coo_point_dist_tmp[2] = coo_point_dist[2] ;

    /* Deuxième rotation d'axe (Oy) et d'angle (Oz', proj normale sur Ox'z) */

    vect1[0] =  0. ;
    vect1[1] =  0. ;
    vect1[2] =  1. ;

    vect2[0] =
      sqrt(normale_fac[0]*normale_fac[0] + normale_fac[1]*normale_fac[1]) ;
    vect2[1] = 0. ;
    vect2[2] = normale_fac[2] ;

    computeVectorProduct(prod_vect, vect1, vect2) ;

    cost = computeCrossProduct(vect1, vect2) / computeNorm(vect2) ;

    if (prod_vect[2] > 0.)
      sint =  computeNorm(prod_vect) / computeNorm(vect2) ;
    else
      sint = -computeNorm(prod_vect) / computeNorm(vect2) ;


    for (isom = 0 ; isom < nbr_som_fac ; isom++) {

      coo_som_fac[3*isom] =
         cost*coo_som_fac_tmp[3*isom] + sint*coo_som_fac_tmp[3*isom + 2] ;
      coo_som_fac[3*isom+1] = coo_som_fac_tmp[3*isom+1] ;
      coo_som_fac[3*isom+2] = 0. ;

    }

    coo_point_dist[0] =
      cost*coo_point_dist_tmp[0] + sint*coo_point_dist_tmp[2] ;
    coo_point_dist[1] = coo_point_dist_tmp[1] ;
    coo_point_dist[2] = 0. ;


    delete [] (coo_som_fac_tmp) ;

  }
  else {

    /* On écrase seulement la coordonnée z du sommet, en intervertissant
       éventuellement les coordonnées dans le plan de projection (Oxy).  */

    if (normale_fac[2] > 0.) {
      for (isom = 0 ; isom < nbr_som_fac ; isom++)
        coo_som_fac[3*isom+2] = 0. ;

      coo_point_dist[2] = 0. ;
    }
    else {
      for (isom = 0 ; isom < nbr_som_fac ; isom++) {
        coo_tmp = coo_som_fac[3*isom] ;
        coo_som_fac[3*isom] = coo_som_fac[3*isom+1] ;
        coo_som_fac[3*isom+1] = coo_tmp ;
        coo_som_fac[3*isom+2] = 0. ;
      }
      coo_tmp = coo_point_dist[0] ;
      coo_point_dist[0] = coo_point_dist[1] ;
      coo_point_dist[1] = coo_tmp ;
      coo_point_dist[2] = 0. ;
    }
  }
}

///
/// \brief Compute Mean Values
///
///

void LocationToLocalMesh::compute2DMeanValues()
{
  /* Boucle sur les points distants */

  const int n_dist_points = fvm::fvm_locator_get_n_dist_points(_fvmLocator);
  const fvm::fvm_lnum_t *dist_locations = fvm::fvm_locator_get_dist_locations(_fvmLocator);
  const fvm::fvm_coord_t *dist_coords = fvm::fvm_locator_get_dist_coords(_fvmLocator);
  fvm::fvm_coord_t coo_point_dist[3];

  /* Tableaux locaux */

  const double eps = 1e-15;
  std::vector <double> coo_som_fac;
  std::vector <double> s; 
  std::vector <double> dist;
  std::vector <double> aire;
  std::vector <double> proScal;

  int tailleDistBarCoords;

  _barycentricCoordinatesIndex = new std::vector <int> (n_dist_points + 1);

  tailleDistBarCoords = 4 * n_dist_points;
  _barycentricCoordinates = new std::vector <double> (tailleDistBarCoords);

  std::vector <int>& nDistBarCoords = *_barycentricCoordinatesIndex;
  std::vector <double>& distBarCoords = *_barycentricCoordinates;

  const int *meshConnectivityIndex = _supportMesh->getEltConnectivityIndex();
  const int *meshConnectivity = _supportMesh->getEltConnectivity();
  const double *meshVertexCoords = _supportMesh->getVertexCoords();

  nDistBarCoords[0] = 0;

  for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {
    //bft_printf("-- Etude du point : %d \n", ipoint);

    /* Initialisation - Copie locale */

    int isOnEdge = 0;
    int isVertex = 0;
    int ielt = dist_locations[ipoint] - 1;
    int nbr_som_fac =  meshConnectivityIndex[ielt+1] - 
                       meshConnectivityIndex[ielt];
    coo_point_dist[0] = dist_coords[3*ipoint];
    coo_point_dist[1] = dist_coords[3*ipoint + 1];
    coo_point_dist[2] = dist_coords[3*ipoint + 2];

    if (ipoint == 0) {
      coo_som_fac.resize(3 * nbr_som_fac);
      s.resize(3 * nbr_som_fac);
      dist.resize(nbr_som_fac);
      aire.resize(nbr_som_fac);
      proScal.resize(nbr_som_fac);
    }
    else
      if (proScal.size() < nbr_som_fac) {
        coo_som_fac.resize(3 * nbr_som_fac);
        s.resize(3 * nbr_som_fac);
        dist.resize(nbr_som_fac);
        aire.resize(nbr_som_fac);
        proScal.resize(nbr_som_fac);
      }

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      coo_som_fac[3*isom]   = meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)];
      coo_som_fac[3*isom+1] = meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+1];
      coo_som_fac[3*isom+2] = meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+2];
    }

    /* Projection sur un plan moyen */

    midplaneProjection(nbr_som_fac, &coo_som_fac[0], coo_point_dist);

    /* Calcul des coordonnnees barycentriques */

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      s[3*isom]   = coo_som_fac[3*isom]   - coo_point_dist[0];
      s[3*isom+1] = coo_som_fac[3*isom+1] - coo_point_dist[1];
      s[3*isom+2] = coo_som_fac[3*isom+2] - coo_point_dist[2];
      dist[isom] = sqrt(s[3*isom]*s[3*isom] +
                        s[3*isom+1]*s[3*isom+1] +
                        s[3*isom+2]*s[3*isom+2]);
    }

    int currentVertex;
    for (int isom = 0; isom < nbr_som_fac; isom++) {
      if (isom != (nbr_som_fac - 1)) {
        aire[isom] = s[3*isom]*s[3*(isom+1)+1] - s[3*(isom+1)]*s[3*isom+1];
        proScal[isom] = s[3*isom] * s[3*(isom+1)] +
                        s[3*isom+1] * s[3*(isom+1)+1] +
                        s[3*isom+2] * s[3*(isom+1)+2];
      }
      else {
        aire[isom] = s[3*isom]*s[1] - s[0]*s[3*isom+1];
        proScal[isom] = s[3*isom] * s[0] +
                        s[3*isom+1] * s[1] +
                        s[3*isom+2] * s[2];
      }
      if (dist[isom] <= eps) {
        isVertex = 1;
        currentVertex = isom;
        break;
      }
      /* faire un test avec eps pour proScal */
      else if (aire[isom] <= eps && proScal[isom] < 0.) {
        isOnEdge = 1;
        currentVertex = isom;
        break;
      }
    }

    /* Mise a jour de la taille du tableau de stockage des coordonnees barycentriques */

    nDistBarCoords[ipoint+1] = nDistBarCoords[ipoint] + nbr_som_fac;

    if (distBarCoords.size() <= nDistBarCoords[ipoint+1]) {
      distBarCoords.resize(2 * distBarCoords.size());
    }

    /* Le point distant est un sommet */

    if (isVertex) {

      for (int isom = 0; isom < nbr_som_fac; isom++)
        distBarCoords[nDistBarCoords[ipoint]+isom] = 0.;
      distBarCoords[nDistBarCoords[ipoint]+currentVertex] = 1.;
    }

    /* Le point distant est sur arete */

    else if (isOnEdge) {
      for (int isom = 0; isom < nbr_som_fac; isom++)
        distBarCoords[nDistBarCoords[ipoint]+isom] = 0.;

      int nextPoint;
      if (currentVertex == (nbr_som_fac - 1))
        nextPoint = 0;
      else
        nextPoint = currentVertex + 1;

      distBarCoords[nDistBarCoords[ipoint]+currentVertex] = dist[nextPoint]     / (dist[nextPoint]+dist[currentVertex]);
      distBarCoords[nDistBarCoords[ipoint]+nextPoint]     = dist[currentVertex] / (dist[nextPoint]+dist[currentVertex]);

    }

    /* Cas general */

    else {
      double sigma = 0;
      for (int isom = 0; isom < nbr_som_fac; isom++) {
        double coef = 0.;
        int previousVertex;
        int nextVertex;
        if (isom != 0)
          previousVertex = isom - 1;
        else
          previousVertex = nbr_som_fac - 1;
        if (isom < nbr_som_fac - 1)
          nextVertex = isom + 1;
        else
          nextVertex = 0;
        if (LocationToLocalMesh::abs(aire[previousVertex]) > eps)
          coef += (dist[previousVertex] - proScal[previousVertex]/dist[isom]) / aire[previousVertex];
        // BUG: verifier le test de calcul du coeff
        if (LocationToLocalMesh::abs(aire[isom]) > eps)
          coef += (dist[nextVertex] - proScal[isom]/dist[isom]) / aire[isom];
        sigma += coef;
        distBarCoords[nDistBarCoords[ipoint]+isom] = coef;
      }
      if (LocationToLocalMesh::abs(sigma) >= eps ) {
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          distBarCoords[nDistBarCoords[ipoint]+isom] /= sigma;
        }
      }
      else {
        double abs_sigma = LocationToLocalMesh::abs(sigma);
        printf("Warning : mise à NAN %f %f", abs_sigma,  eps);
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          distBarCoords[nDistBarCoords[ipoint]+isom] = NAN;
        }
      }
    }

    if (1 == 1) {
      bft::bft_printf("coord %i :", ipoint);
      bft::bft_printf(" %12.5e %12.5e %12.5e", dist_coords[3*ipoint], 
                      dist_coords[3*ipoint+1], 
                      dist_coords[3*ipoint+2] );
      bft::bft_printf("\n");

      bft::bft_printf("coo b %i :", ipoint);
      for (int isom = 0; isom < nbr_som_fac; isom++) {
        bft::bft_printf(" %f", distBarCoords[nDistBarCoords[ipoint]+isom]);
      }
      bft::bft_printf("\n");
    }
  }

  coo_som_fac.clear();
  s.clear();
  aire.clear();
  dist.clear();
  proScal.clear();

  distBarCoords.resize(nDistBarCoords[n_dist_points]);

}



} // Namespace CWIPI

