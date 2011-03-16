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
#include "coo_baryc.h"


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
         compute3DMeanValues();
        //TODO: calcul des coord barycentriques 3D
        //int nPoints;        
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

void LocationToLocalMesh::compute3DMeanValues()
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
  std::vector <double> proScal;
  std::vector <double> angle;
  std::vector <double> normale;
  std::vector <double> distBarCoordsTmp;
  std::vector <int> localIndVertex(_supportMesh->getNVertex());
    
  int tailleDistBarCoords;
  int isOnFace;
  int ind_fac = 0;
  int ind_fac_som = 0;

  const int eltStd = _supportMesh->getNElts() - _supportMesh->getNPolyhedra();
  printf("Nelt  eltSrd   Npolyhedra   %d   %d   %d  \n",_supportMesh->getNElts(),eltStd,_supportMesh->getNPolyhedra());
  printf("nbr vertex   %d \n",_supportMesh->getNVertex());
  _barycentricCoordinatesIndex = new std::vector <int> (n_dist_points + 1);

  tailleDistBarCoords = 4 * n_dist_points;
  _barycentricCoordinates = new std::vector <double> (tailleDistBarCoords);

  std::vector <int>& nDistBarCoords = *_barycentricCoordinatesIndex;
  std::vector <double>& distBarCoords = *_barycentricCoordinates;

  const int *meshConnectivityIndex;
  const int *meshConnectivity;
  const double *meshVertexCoords = _supportMesh->getVertexCoords();
  const int *face_index;
  const int *cell_to_face_connectivity;
  const int *face_connectivity;
  const int *face_connectivity_index;

  int *face_index_Tmp;
  int *cell_to_face_connectivity_Tmp;
  int *face_connectivity_Tmp;
  int *face_connectivity_index_Tmp;

  nDistBarCoords[0] = 0;

  printf("npoint  %d\n",n_dist_points);
  for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {
    
    int ielt = dist_locations[ipoint] - 1;
    isOnFace = 0;

    coo_point_dist[0] = dist_coords[3 * ipoint];
    coo_point_dist[1] = dist_coords[3 * ipoint + 1];
    coo_point_dist[2] = dist_coords[3 * ipoint + 2];

    int nbr_face;
    int nbr_som_fac;      
    int nbr_som;

    if(ielt < eltStd){ 
      
      printf("**********dodecaedre********** \n");
      meshConnectivityIndex = _supportMesh->getEltConnectivityIndex();
      meshConnectivity = _supportMesh->getEltConnectivity();

      nbr_som = meshConnectivityIndex[ielt+1] - meshConnectivityIndex[ielt];
      nbr_face = 4;
      nbr_som_fac = 3;      
            
      face_connectivity_index_Tmp = (int* ) malloc((nbr_face + 1) * sizeof(int) );
      face_connectivity_Tmp = (int* ) malloc(nbr_som_fac * nbr_face * sizeof(int));
      cell_to_face_connectivity_Tmp = (int* ) malloc(nbr_face * sizeof(int));

      face_connectivity_Tmp[0]  = 1;
      face_connectivity_Tmp[1]  = 3;
      face_connectivity_Tmp[2]  = 2;
      
      face_connectivity_Tmp[3]  = 1;
      face_connectivity_Tmp[4]  = 2;
      face_connectivity_Tmp[5]  = 4;
      
      face_connectivity_Tmp[6]  = 1;
      face_connectivity_Tmp[7]  = 4;
      face_connectivity_Tmp[8]  = 3;
      
      face_connectivity_Tmp[9]  = 2;
      face_connectivity_Tmp[10] = 3;
      face_connectivity_Tmp[11] = 4;
      
      for(int i = 0; i < nbr_face ; i++)
        cell_to_face_connectivity_Tmp[i] = i+1;
      
      for (int i = 0; i < nbr_face ; i++)
        face_connectivity_index_Tmp[i] = 3*i;
      
      face_connectivity = face_connectivity_Tmp;
      cell_to_face_connectivity = cell_to_face_connectivity_Tmp;
      face_connectivity_index = face_connectivity_index_Tmp;
 
       }
    else {
      //printf("ielt   %d \n",ielt);
      printf("**********octaedre********** \n");
      ielt -= eltStd;

      //printf("ielt - eltStd   %d \n",ielt);

      meshConnectivityIndex = &(_supportMesh->getPolyhedraCellToVertexConnectivityIndex()[0]);
      meshConnectivity = &(_supportMesh->getPolyhedraCellToVertexConnectivity()[0]);
      face_index = _supportMesh->getPolyhedraFaceIndex();     
      cell_to_face_connectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();     
      face_connectivity = _supportMesh->getPolyhedraFaceConnectivity();     
      face_connectivity_index = _supportMesh->getPolyhedraFaceConnectivityIndex();
     
      nbr_som = meshConnectivityIndex[ielt + 1] - meshConnectivityIndex[ielt];
      nbr_face = face_index[ielt + 1] - face_index[ielt];
      nbr_som_fac = 3;     
      ind_fac = face_index[ielt];
      //printf("meshConnectivityIndex    %d   %d  \n",meshConnectivityIndex[0],meshConnectivityIndex[1]);
      }
  
    for (int isom = 0 ; isom < nbr_som ; isom++){
      localIndVertex[meshConnectivity[meshConnectivityIndex[ielt] + isom]-1] = isom ;     
    }
        
    nDistBarCoords[ipoint + 1] = nDistBarCoords[ipoint] + nbr_som;

    if (distBarCoords.size() <= nDistBarCoords[ipoint + 1])
      distBarCoords.resize(2 * distBarCoords.size());
 
    if(ipoint == 0) {
      angle.resize(nbr_som_fac); // danger lorsque nbr_som_fac != 3
      normale.resize(3 * nbr_som_fac); // danger lorsque nbr_som_fac != 3
      distBarCoordsTmp.resize(nbr_som);
      coo_som_fac.resize(3 * nbr_som);        
      dist.resize(nbr_som);        
      s.resize(3 * nbr_som);

      }
    else
      if (dist.size() <  nbr_som_fac){
        coo_som_fac.resize(3 * nbr_som_fac);        
        dist.resize(nbr_som_fac);
        s.resize(3 * nbr_som);
      } 


    double sigma = 0;


    for (int isom = 0; isom < nbr_som; isom++)
      distBarCoordsTmp[isom] = 0; 
    
    for(int i = 0 ;i < _supportMesh->getNVertex() ;i++){
      printf("meshVertexCoords   %f   %f   %f  \n",meshVertexCoords[3 *i],meshVertexCoords[3 *i+1],meshVertexCoords[3 *i+2]);
    } 

    //Verifier que le sommet appartient ni a une face, arrete ou egal a un sommet
     for (int isom = 0; isom < nbr_som ; isom++){
       // printf("(meshConnectivityIndex[ielt] + isom)   %d \n",(meshConnectivityIndex[ielt] + isom));
       
       coo_som_fac[3 * isom]     = meshVertexCoords[3 * (meshConnectivity[meshConnectivityIndex[ielt] + isom]-1)];
       coo_som_fac[3 * isom + 1] = meshVertexCoords[3 * (meshConnectivity[meshConnectivityIndex[ielt] + isom]-1) + 1];
       coo_som_fac[3 * isom + 2] = meshVertexCoords[3 * (meshConnectivity[meshConnectivityIndex[ielt] + isom]-1) + 2];

       // printf("vertex coord  %d   %f   %f   %f \n", isom, coo_som_fac[3 * isom], coo_som_fac[3 * isom+1], coo_som_fac[3 * isom+2]);
       
       s[3 * isom]     = coo_som_fac[3 * isom]     - coo_point_dist[0];
       s[3 * isom + 1] = coo_som_fac[3 * isom + 1] - coo_point_dist[1];
       s[3 * isom + 2] = coo_som_fac[3 * isom + 2] - coo_point_dist[2];

       
       //printf("s %d   %f   %f   %f \n", isom, s[3 * isom], s[3 * isom+1], s[3 * isom+2]);

       dist[isom] = sqrt(s[3 * isom]*s[3*isom] +
                         s[3 * isom + 1]*s[3*isom + 1] +
                         s[3 * isom + 2]*s[3*isom + 2]);           
     }
    
    for(int iface = 0; iface < nbr_face ; iface++){
      double det;

      ind_fac_som = face_connectivity_index[cell_to_face_connectivity[ind_fac + iface] - 1];

      const int i = localIndVertex[face_connectivity[ind_fac_som] - 1];
      const int j = localIndVertex[face_connectivity[ind_fac_som + 1] - 1];
      const int k = localIndVertex[face_connectivity[ind_fac_som + 2] - 1];

      // printf("--- sommet   %d  %d  %d \n",i,j,k);
      /*   printf("s  i  j   k  %f  %f  %f \n %f  %f  %f \n %f  %f  %f \n",
             s[3*i],s[3*i+1],s[3*i+2],
             s[3*j],s[3*j+1],s[3*j+2],
             s[3*k],s[3*k+1],s[3*k+2]);*/


      det =  (1./6) * (  ( (s[3*i+1] * s[3*j+2] - s[3*j+1] * s[3*i+2]) * s[3*k])
                       + ( (s[3*j]  * s[3*i+2] - s[3*i] * s[3*j+2]) * s[3*k+1])
                       + ( (s[3*i]  * s[3*j+2] - s[3*j] * s[3*i+1]) * s[3*k+2]) );

      // printf("determinant   %f \n",det);
      if(abs(det) < eps){

        double aireTri_ijk;
        double aireTri_ijv;
        double aireTri_ikv;
        double aireTri_jkv;
        double coo_ijx;
        double coo_ijy;
        double coo_ijz;
        double coo_ikx;
        double coo_iky;
        double coo_ikz;
       
        coo_ijx = coo_som_fac[3*j]   - coo_som_fac[3*i];
        coo_ijy = coo_som_fac[3*j+1] - coo_som_fac[3*i+1];
        coo_ijz = coo_som_fac[3*j+2] - coo_som_fac[3*i+2];
        coo_ikx = coo_som_fac[3*k]   - coo_som_fac[3*i];
        coo_iky = coo_som_fac[3*k+1] - coo_som_fac[3*i+1];
        coo_ikz = coo_som_fac[3*k+2] - coo_som_fac[3*i+2];
        
        aireTri_ijk = sqrt(  ( coo_ijy * coo_ikz - coo_ijz * coo_iky ) * ( coo_ijy * coo_ikz - coo_ijz * coo_iky )
                           + ( coo_ijz * coo_ikx - coo_ijx * coo_ikz ) * ( coo_ijz * coo_ikx - coo_ijx * coo_ikz )
                           + ( coo_ijx * coo_iky - coo_ijy * coo_ikx ) * ( coo_ijx * coo_iky - coo_ijy * coo_ikx ));

        aireTri_ijv = sqrt(  ( s[3*i+1] * s[3*j+2] - s[3*i+2] * s[3*j+1] ) * ( s[3*i+1] * s[3*j+2] - s[3*i+2] * s[3*j+1] )
                           + ( s[3*i+2] * s[3*j]   - s[3*i]   * s[3*j+2] ) * ( s[3*i+2] * s[3*j]   - s[3*i]   * s[3*j+2] )
                           + ( s[3*i]   * s[3*j+1] - s[3*i+1] * s[3*j] )   * ( s[3*i]   * s[3*j+1] - s[3*i+1] * s[3*j] ));

        aireTri_ikv = sqrt(  ( s[3*i+1] * s[3*k+2] - s[3*i+2] * s[3*k+1] ) * ( s[3*i+1] * s[3*k+2] - s[3*i+2] * s[3*k+1] )
                           + ( s[3*i+2] * s[3*k]   - s[3*i]   * s[3*k+2] ) * ( s[3*i+2] * s[3*k]   - s[3*i]   * s[3*k+2] )
                           + ( s[3*i]   * s[3*k+1] - s[3*i+1] * s[3*k] )   * ( s[3*i]   * s[3*k+1] - s[3*i+1] * s[3*k] ));

        aireTri_jkv = sqrt(  ( s[3*j+1] * s[3*k+2] - s[3*j+2] * s[3*k+1] ) * ( s[3*j+1] * s[3*k+2] - s[3*j+2] * s[3*k+1] )
                           + ( s[3*j+2] * s[3*k]   - s[3*j]   * s[3*k+2] ) * ( s[3*j+2] * s[3*k]   - s[3*j]   * s[3*k+2] )
                           + ( s[3*j]   * s[3*k+1] - s[3*j+1] * s[3*k] )   * ( s[3*j]   * s[3*k+1] - s[3*j+1] * s[3*k] ));
        
        
        /* printf("aire ijk  %f  \n",aireTri_ijk);
        printf("aire ijv  %f  \n",aireTri_ijv);
        printf("aire ikv  %f  \n",aireTri_ikv);
        printf("aire jkv  %f  \n",aireTri_jkv);*/

        for(int isom = 0 ; isom < nbr_som ; isom++)
          distBarCoords[nDistBarCoords[ipoint]+isom] = 0.;
          
        distBarCoords[ nDistBarCoords[ipoint] + i ] = aireTri_jkv / aireTri_ijk ;
        distBarCoords[ nDistBarCoords[ipoint] + j ] = aireTri_ikv / aireTri_ijk ;
        distBarCoords[ nDistBarCoords[ipoint] + k ] = aireTri_ijv / aireTri_ijk ;

        isOnFace = 1;

        break;

      }

    }

    if(!isOnFace){
      int ind_som;
        // normalisation de s (on sait que s est non nul) 

        for(int isom = 0 ; isom < nbr_som ; isom++){
          // printf("s ind %d  %f  %f  %f\n",isom,s[3*isom],s[3*isom+1],s[3*isom+2]);
          s[3 * isom]     /= dist[isom];
          s[3 * isom + 1] /= dist[isom];
          s[3 * isom + 2] /= dist[isom];       
        }              

      for(int iface = 0; iface < nbr_face ; iface++){        
        
        ind_fac_som = face_connectivity_index[cell_to_face_connectivity[ind_fac + iface] - 1] ;
      
        /*for(int isom = 0 ; isom < nbr_som_fac ; isom++){
          printf("nface %d  %d \n",nbr_som_fac*iface+isom,face_connectivity[nbr_som_fac*iface+isom]);
          printf("mesh coord   %f   %f   %f \n",
          meshVertexCoords[3*(meshConnectivityIndex[ielt] + face_connectivity[nbr_som_fac*iface + isom] - 1)],
          meshVertexCoords[3*(meshConnectivityIndex[ielt] + face_connectivity[nbr_som_fac*iface + isom] - 1) + 1],
          meshVertexCoords[3*(meshConnectivityIndex[ielt] + face_connectivity[nbr_som_fac*iface + isom] - 1) + 2]);                           
          }*/

        for(int isom = 0; isom < nbr_som_fac; isom++){
                   
          int isuiv;
          int iprec;
          double prod_scal;
          double mod;                    

          if(isom ==0)
            iprec = localIndVertex[face_connectivity[ ind_fac_som +  nbr_som_fac - 1 ] - 1]; // se place a la ieme face sur le dernier element
          else 
            iprec = localIndVertex[face_connectivity[ind_fac_som + isom - 1 ] - 1]; //prend l'indice du sommet juste precedent
          
          if (isom == nbr_som_fac - 1)
            isuiv = localIndVertex[face_connectivity[ind_fac_som] - 1];
          else 
            isuiv = localIndVertex[face_connectivity[ind_fac_som + isom + 1] - 1];       
                                   
          prod_scal = s[3 * iprec]     * s[3 * isuiv]
                    + s[3 * iprec + 1] * s[3 * isuiv + 1]
                    + s[3 * iprec + 2] * s[3 * isuiv + 2];

          angle[isom] = acos(prod_scal); //s est de norme 1         

          normale[3 * isom]     = s[3 * iprec + 1] * s[3 * isuiv + 2] - s[3 * iprec + 2] * s[3*isuiv + 1];        
          normale[3 * isom + 1] = s[3 * iprec + 2] * s[3 * isuiv]     - s[3 * iprec]     * s[3 * isuiv + 2];
          normale[3 * isom + 2] = s[3 * iprec]     * s[3 * isuiv + 1] - s[3 * iprec + 1] * s[3 * isuiv];        
                             
          mod = sqrt(normale[3 * isom]*normale[3 * isom]
                     + normale[3 * isom + 1]*normale[3 * isom + 1]
                     + normale[3 * isom + 2]*normale[3 * isom + 2]);
          
          normale[3 * isom]     /= mod ;
          normale[3 * isom + 1] /= mod ;
          normale[3 * isom + 2] /= mod ;
            
        }           
        
        for(int isom = 0; isom < nbr_som_fac; isom++){
           
          int isuiv;
          int iprec;
          double ps_nij_njk; //a ameliorer
          double ps_nki_njk; //a ameliorer
          double ps_ei_njk;  //a ameliorer          

          if(isom ==0)
            iprec =  nbr_som_fac - 1; // se place a la ieme face sur le dernier element
          else 
            iprec = isom - 1; //prend l'indice du sommet juste precedent
          
          if (isom == nbr_som_fac - 1)
            isuiv = 0;
          else 
            isuiv = isom + 1;            
                  
          ps_nij_njk = normale[3 * isom] * normale[3 * isuiv]
            + normale[3 * isom + 1] * normale[3 * isuiv + 1]
            + normale[3 * isom + 2] * normale[3 * isuiv + 2];
          
          ps_nki_njk = normale[3 * isom] * normale[3 * iprec]
            + normale[3 * isom + 1] * normale[3 * iprec + 1]
            + normale[3 * isom + 2] * normale[3 * iprec + 2];
          
          ps_ei_njk = s[3 * (localIndVertex[face_connectivity[ind_fac_som + isom] - 1])] * normale[3 * isom] 
            + s[3 * (localIndVertex[face_connectivity[ind_fac_som + isom] - 1]) + 1] * normale[3 * isom + 1]
            + s[3 * (localIndVertex[face_connectivity[ind_fac_som + isom] - 1]) + 2] * normale[3 * isom + 2];       
          
          distBarCoordsTmp[localIndVertex[(face_connectivity[ind_fac_som + isom] - 1)]] += (angle[isom] + angle[isuiv] * ps_nij_njk + angle[iprec] * ps_nki_njk) / (2 * ps_ei_njk);           
          //printf("(face_connectivity[ind_fac_som + isom] - 1)   %d \n",(face_connectivity[ind_fac_som + isom] - 1));
        }
          
      }
      
      for(int isom = 0; isom<nbr_som;isom++){
        // printf("isom distbarcoordtm   %d   %f   \n",isom,distBarCoordsTmp[isom]);
        sigma += distBarCoordsTmp[isom];      
      }      
      //printf("sigma   %f   \n",sigma);
      
      for(int isom = 0; isom<nbr_som;isom++)
        distBarCoords[nDistBarCoords[ipoint] + isom] = distBarCoordsTmp[isom]/sigma;
      
      
    }
    
    for(int isom = 0; isom<nbr_som;isom++)
      printf("sommet   %f %f %f  coord bar final   %f \n",coo_som_fac[3*isom],coo_som_fac[3*isom+1],coo_som_fac[3*isom+2],distBarCoords[nDistBarCoords[ipoint] + isom]);
    
    
    if( 1== 1){
      std::vector <double> test(3);
      
      for(int i =0;i<3;i++)
        test[i] = 0;
      
      for(int isom =0;isom<nbr_som;isom++){
        
        test[0] += distBarCoords[nDistBarCoords[ipoint] + isom] * coo_point_dist[0];
        test[1] += distBarCoords[nDistBarCoords[ipoint] + isom] * coo_point_dist[1];
        test[2] += distBarCoords[nDistBarCoords[ipoint] + isom] * coo_point_dist[2];
        
      }
      for(int i =0;i<3;i++)
        printf("test   %f  coord  %f \n",test[i],coo_point_dist[i]);

    }
   
    if (1 == 1) {
      bft::bft_printf("coord %i :", ipoint);
      bft::bft_printf(" %12.5e %12.5e %12.5e", dist_coords[3 * ipoint], 
                      dist_coords[3 * ipoint + 1], 
                      dist_coords[3 * ipoint + 2] );
      bft::bft_printf("\n");
      
      bft::bft_printf("coo b %i :", ipoint);
      for (int isom = 0; isom < nbr_som; isom++) 
        bft::bft_printf(" %f", distBarCoords[nDistBarCoords[ipoint] + isom]);
      
      bft::bft_printf("\n");
      
      }
    
    if(nbr_som == 4){    
      free(cell_to_face_connectivity_Tmp);
      free(face_connectivity_index_Tmp);
      free(face_connectivity_Tmp);
    }
    
    
  }
  coo_som_fac.clear();
  s.clear();
  dist.clear();
  proScal.clear();
  angle.clear();
  distBarCoordsTmp.clear();
  localIndVertex.clear();
}

} // Namespace CWIPI

