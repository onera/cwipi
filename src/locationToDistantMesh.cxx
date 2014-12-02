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
/*
 * DistantLocation.cxx
 *
 *  Created on: Oct 16, 2009
 *      Author: equemera
 */


#include <mpi.h>
#include <fvmc_locator.h>

#include "locationToDistantMesh.hxx"

namespace cwipi
{

///
/// \brief Default constructor
///

LocationToDistantMesh::LocationToDistantMesh(const bool isCoupledRank,
                                             const cwipi_coupling_type_t couplingType,
                                             const ApplicationProperties& localApplicationProperties)
: _isCoupledRank(isCoupledRank), _couplingType(couplingType),
  _localApplicationProperties(localApplicationProperties)
{
  _coordsPointsToLocate = NULL;
  _nPointsToLocate = 0;
  _nUnlocatedPoint = 0 ;
  _nLocatedPoint = 0;
  _unlocatedPoint = NULL;
  _locatedPoint = NULL;
  _elementContaining = NULL;
  _elementContainingBarycentricCoordinates = NULL;
  _elementContainingMPIrankContaining = NULL;
  _elementContainingNVertex = NULL;
  _elementContainingVertex = NULL;
  _elementContainingVertexCoords = NULL;
  _locationInfo = CWIPI_BASIC_INFO;
  _toLocate = true;

}

///
/// \brief Default destructor
///

LocationToDistantMesh::~LocationToDistantMesh()
{
  clear();
}

///
/// \brief Set points to locate
///
///   @param nPointsToLocate      Number of points to locate
///   @param coordsPointsToLocate Coordinates of points to locate
///

void LocationToDistantMesh::setpointsToLocate(int nPointsToLocate, double *coordsPointsToLocate)
{
  _toLocate = true;
  if (_isCoupledRank) {
    _nPointsToLocate = nPointsToLocate;
    _coordsPointsToLocate = coordsPointsToLocate;
  }
}

///
/// \brief Set points to locate
///
///   @param nPointsToLocate      Number of points to locate
///   @param coordsPointsToLocate Coordinates of points to locate
///

void LocationToDistantMesh::synchronize()
{
  const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

  int rootRank;

  if (_isCoupledRank) {
    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      MPI_Comm_rank(localComm, &rootRank);

      assert(rootRank == 0);

      MPI_Bcast(&_nPointsToLocate, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast(&_nUnlocatedPoint, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast(&_nLocatedPoint, 1, MPI_INT, rootRank, localComm );

      MPI_Bcast(_locatedPoint,
                _nLocatedPoint,
                MPI_INT,
                rootRank,
                localComm );

      MPI_Bcast(_unlocatedPoint, _nUnlocatedPoint, MPI_INT, rootRank, localComm );

      _toLocate = false;

      if (_locationInfo == CWIPI_DISTANT_MESH_INFO) {


      }
    }
  }

  else {
    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      rootRank = 0;

      _toLocate = false;

      MPI_Bcast(&_nPointsToLocate, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast(&_nUnlocatedPoint, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast(&_nLocatedPoint, 1, MPI_INT, rootRank, localComm );

      if (_locatedPoint != NULL)
        delete [] _locatedPoint;

      if (_unlocatedPoint != NULL)
        delete [] _unlocatedPoint ;

      _locatedPoint = new int[_nLocatedPoint];
      _unlocatedPoint = new int[_nUnlocatedPoint];

      MPI_Bcast(_locatedPoint,
                _nLocatedPoint,
                MPI_INT,
                rootRank,
                localComm);
      MPI_Bcast( _unlocatedPoint, _nUnlocatedPoint, MPI_INT, rootRank, localComm );

      if (_locationInfo == CWIPI_DISTANT_MESH_INFO) {

        if (_elementContaining != NULL)
          delete [] _elementContaining;

        if (_elementContainingNVertex != NULL)
          delete [] _elementContainingNVertex ;

        if (_elementContainingBarycentricCoordinates != NULL)
          delete [] _elementContainingBarycentricCoordinates ;

        if (_elementContainingVertex != NULL)
          delete [] _elementContainingVertex ;

        if (_elementContainingVertexCoords != NULL)
          delete [] _elementContainingVertexCoords ;

        _elementContaining = new int [_nLocatedPoint];
        _elementContainingNVertex = new int [_nLocatedPoint+1];

        MPI_Bcast(_elementContaining,
                  _nLocatedPoint,
                  MPI_INT,
                  rootRank,
                  localComm);

        MPI_Bcast(_elementContainingNVertex,
                  _nLocatedPoint + 1,
                  MPI_INT,
                  rootRank,
                  localComm);

        _elementContainingBarycentricCoordinates = new double [_elementContainingNVertex[_nLocatedPoint]];
        _elementContainingVertex = new int [_elementContainingNVertex[_nLocatedPoint]];
        _elementContainingVertexCoords = new double [3 * _elementContainingNVertex[_nLocatedPoint]];

        MPI_Bcast(_elementContainingBarycentricCoordinates,
                  _elementContainingNVertex[_nLocatedPoint],
                  MPI_DOUBLE,
                  rootRank,
                  localComm);

        MPI_Bcast(_elementContainingVertex,
                  _elementContainingNVertex[_nLocatedPoint],
                  MPI_INT,
                  rootRank,
                  localComm);

        MPI_Bcast(_elementContainingVertexCoords,
                  3 *_elementContainingNVertex[_nLocatedPoint],
                  MPI_DOUBLE,
                  rootRank,
                  localComm);

      }
    }
  }
}

///
/// 
///
size_t LocationToDistantMesh::locationSize()
{
  size_t il_size = 0;
  il_size += 3 * sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingNVertex != NULL) 
    il_size += (_nLocatedPoint+1) * sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingVertex != NULL) 
    il_size += _elementContainingNVertex[_nLocatedPoint] * sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingBarycentricCoordinates != NULL) 
    il_size += _elementContainingNVertex[_nLocatedPoint] * sizeof(double);
  
  il_size += sizeof(int);
  if(_elementContainingMPIrankContaining != NULL) 
    il_size += _nLocatedPoint *sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingVertexCoords != NULL) 
    il_size += (3 * _elementContainingNVertex[_nLocatedPoint]) * sizeof(double);

  il_size += sizeof(int);
  if (_elementContaining != NULL) 
    il_size += _nLocatedPoint * sizeof(int);
  
  return il_size;
  
}

void LocationToDistantMesh::packLocation(unsigned char *buff)
{
  int s;
  size_t cur_pos;
  void *p;
  p = (void *)buff;

  p = mempcpy(p,(void *)&_nPointsToLocate,sizeof(int));
  p = mempcpy(p,(void *)&_nLocatedPoint,sizeof(int));
  p = mempcpy(p,(void *)&_nUnlocatedPoint,sizeof(int));
  
  // pour chaque tableau, on commence par stoker sa taille 
  // pour pouvoir l'allouer si nécessaire à la lecture
  
  if (_elementContainingNVertex != NULL) {
    s = _nLocatedPoint+1;
    p = mempcpy(p,(void *)&s, sizeof(int));
    p = mempcpy(p,(void *) _elementContainingNVertex,s*sizeof(int));
  } else {
    s = 0;
    p = mempcpy(p,(void *)&s, sizeof(int));
  }

  if (_elementContainingVertex != NULL) {
    s = _elementContainingNVertex[_nLocatedPoint];
    p = mempcpy(p,(void *)&s, sizeof(int));
    p = mempcpy(p,(void *) _elementContainingVertex,s*sizeof(int));
  } else {
    s = 0;
    p = mempcpy(p,(void *)&s, sizeof(int));
  }

  if (_elementContainingBarycentricCoordinates != NULL) {
    s = _elementContainingNVertex[_nLocatedPoint];
    p = mempcpy(p,(void *)&s, sizeof(int));
    p = mempcpy(p,(void *) _elementContainingBarycentricCoordinates,s*sizeof(double));
  } else {
    s = 0;
    p = mempcpy(p,(void *)&s, sizeof(int));
  }

  // attention : alloué dans locationToLocalMesh
  if(_elementContainingMPIrankContaining != NULL) {
    s = _nLocatedPoint;
    p = mempcpy(p,(void *)&s, sizeof(int));
    p = mempcpy(p,(void *) _elementContainingMPIrankContaining,s*sizeof(int));
  } else {
    s = 0;
    p = mempcpy(p,(void *)&s, sizeof(int));
  }
  
  if (_elementContainingVertexCoords != NULL) {
    s = (3 * _elementContainingNVertex[_nLocatedPoint]);
    p = mempcpy(p,(void *)&s, sizeof(int));
    p = mempcpy(p,(void *)_elementContainingVertexCoords,s*sizeof(double));
  } else {
    s = 0;
    p = mempcpy(p,(void *)&s, sizeof(int));
  }

  if (_elementContaining != NULL) {
    s = _nLocatedPoint;
    p = mempcpy(p,(void *)&s, sizeof(int));
    p = mempcpy(p,(void *)_elementContaining,_nLocatedPoint*sizeof(int));
  } else {
    s = 0;
    p = mempcpy(p,(void *)&s, sizeof(int));
  }
}

void LocationToDistantMesh::unpackLocation(unsigned char *buff)
{
  int s;
  size_t cur_pos;
  
  cur_pos = 0;

  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_nPointsToLocate,sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_nLocatedPoint,sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_nUnlocatedPoint,sizeof(int));
 
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingNVertex != NULL) delete [] _elementContainingNVertex;
    _elementContainingNVertex = new int[s];
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)_elementContainingNVertex,s*sizeof(int));
  }
    
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingVertex != NULL) delete [] _elementContainingVertex;
    _elementContainingVertex = new int[s];
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *) _elementContainingVertex,s*sizeof(int));
  }

  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingBarycentricCoordinates != NULL) delete [] _elementContainingBarycentricCoordinates;
    _elementContainingBarycentricCoordinates = new double[s];
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *) _elementContainingBarycentricCoordinates,s*sizeof(double));
  } 

  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if(_elementContainingMPIrankContaining != NULL) delete [] _elementContainingMPIrankContaining;
    _elementContainingMPIrankContaining = new int[s];
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *) _elementContainingMPIrankContaining,s*sizeof(int));
  } 
 
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingVertexCoords != NULL) delete [] _elementContainingVertexCoords;
    _elementContainingVertexCoords = new double[s] ;
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)_elementContainingVertexCoords,s*sizeof(double));
  } 
  
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContaining != NULL) delete [] _elementContaining;
    _elementContaining = new int[s];
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)_elementContaining,_nLocatedPoint*sizeof(int));
  } 
  _toLocate = false;
}




///
/// \brief Clear location
///

void LocationToDistantMesh::clear()
{
  if (!_isCoupledRank && _couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {
    if (_locatedPoint != NULL) {
      delete []  _locatedPoint;
      _locatedPoint = NULL;
    }
    

    if (_unlocatedPoint != NULL) {
      delete []  _unlocatedPoint;
      _unlocatedPoint = NULL;
    }

  }

  if (_elementContainingBarycentricCoordinates != NULL)
    delete [] _elementContainingBarycentricCoordinates;
  if (_elementContainingMPIrankContaining != NULL)
    delete [] _elementContainingMPIrankContaining;
  if (_elementContainingNVertex != NULL)
    delete [] _elementContainingNVertex;
  if (_elementContainingVertex != NULL)
    delete [] _elementContainingVertex;
  if (_elementContainingVertexCoords != NULL)
    delete [] _elementContainingVertexCoords;
  if (_elementContaining != NULL)
    delete [] _elementContaining;

  _elementContainingBarycentricCoordinates = NULL;
  _elementContainingMPIrankContaining = NULL;
  _elementContainingNVertex = NULL;
  _elementContainingVertex = NULL;
  _elementContainingVertexCoords = NULL;
  _elementContaining = NULL;

}

}
