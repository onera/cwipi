/*
 * DistantLocation.cxx
 *
 *  Created on: Oct 16, 2009
 *      Author: equemera
 */

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
/// \brief Clear location
///

void LocationToDistantMesh::clear()
{
  if (!_isCoupledRank && _couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {
    if (_locatedPoint != NULL)
      delete []  _locatedPoint;

    if (_unlocatedPoint != NULL)
      delete []  _unlocatedPoint;
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

}

}
