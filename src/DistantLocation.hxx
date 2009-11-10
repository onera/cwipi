/*
 * LocationToDistantMesh.hxx
 *
 *  Created on: Oct 16, 2009
 *      Author: equemera
 */

#ifndef LOCATIONTODISTANTMESH_HXX_
#define LOCATIONTODISTANTMESH_HXX_

#include <string>
#include <vector>

#include "cwipi.h"

namespace cwipi
{

///
/// \brief This Class realize the location of point in a distant mesh
///
/// This class must be used with symetric class LocationToLocalMesh

class LocationToDistantMesh
{
public:

  ///
  /// \brief Default constructor
  ///

  LocationToDistantMesh();

  ///
  /// \brief Default destructor
  ///

  virtual ~LocationToDistantMesh();

  ///
  /// \brief Set information type to receive with location result
  ///
  ///   @param information type to receive with location result
  ///

  void setLocatedPointsInfo(cwipi_located_point_info_t info);

  ///
  /// \brief Set points to locate
  ///
  ///   @param nPointsToLocate      Number of points to locate
  ///   @param coordsPointsToLocate Coordinates of points to locate
  ///

  void setpointsToLocate(int nPointsToLocate, double *coordsPointsToLocate);

  ///
  /// \brief Points location
  ///

  void locate();

  ///
  /// \brief Set information to exchange
  ///



private :
  int     _nPointsToLocate;                       ///< Number of points to locate
  double *_coordsPointsToLocate;                  ///< Coordinates of points to locate
  int     _nUnlocatedPoint;                       ///< Number of unlocated points
  int     _nLocatedPoint;                         ///< Number of located points
  int    *_unlocatedPoint;                        ///< Unlocated points
  int    *_locatedPoint;                          ///< Located points

  int    *_nElementContainingVertex;              ///< Number of vertex of element containing (size : nLocatedpoint + 1)
  int    *_elementContainingVertex;               ///< Vertices of element containing (size : _nElementContainingVertex[nLocatedpoint])
  double *_elementContainingVertexCoords;         ///< Vertex coordinates of element containing (size : 3*_nElementContainingVertex[nLocatedpoint])
  double *_barycentricCoordinates;                ///< Barycentric coordinates (size : _nElementContainingVertex[nLocatedpoint])
  int    *_MPIrankContaining;                     ///< MPI rank containing (size : nLocatedpoint)
  cwipi_located_point_info_t located_points_info; ///< received information for located points
  LocationToLocalMesh *locationToLocalMesh;       ///< Symetric location (distant mesh points in local mesh)
};

}

#endif /* LOCATIONTODISTANTMESH_HXX_ */
