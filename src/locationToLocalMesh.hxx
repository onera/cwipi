#ifndef LOCATIONTOLOCALMESH_HXX_
#define LOCATIONTOLOCALMESH_HXX_
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

#define MPICH_IGNORE_CXX_SEEK 1

#include <string>
#include <vector>

#include <fvm_locator.h>
#include <fvm_nodal.h>
#include <bft_error.h>

#include "mesh.hxx"
#include "cwipi.h"

namespace cwipi
{

class LocationToDistantMesh;
class ApplicationProperties;

class LocationToLocalMesh
{
public:

  ///
  /// \brief Constructor
  ///
  ///   @param supportMesh                                 Mesh where distant points are localized
  ///   @param solverType                                  Solver type of current application
  ///   @param tolerance                                   Geometric tolerance for localization
  ///   @param couplingComm                                Coupling MPI communicator
  ///   @param coupledApplicationNRankCouplingComm         Rank number of coupled application
  ///   @param coupledApplicationBeginningRankCouplingComm Beginning rank of coupled application
  ///   @param isCoupledRank                               Indicate if current MPI rank is a coupled rank
  ///   @param entitiesDim                                 Geometric entities dimension
  ///   @param localApplicationProperties                  Local application properties
  ///   @param locationToDistantMesh                       Information about local points location in the distant mesh
  ///

  LocationToLocalMesh(
                      const cwipi_solver_type_t  &solverType,
                      const double &tolerance,
                      const MPI_Comm& couplingComm,
                      const int &coupledApplicationNRankCouplingComm,
                      const int &coupledApplicationBeginningRankCouplingComm,
                      const bool isCoupledRank,
                      const int entitiesDim,
                      const ApplicationProperties& localApplicationProperties,
                      LocationToDistantMesh &locationToDistantMesh);

  ///
  /// \brief Default destructor
  ///

  virtual ~LocationToLocalMesh();

  ///
  /// \brief distant points location in the local mesh
  ///

  void locate();

  ///
  /// \brief Get barycentric coordinates index  of distant points in the local mesh (size = nDistantPoint + 1)
  ///

  inline const int *getBarycentricCoordinatesIndex() const;

  ///
  /// \brief Get barycentric coordinates of distant points in the local mesh (size = BarycentricCoordinatesIndex[nDistantPoint])
  ///

  inline const double *getBarycentricCoordinates() const;

  ///
  /// \brief Return location result (size = nDistantpoint)
  ///

  inline const int *getLocation() const;

  ///
  /// \brief Return number of located distant point
  ///

  inline int getNLocatedDistantPoint() const;

  ///
  /// \brief Return number coordinates of located distant point
  ///

  inline const double *getPointCoordinates() const;

  ///
  /// \brief Return fvm locator
  ///

  inline fvm::fvm_locator_t *getFVMLocator() const;

  ///
  /// \brief Exchange field on vertices of cells that contain each located points
  ///
  ///   @param [in]  sendingField    Field defined on local mesh vertices
  ///   @param [out] receivingField  Field defined on vertices of distant
  ///                                elements that contain each located point
  ///   @param [in]  stride          Number of field component
  ///

  void exchangeCellVertexFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride);

  ///
  /// \brief Exchange field on cells that contain each located points
  ///
  ///   @param [in]  sendingField    Field defined on local mesh vertices
  ///   @param [out] receivingField  Field defined on vertices of distant
  ///                                elements that contain each located point
  ///   @param [in]  stride          Number of field component
  ///

  void exchangeCellCenterFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride);

  ///
  /// \brief Set support mesh
  ///
  ///   @param [in]      supportMesh  location support mesh
  ///

  inline void setSupportMesh(Mesh *supportMesh);

private :

  Mesh                       *_supportMesh;                                 ///< Mesh where distant points are localized
  const cwipi_solver_type_t  &_solverType;                                  ///< Solver type of current application
  const double               &_tolerance;                                   ///< Geometric tolerance for localization
  const MPI_Comm             &_couplingComm;                                ///< Coupling MPI communicator
  const int                  &_coupledApplicationNRankCouplingComm;         ///< Rank number of coupled application
  const int                  &_coupledApplicationBeginningRankCouplingComm; ///< Beginning rank of coupled application
  const bool                  _isCoupledRank;                               ///< Indicate if current MPI rank is a coupled rank
  const int                   _entitiesDim;                                 ///< Geometric entities dimension
  const ApplicationProperties&_localApplicationProperties;                  ///< Application properties
  LocationToDistantMesh      &_locationToDistantMesh;                       ///< Information about local points location in the distant mesh

  fvm::fvm_locator_t              *_fvmLocator;                                  ///< fvm structure that build the location
  int                        *_barycentricCoordinatesIndex;                 ///< Barycentric coordinates for each
  double                     *_barycentricCoordinates;                      ///< Barycentric coordinates associated to the element that contains each located distant point
  int                         _nDistantPoint;                               ///< Number of distant points located in the local mesh
  int                        *_location;                                    ///< Local elements that contain distant points
  std::vector <int>          *_nVertex;                                     ///< Vertices number of local elements that contain distant points
  bool                        _toLocate;                                    ///< Status to activate location
  int                         _maxElementContainingNVertex;                 ///< Maximum number of vertices of elements that contain located point
};

///
/// \brief Get barycentric coordinates index  of distant points in the local mesh (size = nDistantPoint + 1)
///

const int *LocationToLocalMesh::getBarycentricCoordinatesIndex() const
{
  if (_toLocate)
    bft::bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _barycentricCoordinatesIndex;
}

///
/// \brief Get barycentric coordinates of distant points in the local mesh (size = BarycentricCoordinatesIndex[nDistantPoint])
///

const double *LocationToLocalMesh::getBarycentricCoordinates() const
{
  if (_toLocate)
    bft::bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _barycentricCoordinates;
}

///
/// \brief Return location result (size = nDistantpoint)
///

const int *LocationToLocalMesh::getLocation() const
{
  if (_toLocate)
    bft::bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _location;
}

///
/// \brief Return number of located distant point
///

int LocationToLocalMesh::getNLocatedDistantPoint() const
{
  if (_toLocate)
    bft::bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _nDistantPoint;
}

///
/// \brief Return number coordinates of located distant point
///

const double *LocationToLocalMesh::getPointCoordinates() const
{
  if (_toLocate)
    bft::bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return fvm::fvm_locator_get_dist_coords(_fvmLocator);
}

///
/// \brief Return fvm locator
///

fvm::fvm_locator_t *LocationToLocalMesh::getFVMLocator() const
{
  if (_toLocate)
    bft::bft_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _fvmLocator;
}

///
/// \brief Set support mesh
///
///   @param [in]      supportMesh  location support mesh
///

inline void LocationToLocalMesh::setSupportMesh(Mesh *supportMesh)
{
  _toLocate = true;
  _supportMesh = supportMesh;
}

} // Namespace cwipi

#endif /* LOCALLOCATION_HXX_ */
