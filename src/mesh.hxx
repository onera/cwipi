#ifndef __COUPLING_MESH_H__
#define __COUPLING_MESH_H__
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

#include <fvm_nodal.h>
#include <vector>
#include <mpi.h>
#include <bft_error.h>

namespace cwipi {

  class Mesh {

  public:
    Mesh(const MPI_Comm &localComm,
         const int nDim,
         const int nVertex,
         const int nElts,
         double* coords,
         int *eltConnectivityIndex,
         int *eltConnectivity);

    Mesh(const MPI_Comm &localComm,
         fvm::fvm_nodal_t* fvm_nodal);

    virtual ~Mesh();

    void addPolyhedra(const int nElt,
                      int *faceIndex,
                      int *cellToFaceConnectivity,
                      int *faceConnectivityIndex,
                      int *faceConnectivity);

    inline const int& getNVertex() const;

    inline const double* getVertexCoords() const;

    inline fvm::fvm_nodal_t& getFvmNodal() const;

    inline fvm::fvm_nodal_t& getFvmNodal();

    inline const int& getNElts() const;

    inline const int& getNPolyhedra() const;

    inline const int* getEltConnectivityIndex() const;

    inline const int* getEltConnectivity() const;

    inline const std::vector<double>& getVolume();

    inline const std::vector<double>& getCellCenterCoords();

    inline const std::vector<double>& getNormalFace();    

    inline const int *getPolyhedraFaceIndex() const;

    inline const int *getPolyhedraCellToFaceConnectivity() const;

    inline const int *getPolyhedraFaceConnectivityIndex() const;

    inline const int *getPolyhedraFaceConnectivity() const;

    inline const std::vector<int>& getPolyhedraCellToVertexConnectivity();

    inline const std::vector<int>& getPolyhedraCellToVertexConnectivityIndex();

    void update();

  private :
    Mesh();

    Mesh(const Mesh&);

    Mesh& operator=(const Mesh&);

  protected :

    void _computeCellCenterCoordsWithVertex(const int i,
                                            const int nCurrentEltVertex,
                                            const int index,
                                            const int *eltConnectivity,
                                            std::vector<double> *cellCenterCoords) ;

    void _computeMeshProperties1D() ;

    void _computeMeshProperties2D(const int  nElts,
                                  const int *faceConnectivityIndex,
                                  const int *faceConnectivity,
                                  std::vector<double> *faceNormal,
                                  std::vector<double> *faceSurface,
                                  std::vector<double> *faceCenter);

    void _computeMeshProperties3D();

    void _computeMeshProperties();

    void _finalizeNodal();

    void _computeMeshPolyhedraProperties();

  private:
    // TODO: renommer _nDim par entitesDim
    const MPI_Comm & _localComm;
    int            _nDim;
    int           _nVertex;
    int           _nElts;
    int           _nPolyhedra;
    double       *_coords;
    int          *_polygonIndex;

    int          *_eltConnectivityIndex;
    int          *_eltConnectivity;

    int          *_polyhedraFaceIndex;
    int          *_polyhedraCellToFaceConnectivity;
    int          *_polyhedraFaceConnectivityIndex;
    int          *_polyhedraFaceConnectivity;
    bool         _isNodalFinalized;
    std::vector<int>          *_polyhedraCellToVertexConnectivity;
    std::vector<int>          *_polyhedraCellToVertexConnectivityIndex;

    std::vector<double>  *_cellCenterCoords;
    std::vector<double>  *_cellVolume;
    std::vector<double>  *_normalFace;

    fvm::fvm_nodal_t *_fvmNodal;
  };

  const int& Mesh::getNVertex()  const
  {
    return _nVertex;
  }

  const double* Mesh::getVertexCoords()  const
  {
    return _coords;
  }

  fvm::fvm_nodal_t& Mesh::getFvmNodal()
  {
    _finalizeNodal();
    return *_fvmNodal;
  }

  fvm::fvm_nodal_t& Mesh::getFvmNodal() const 
  {
    if (! _isNodalFinalized) {
      bft::bft_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);
    }
    return *_fvmNodal;
  }

  const int& Mesh::getNElts() const
  {
    return _nElts;
  }

  const int& Mesh::getNPolyhedra() const
  {
    return _nPolyhedra;
  }

  const int* Mesh::getEltConnectivityIndex() const
  {
    return _eltConnectivityIndex;
  }

  const int* Mesh::getEltConnectivity() const
  {
    return _eltConnectivity;
  }

  const std::vector<double>& Mesh::getVolume()
  {
    if (_cellVolume == NULL)
      _computeMeshProperties();
    return *_cellVolume;
  }

  const std::vector<double>& Mesh::getCellCenterCoords()
  {
    if (_cellCenterCoords == NULL)
      _computeMeshProperties();

    return *_cellCenterCoords;
  }

  const std::vector<double>& Mesh::getNormalFace()
  {
    if (_normalFace == NULL)
      _computeMeshProperties();

    return *_normalFace;
  }
  

  const int *Mesh::getPolyhedraFaceIndex() const
  {
    return _polyhedraFaceIndex;
  }

  const int *Mesh::getPolyhedraCellToFaceConnectivity() const
  {
    return _polyhedraCellToFaceConnectivity;
  }

  const int *Mesh::getPolyhedraFaceConnectivityIndex() const
  {
    return _polyhedraFaceConnectivityIndex;
  }

  const int *Mesh::getPolyhedraFaceConnectivity() const
  {
    return _polyhedraFaceConnectivity;
  }

  const std::vector<int>& Mesh::getPolyhedraCellToVertexConnectivity()
  {
    if (_polyhedraCellToVertexConnectivity == NULL)
      _computeMeshPolyhedraProperties();
    
    return *_polyhedraCellToVertexConnectivity;
  }
  
  const std::vector<int>& Mesh::getPolyhedraCellToVertexConnectivityIndex()
  {
    if (_polyhedraCellToVertexConnectivityIndex == NULL)
      _computeMeshPolyhedraProperties();
    
    return *_polyhedraCellToVertexConnectivityIndex;
  }


}

#endif //__COUPLING_MESH_H__
