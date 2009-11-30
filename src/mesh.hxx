#ifndef __COUPLING_MESH_H__
#define __COUPLING_MESH_H__

#include <fvm_nodal.h>
#include <vector>
#include <mpi.h>

namespace cwipi {

  class Mesh {

  public:
    Mesh(const MPI_Comm &localComm,
         const int nDim,
         const int nVertex,
         const int nElts,
         const double* coords,
         int *eltConnectivityIndex,
         int *eltConnectivity);

    virtual ~Mesh();

    void addPolyhedra(const int nElt,
                      int *faceIndex,
                      int *cellToFaceConnectivity,
                      int *faceConnectivityIndex,
                      int *faceConnectivity);

    inline const int& getNVertex() const;

    inline const double* getVertexCoords() const;

    inline fvm_nodal_t& getFvmNodal() const;

    inline const int& getNElts() const;

    inline const int& getNPolyhedra() const;

    inline const int* getEltConnectivityIndex() const;

    inline const int* getEltConnectivity() const;

    inline const std::vector<double>& getVolume();

    inline const std::vector<double>& getCellCenterCoords();

    inline const int *getPolyhedraFaceIndex() const;

    inline const int *getPolyhedraCellToFaceConnectivity() const;

    inline const int *getPolyhedraFaceConnectivityIndex() const;

    inline const int *getPolyhedraFaceConnectivity() const;

    void update();

  private:
    Mesh();

    Mesh(const Mesh&);

    Mesh& operator=(const Mesh&);

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

  private:
    // TODO: renommer _nDim par entitesDim
    const MPI_Comm & _localComm;
    const int     _nDim;
    const int     _nVertex;
    int           _nElts;
    int           _nPolyhedra;
    const double *_coords;
    int          *_eltConnectivityIndex;
    int          *_polygonIndex;
    int          *_eltConnectivity;
    int          *_polyhedraFaceIndex;
    int          *_polyhedraCellToFaceConnectivity;
    int          *_polyhedraFaceConnectivityIndex;
    int          *_polyhedraFaceConnectivity;
    std::vector<double>  *_cellCenterCoords;
    std::vector<double>  *_cellVolume;
    fvm_nodal_t *_fvmNodal;
  };

  const int& Mesh::getNVertex()  const
  {
    return _nVertex;
  }

  const double* Mesh::getVertexCoords()  const
  {
    return _coords;
  }

  fvm_nodal_t& Mesh::getFvmNodal() const
  {
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

}


#endif //__COUPLING_MESH_H__
