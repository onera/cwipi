#ifndef __COUPLING_MESH_H__
#define __COUPLING_MESH_H__

#include <fvm_nodal.h>
#include <vector>

namespace couplings {

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

    inline const int& getNVertex();

    inline const double* getVertexCoords();

    inline fvm_nodal_t& getFvmNodal();

    inline const int& getNElts();

    inline const int& getNPolyhedra();

    inline const int* getEltConnectivityIndex();

    inline const int* getEltConnectivity();

    inline const std::vector<double>& getVolume();

    inline const std::vector<double>& getCellCenterCoords();

    inline const int *getPolyhedraFaceIndex();

    inline const int *getPolyhedraCellToFaceConnectivity();

    inline const int *getPolyhedraFaceConnectivityIndex();

    inline const int *getPolyhedraFaceConnectivity();

    void update();

  private:
    Mesh();

    Mesh(const Mesh&);

    Mesh& operator=(const Mesh&);

    void _computeCellCenterCoordsWithVertex(const int i,
                                            const int nCurrentEltVertex,
                                            const int index,
                                            const int *eltConnectivity,
                                            std::vector<double> *cellCenterCoords);

    void _computeMeshProperties1D();

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

  const int& Mesh::getNVertex()
  {
    return _nVertex;
  }

  const double* Mesh::getVertexCoords()
  {
    return _coords;
  }

  fvm_nodal_t& Mesh::getFvmNodal()
  {
    return *_fvmNodal;
  }

  const int& Mesh::getNElts()
  {
    return _nElts;
  }

  const int& Mesh::getNPolyhedra()
  {
    return _nPolyhedra;
  }

  const int* Mesh::getEltConnectivityIndex()
  {
    return _eltConnectivityIndex;
  }

  const int* Mesh::getEltConnectivity()
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

  const int *Mesh::getPolyhedraFaceIndex()
  {
    return _polyhedraFaceIndex;
  }

  const int *Mesh::getPolyhedraCellToFaceConnectivity()
  {
    return _polyhedraCellToFaceConnectivity;
  }

  const int *Mesh::getPolyhedraFaceConnectivityIndex()
  {
    return _polyhedraFaceConnectivityIndex;
  }

  const int *Mesh::getPolyhedraFaceConnectivity()
  {
    return _polyhedraFaceConnectivity;
  }

}


#endif //__COUPLING_MESH_H__
