#include <cassert>
#include <cmath>

#include <iostream>

#include <mpi.h>

#include <bft_error.h>
#include <bft_printf.h>

#include <fvm_nodal_append.h>
#include <fvm_nodal_order.h>

#include "mesh.hxx"
#include "quickSort.h"

namespace couplings {


  Mesh::Mesh(const MPI_Comm &localComm,
             const int nDim,
             const int nVertex,
             const int nElts,
             const double* coords,
             int *eltConnectivityIndex,
             int *eltConnectivity
             )
    : _localComm(localComm),
      _nDim(nDim), _nVertex(nVertex), _nElts(nElts), _nPolyhedra(0), _coords(coords),
      _eltConnectivityIndex(eltConnectivityIndex), _eltConnectivity(eltConnectivity),
      _polyhedraFaceIndex(NULL), _polyhedraCellToFaceConnectivity(NULL),
      _polyhedraFaceConnectivityIndex(NULL), _polyhedraFaceConnectivity(NULL), _cellCenterCoords(NULL),
      _cellVolume(NULL), _fvmNodal(NULL), _polygonIndex(NULL)

  {
    //
    // Check dim

    if (_nDim > 3 || _nDim < 1)
      bft_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);

    //
    // Check order

    bool sorted = true;

    int nbTriangle   = 0;
    int nbQuadrangle = 0;
    int nbPoly       = 0;

    int nbTetra    = 0;
    int nbPyramid  = 0;
    int nbPrism    = 0;
    int nbHexaedra = 0;

    if (_nDim > 1) {

      if (_nDim == 2) {

        for (int i = 0; i < _nElts; i++) {
          int nCurrentEltVertex = eltConnectivityIndex[i+1] - eltConnectivityIndex[i];
          if (nCurrentEltVertex == 3) {
            if (nbQuadrangle != 0 ||
                nbPoly       != 0)
              sorted = false;
            ++nbTriangle;
          }

          else if (nCurrentEltVertex == 4) {
            if (nbPoly != 0)
              sorted = false;
            ++nbQuadrangle;
          }

          else if (nCurrentEltVertex > 4) {
            ++nbPoly;
          }

          else
            bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }
      }

      else if (_nDim == 3) {

        for (int i = 0; i < _nElts; i++) {
          int nCurrentEltVertex = eltConnectivityIndex[i+1] - eltConnectivityIndex[i];
          if (nCurrentEltVertex == 4) {
            if (nbPyramid  != 0  ||
                nbPrism    != 0  ||
                nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbTetra;
          }

          else if (nCurrentEltVertex == 5) {
            if (nbPrism    != 0  ||
                nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbPyramid;
          }

          else if (nCurrentEltVertex == 6) {
            if (nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbPrism;
          }

          else if (nCurrentEltVertex == 8) {
            if (nbPoly     != 0)
              sorted = false;
            ++nbHexaedra;
          }

          else if (nCurrentEltVertex > 8) {
            bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
            ++nbPoly;
          }

          else
            bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }
      }
    }

    //
    // Sorting

    if (!sorted) {

      switch (_nDim) {

      case 1 :
        bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Bug for edges\n");
        break;

      case 2 :
        bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : triangle, quadrangle\n");
        break;

      case 3 :
        bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : tetraedra, pyramid, prism, hexaedra\n");
        break;

      default :
        bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "unknown dimension : %i\n", _nDim);
        break;
      }
    }

    //
    //  fvm_nodal building

    _fvmNodal = fvm_nodal_create("Mesh", 3);

    switch (_nDim) {

    case 1 :
      fvm_nodal_append_shared(_fvmNodal,
                              _nElts,
                              FVM_EDGE,
                              NULL,
                              NULL,
                              NULL,
                              _eltConnectivity,
                              NULL);
      break;

    case 2 :
      if (nbTriangle != 0)
        fvm_nodal_append_shared(_fvmNodal,
                                nbTriangle,
                                FVM_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      if (nbQuadrangle != 0)
        fvm_nodal_append_shared(_fvmNodal,
                                nbQuadrangle,
                                FVM_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 3*nbTriangle,
                                NULL);

      if (nbPoly != 0) {
        _polygonIndex = new int[nbPoly+1];
        for(int i = 0; i < nbPoly+1; i++) {
          _polygonIndex[i] = _eltConnectivityIndex[nbTriangle+nbQuadrangle+i]-_eltConnectivityIndex[nbTriangle+nbQuadrangle];
        }

        fvm_nodal_append_shared(_fvmNodal,
                                nbPoly,
                                FVM_FACE_POLY,
                                NULL,
                                NULL,
                                _polygonIndex,
                                _eltConnectivity + 3*nbTriangle + 4*nbQuadrangle,
                                NULL);
      }
      break;

    case 3 :
      if (nbTetra != 0)
        fvm_nodal_append_shared(_fvmNodal,
                                nbTetra,
                                FVM_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      if (nbPyramid != 0)
        fvm_nodal_append_shared(_fvmNodal,
                                nbPyramid,
                                FVM_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra,
                                NULL);

      if (nbPrism != 0)
        fvm_nodal_append_shared(_fvmNodal,
                                nbPrism,
                                FVM_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra + 5*nbPyramid,
                                NULL);

      if (nbHexaedra != 0)
        fvm_nodal_append_shared(_fvmNodal,
                                nbHexaedra,
                                FVM_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity +  4*nbTetra + 5*nbPyramid + 6*nbPrism,
                                NULL);
      break;
    }

    //
    // Shared vertices

    fvm_nodal_set_shared_vertices(_fvmNodal, _coords);

    //
    // Order Fvm_nodal

    int localCommSize = 0;
    int *allNElts     = new int[localCommSize];
    int localRank = 0;
    unsigned int *globalEltNum = new unsigned int[_nElts];

    MPI_Comm_size(_localComm, &localCommSize);
    MPI_Comm_rank(_localComm, &localRank);

    MPI_Allgather((void *) const_cast<int*> (&_nElts),
                  1,
                  MPI_INT,
                  allNElts,
                  1,
                  MPI_INT,
                  localComm);

    int nGlobal = 0;
    for(int i = 0; i < localRank; i++)
      nGlobal += allNElts[i];

    for(int i = 0; i < nElts; i++)
      globalEltNum[i] = nGlobal + i + 1;

    switch (_nDim) {
      case 2 :
        fvm_nodal_order_faces(_fvmNodal, globalEltNum);
        break;
      case 3 :
        fvm_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
    }

    fvm_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

    //
    // global vertex num

    unsigned int *globalVertexNum = new unsigned int[_nVertex];

    MPI_Allgather((void *) const_cast<int*> (&_nVertex),
                  1,
                  MPI_INT,
                  allNElts,
                  1,
                  MPI_INT,
                  _localComm);

    nGlobal = 0;
    for(int i = 0; i < localRank; i++)
      nGlobal += allNElts[i];

    for(int i = 0; i < _nVertex; i++)
      globalVertexNum[i] = nGlobal + i + 1;

    fvm_nodal_order_vertices(_fvmNodal, globalVertexNum);
    fvm_nodal_init_io_num(_fvmNodal, globalVertexNum, 0);

    delete[] globalVertexNum;
    delete[] allNElts;

#if defined(DEBUG) && 0
    fvm_nodal_dump(_fvmNodal);
#endif

  }


  Mesh::~Mesh()
  {
    delete _cellCenterCoords;
    delete _cellVolume;
    delete[] _polygonIndex;
    fvm_nodal_destroy(_fvmNodal);
  }


  void Mesh::addPolyhedra(const int nElt,
                          int *faceIndex,
                          int *cellToFaceConnectivity,
                          int *faceConnectivityIndex,
                          int *faceConnectivity)
  {

    if (_fvmNodal == NULL)
      bft_error(__FILE__, __LINE__, 0, "No mesh to add element\n");

    _nPolyhedra += nElt;
    _nElts += nElt;

    _polyhedraFaceIndex              = faceIndex;
    _polyhedraCellToFaceConnectivity = cellToFaceConnectivity;
    _polyhedraFaceConnectivityIndex  = faceConnectivityIndex;
    _polyhedraFaceConnectivity       = faceConnectivity;

    fvm_nodal_append_shared(_fvmNodal,
                              nElt,
                              FVM_CELL_POLY,
                              faceIndex,
                              cellToFaceConnectivity,
                              faceConnectivityIndex,
                              faceConnectivity,
                              NULL);

    if (_cellCenterCoords != NULL || _cellVolume != NULL)
      _computeMeshProperties();

    int localCommSize = 0;
    int *allNElts     = new int[localCommSize];
    int localRank = 0;
    unsigned int *globalEltNum = new unsigned int[_nElts];

    MPI_Comm_size(_localComm, &localCommSize);
    MPI_Comm_rank(_localComm, &localRank);

    MPI_Allgather((void *) const_cast<int*> (&_nElts),
                  1,
                  MPI_INT,
                  allNElts,
                  1,
                  MPI_INT,
                  _localComm);

    int nGlobal = 0;
    for(int i = 0; i < localRank; i++)
      nGlobal += allNElts[i];

    for(int i = 0; i < _nElts; i++)
      globalEltNum[i] = nGlobal + i + 1;

    switch (_nDim) {
      case 3 :
        fvm_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
      default :
        bft_error(__FILE__, __LINE__,0, "Polyhedra is 3D a element !\n");
        break;
    }

    fvm_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

#if defined(DEBUG) && 0
    fvm_nodal_dump(_fvmNodal);
#endif

  }


  void Mesh::update()
  {
    if (_cellCenterCoords != NULL || _cellVolume != NULL)
      _computeMeshProperties();
  }


  void Mesh::_computeCellCenterCoordsWithVertex(const int i,
                                                const int nCurrentEltVertex,
                                                const int index,
                                                const int *eltConnectivity,
                                                std::vector<double> *cellCenterCoords)
  {
    assert (_cellCenterCoords != NULL);

    std::vector<double> &refCellCenterCoords = *cellCenterCoords;

    refCellCenterCoords[3*i]   = 0.;
    refCellCenterCoords[3*i+1] = 0.;
    refCellCenterCoords[3*i+2] = 0.;

    for (int j = 0; j < nCurrentEltVertex ; j++) {
      int ivertex = eltConnectivity[index+j] - 1;
      refCellCenterCoords[3*i]   += _coords[3*ivertex];
      refCellCenterCoords[3*i+1] += _coords[3*ivertex+1];
      refCellCenterCoords[3*i+2] += _coords[3*ivertex+2];
    }
    refCellCenterCoords[3*i]   /= nCurrentEltVertex;
    refCellCenterCoords[3*i+1] /= nCurrentEltVertex;
    refCellCenterCoords[3*i+2] /= nCurrentEltVertex;
  }


  void Mesh::_computeMeshProperties1D()
  {
    int nCurrentEltVertex;
    int index;
    std::vector<double> &refCellVolume = *_cellVolume;
    for (int i = 0; i < _nElts ; i++) {
      nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
      assert(nCurrentEltVertex == 2);
      index = _eltConnectivityIndex[i];
      _computeCellCenterCoordsWithVertex(i, nCurrentEltVertex, index,
                                         _eltConnectivity, _cellCenterCoords);
      int index = _eltConnectivityIndex[i];
      int pt1 = _eltConnectivity[index]   - 1;
      int pt2 = _eltConnectivity[index+1] - 1;
      refCellVolume[i] = sqrt((_coords[3*pt2]-_coords[3*pt1])*(_coords[3*pt2]-_coords[3*pt1])+
                              (_coords[3*pt2+1]-_coords[3*pt1+1])*(_coords[3*pt2+1]-_coords[3*pt1+1])+
                              (_coords[3*pt2+2]-_coords[3*pt1+2])*(_coords[3*pt2+2]-_coords[3*pt1+2]));
    }
  }

  void Mesh::_computeMeshProperties2D(const int  nElts,
                                      const int *faceConnectivityIndex,
                                      const int *faceConnectivity,
                                      std::vector<double> *faceNormal,
                                      std::vector<double> *faceSurface,
                                      std::vector<double> *faceCenter)

  {
    int nCurrentEltVertex;
    double v1[3];
    double v2[3];
    double barycentre[3];
    std::vector <double> triNormal;
    std::vector <double> triBarycentre;
    int index;
    std::vector<double> &refFaceSurface = *faceSurface;
    std::vector<double> &refFaceNormal  = *faceNormal;

    int maxCurrentEltVertex  = 0;
    double surftot = 0.;

    for (int i = 0; i < nElts ; i++) {
      nCurrentEltVertex = faceConnectivityIndex[i+1] - faceConnectivityIndex[i];
      if (nCurrentEltVertex > maxCurrentEltVertex)
        maxCurrentEltVertex = nCurrentEltVertex;
    }

    if (maxCurrentEltVertex > 3) {
      triNormal.resize(3*maxCurrentEltVertex, 0.);
      triBarycentre.resize(3*maxCurrentEltVertex, 0.);
    }

    for (int i = 0; i < nElts ; i++) {
      nCurrentEltVertex = faceConnectivityIndex[i+1] - faceConnectivityIndex[i];
      index = faceConnectivityIndex[i];
      _computeCellCenterCoordsWithVertex(i, nCurrentEltVertex, index,
                                         faceConnectivity,  faceCenter);
      if (nCurrentEltVertex == 3) {
        int pt1 = faceConnectivity[index]   - 1;
        int pt2 = faceConnectivity[index+1] - 1;
        int pt3 = faceConnectivity[index+2] - 1;

        v1[0] = _coords[3*pt2]   - _coords[3*pt1];
        v1[1] = _coords[3*pt2+1] - _coords[3*pt1+1];
        v1[2] = _coords[3*pt2+2] - _coords[3*pt1+2];

        v2[0] = _coords[3*pt3]   - _coords[3*pt1];
        v2[1] = _coords[3*pt3+1] - _coords[3*pt1+1];
        v2[2] = _coords[3*pt3+2] - _coords[3*pt1+2];

        refFaceNormal[3*i]   = 0.5 * (v1[1]*v2[2]-v1[2]*v2[1]);
        refFaceNormal[3*i+1] = 0.5 * (v1[2]*v2[0]-v1[0]*v2[2]);
        refFaceNormal[3*i+2] = 0.5 * (v1[0]*v2[1]-v1[1]*v2[0]);

        refFaceSurface[i] = sqrt(refFaceNormal[3*i]*refFaceNormal[3*i]+
                                 refFaceNormal[3*i+1]*refFaceNormal[3*i+1]+
                                 refFaceNormal[3*i+2]*refFaceNormal[3*i+2]);
      }

      else if (nCurrentEltVertex > 3) {
        std::vector<double> &refFaceCenter = *faceCenter;
        barycentre[0] = refFaceCenter[3*i];
        barycentre[1] = refFaceCenter[3*i+1];
        barycentre[2] = refFaceCenter[3*i+2];

        refFaceNormal[3*i]   = 0;
        refFaceNormal[3*i+1] = 0;
        refFaceNormal[3*i+2] = 0;

        refFaceCenter[3*i] = 0.;
        refFaceCenter[3*i+1] = 0.;
        refFaceCenter[3*i+2] = 0.;

        for (int k = 0; k < nCurrentEltVertex ; k++) {
          int pt1 = faceConnectivity[index + k] - 1;
          int pt2 = faceConnectivity[index + (k+1)%nCurrentEltVertex] - 1;
          v1[0] = _coords[3*pt1]   - barycentre[0];
          v1[1] = _coords[3*pt1+1] - barycentre[1];
          v1[2] = _coords[3*pt1+2] - barycentre[2];

          v2[0] = _coords[3*pt2]   - barycentre[0];
          v2[1] = _coords[3*pt2+1] - barycentre[1];
          v2[2] = _coords[3*pt2+2] - barycentre[2];

          triBarycentre[3*k]   = (_coords[3*pt1]   + _coords[3*pt2]   + barycentre[0]) /3.;
          triBarycentre[3*k+1] = (_coords[3*pt1+1] + _coords[3*pt2+1] + barycentre[1]) /3.;
          triBarycentre[3*k+2] = (_coords[3*pt1+2] + _coords[3*pt2+2] + barycentre[2]) /3.;

          triNormal[3*k]   = 0.5 * (v1[1]*v2[2]-v1[2]*v2[1]);
          triNormal[3*k+1] = 0.5 * (v1[2]*v2[0]-v1[0]*v2[2]);
          triNormal[3*k+2] = 0.5 * (v1[0]*v2[1]-v1[1]*v2[0]);

          refFaceNormal[3*i]   += triNormal[3*k];
          refFaceNormal[3*i+1] += triNormal[3*k+1];
          refFaceNormal[3*i+2] += triNormal[3*k+2];
        }

        refFaceSurface[i] = 0;
        for (int k = 0; k < nCurrentEltVertex ; k++) {

          double triSurface = sqrt(triNormal[3*k]   * triNormal[3*k]+
                                   triNormal[3*k+1] * triNormal[3*k+1]+
                                   triNormal[3*k+2] * triNormal[3*k+2]);

          if (triNormal[3*k]   * refFaceNormal[0] +
              triNormal[3*k+1] * refFaceNormal[1] +
              triNormal[3*k+2] * refFaceNormal[2] < 0.)

            triSurface *= -1;
          refFaceSurface[i] += triSurface;

          refFaceCenter[3*i]   += triSurface * triBarycentre[3*k];
          refFaceCenter[3*i+1] += triSurface * triBarycentre[3*k+1];
          refFaceCenter[3*i+2] += triSurface * triBarycentre[3*k+2];

        }

        refFaceCenter[3*i]   /= refFaceSurface[i];
        refFaceCenter[3*i+1] /= refFaceSurface[i];
        refFaceCenter[3*i+2] /= refFaceSurface[i];

      }
      surftot += refFaceSurface[i];
    }
  }

  void Mesh::_computeMeshProperties3D()
  {
    std::vector<double> refCellCenterCoords = *_cellCenterCoords;
    std::vector<double> refCellVolume = *_cellVolume;
    std::vector<double> tetraCenterCoords(3,0.);


    int nStandardElement = _nElts - _nPolyhedra;
    if (nStandardElement > 0) {
      std::vector<int>  faceConnectivityIndex(4,0);
      std::vector<int>  faceConnectivity(12,0);
      std::vector<int>  tetraConnec(24,0);
      int ntetra =0;
      std::vector<double> faceSurface(4,0.);
      std::vector<double> faceCenter(3*4,0.);
      std::vector<double> faceNormal(3*4,0.);

      int nFace;

      for (int i = 0; i < nStandardElement ; i++) {
        int nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
        int index = _eltConnectivityIndex[i];

        //
        // Element spliting

        switch(nCurrentEltVertex) {

        case 4 :
          ntetra = 1;
          tetraConnec[0] = _eltConnectivity[index];
          tetraConnec[1] = _eltConnectivity[index+1];
          tetraConnec[2] = _eltConnectivity[index+2];
          tetraConnec[3] = _eltConnectivity[index+3];
          break;

        case 5 :

          ntetra = 2;
          tetraConnec[0] = _eltConnectivity[index];
          tetraConnec[1] = _eltConnectivity[index+1];
          tetraConnec[2] = _eltConnectivity[index+3];
          tetraConnec[3] = _eltConnectivity[index+4];
          tetraConnec[4] = _eltConnectivity[index+1];
          tetraConnec[5] = _eltConnectivity[index+2];
          tetraConnec[6] = _eltConnectivity[index+3];
          tetraConnec[7] = _eltConnectivity[index+4];
          break;

        case 6 :

          ntetra = 3;
          tetraConnec[0] = _eltConnectivity[index];
          tetraConnec[1] = _eltConnectivity[index+1];
          tetraConnec[2] = _eltConnectivity[index+2];
          tetraConnec[3] = _eltConnectivity[index+3];

          tetraConnec[4] = _eltConnectivity[index+3];
          tetraConnec[5] = _eltConnectivity[index+5];
          tetraConnec[6] = _eltConnectivity[index+4];
          tetraConnec[7] = _eltConnectivity[index+1];

          tetraConnec[8] = _eltConnectivity[index+2];
          tetraConnec[9] = _eltConnectivity[index+5];
          tetraConnec[10] = _eltConnectivity[index+3];
          tetraConnec[11] = _eltConnectivity[index+1];
          break;

        case 8 :

          ntetra = 6;
          tetraConnec[0] = _eltConnectivity[index];
          tetraConnec[1] = _eltConnectivity[index+1];
          tetraConnec[2] = _eltConnectivity[index+3];
          tetraConnec[3] = _eltConnectivity[index+4];

          tetraConnec[4] = _eltConnectivity[index+5];
          tetraConnec[5] = _eltConnectivity[index+4];
          tetraConnec[6] = _eltConnectivity[index+7];
          tetraConnec[7] = _eltConnectivity[index+1];

          tetraConnec[8] = _eltConnectivity[index+4];
          tetraConnec[9] = _eltConnectivity[index+3];
          tetraConnec[10] = _eltConnectivity[index+7];
          tetraConnec[11] = _eltConnectivity[index+1];

          tetraConnec[12] = _eltConnectivity[index+1];
          tetraConnec[13] = _eltConnectivity[index+2];
          tetraConnec[14] = _eltConnectivity[index+3];
          tetraConnec[15] = _eltConnectivity[index+6];

          tetraConnec[16] = _eltConnectivity[index+5];
          tetraConnec[17] = _eltConnectivity[index+7];
          tetraConnec[18] = _eltConnectivity[index+6];
          tetraConnec[19] = _eltConnectivity[index+1];

          tetraConnec[20] = _eltConnectivity[index+1];
          tetraConnec[21] = _eltConnectivity[index+3];
          tetraConnec[22] = _eltConnectivity[index+7];
          tetraConnec[23] = _eltConnectivity[index+6];
          break;
        }

        refCellCenterCoords[3*i] = 0.;
        refCellCenterCoords[3*i+1] = 0.;
        refCellCenterCoords[3*i+2] = 0.;
        refCellVolume[i] = 0.;

        for (int j = 0; j < ntetra; j++) {
          double tetraVolume = 0.;

          // Tetraedre case
          // local faces : 1, 3, 2,
          //               1, 2, 4
          //               1, 4, 3
          //               2, 3, 4

          nFace = 4;
          faceConnectivityIndex[1] = 3;
          faceConnectivityIndex[2] = 6;
          faceConnectivityIndex[3] = 9;
          faceConnectivityIndex[4] = 12;
          faceConnectivity[0]      = tetraConnec[4*j+0];
          faceConnectivity[1]      = tetraConnec[4*j+2];
          faceConnectivity[2]      = tetraConnec[4*j+1];
          faceConnectivity[3]      = tetraConnec[4*j+0];
          faceConnectivity[4]      = tetraConnec[4*j+1];
          faceConnectivity[5]      = tetraConnec[4*j+3];
          faceConnectivity[6]      = tetraConnec[4*j+0];
          faceConnectivity[7]      = tetraConnec[4*j+3];
          faceConnectivity[8]      = tetraConnec[4*j+2];
          faceConnectivity[9]      = tetraConnec[4*j+1];
          faceConnectivity[10]     = tetraConnec[4*j+2];
          faceConnectivity[11]     = tetraConnec[4*j+3];

          _computeMeshProperties2D(nFace,
                                   &faceConnectivityIndex[0],
                                   &faceConnectivity[0],
                                   &faceNormal,
                                   &faceSurface,
                                   &faceCenter);

          double cellSurface = 0.;
          for (int k = 0; k < nFace; k++) {
            cellSurface += faceSurface[k];
            tetraCenterCoords[0] += faceNormal[k]   * faceCenter[3*k];
            tetraCenterCoords[1] += faceNormal[k+1] * faceCenter[3*k+1];
            tetraCenterCoords[2] += faceNormal[k+2] * faceCenter[3*k+2];
            tetraVolume += faceNormal[k] * faceCenter[3*k] +
              faceNormal[k+1] * faceCenter[3*k+1] +
                         faceNormal[k+2] * faceCenter[3*k+2];
          }
          tetraVolume *= 1./3.;
          tetraCenterCoords[0] /= cellSurface;
          tetraCenterCoords[1] /= cellSurface;
          tetraCenterCoords[2] /= cellSurface;

          refCellCenterCoords[3*i]   += tetraVolume*tetraCenterCoords[0] ;
          refCellCenterCoords[3*i+1] += tetraVolume*tetraCenterCoords[1];
          refCellCenterCoords[3*i+2] += tetraVolume*tetraCenterCoords[2];
          refCellVolume[i] += tetraVolume;
        }
        refCellCenterCoords[3*i] /= refCellVolume[i];
        refCellCenterCoords[3*i+1] /= refCellVolume[i];
        refCellCenterCoords[3*i+2] /= refCellVolume[i];
      }
    }

    if (_nPolyhedra > 0) {

      //
      // Polyedra splitting

      // Not yet implemented
      bft_error(__FILE__, __LINE__, 0, "Not implemented yet\n");

//       fvm_tesselation_t *fvm_tesselation = fvm_tesselation_create(FVM_CELL_POLY,
//                                                                   _nPolyhedra,
//                                                                   _polyhedraFaceIndex,
//                                                                   _polyhedraCellToFaceConnectivity,
//                                                                   _polyhedraFaceConnectivityIndex,
//                                                                   _polyhedraFaceConnectivity,
//                                                                   null);

//       fvm_tesselation_init(fvm_tesselation, 3,_coords, NULL, NULL);


//       std::vector<int>  tesselationFaceConnectivityIndex;
//       std::vector<int>  tesselationFaceConnectivity;
//       std::vector<int>  triangulateFaceConnectivityIndex;
//       std::vector<int>  triangulateFaceConnectivity;
//       std::vector<int>  faceConnectivityIndex;
//       std::vector<int>  faceConnectivity;
//       int ntetra = 0;
//       std::vector<double> faceSurface(4,0.);
//       std::vector<double> faceCenter(3*4,0.);
//       std::vector<double> faceNormal(3*4,0.);

//       int maxFacePolyhedra = 0;
//       int maxVertexFace = 0;

//       for (int i = 0; i < _nPolyhedra; i++) {

//         int nFacePolyhedra = _polyhedraFaceIndex[i+1] - _polyhedraFaceIndex[i];
//         if (maxFacePolyhedra < nFacePolyhedra)
//           maxFacePolyhedra = nFacePolyhedra;
//         int faceIndex = _polyhedraCellToFaceConnectivity[i];
//         for (int j = 0; j < nFacePolyhedra; j++) {
//           int iface = _polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
//           int nVertexFace += (_polyhedraFaceConnectivityIndex[iface+1] - _polyhedraFaceConnectivityIndex[iface]);
//           if (maxVertexFace < nVertexFace)
//             maxVertexFace = nVertexFace;
//         }
//       }


//       std::vector<int>  triangulateFaceConnectivityIndex(maxFacePolyhedra*(maxVertexFace-2) + 1);
//       std::vector<int>  triangulateFaceConnectivity(3*maxFacePolyhedra*(maxVertexFace-2));
//       std::vector<int>  faceConnectivity(maxVertexFace, 0);

//       for (int i = 0; i < _nPolyhedra; i++) {

//         int nFacePolyhedra = _polyhedraFaceIndex[i+1] - _polyhedraFaceIndex[i];
//         int faceIndex = _polyhedraCellToFaceConnectivity[i];

//         for (int j = 0; j < nFacePolyhedra; j++) {
//           int iface = _polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
//           int nVertexFace += (_polyhedraFaceConnectivityIndex[iface+1] - _polyhedraFaceConnectivityIndex[iface]);
//           if (maxVertexFace < nVertexFace)
//             maxVertexFace = nVertexFace;
//         }
//       }


//       if (maxFacePolyhedra > faceSurface.size()) {
//         faceSurface.resize(maxFacePolyhedra,0.);
//         faceCenter.resize(3*maxFacePolyhedra,0.);
//         faceNormal.resize(3*maxFacePolyhedra,0.);
//         faceConnectivityIndex.resize(maxFacePolyhedra+1);
//       }

//       if  (maxVertexFace > faceConnectivity.size())
//         faceConnectivity.resize(maxVertexFace);


//       fvm_triangulate_state_t* state = fvm_triangulate_state_create(maxVertexFace);

//       std::vector<int> triangulateConnectivityIndex = NULL;
//       std::vector<int> triangulateConnectivity = NULL;

//       for (int i = 0; i < _nPolyhedra; i++) {
//         int nFacePolyhedra = _polyhedraFaceIndex[i+1] - _polyhedraFaceIndex[i];
//         int faceIndex = _polyhedraCellToFaceConnectivity[i];
//         int nVertexFace = 0;
//         faceConnectivityIndex[0] = 0;
//         for (int j = 0; j < nFacePolyhedra; j++) {
//           int iface = _polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
//           bool reorient = false;
//           if (iface < 0) {
//             iface = -iface;
//             reorient = true;
//           }

//           nVertexFace += _polyhedraFaceConnectivityIndex[iface+1] - _polyhedraFaceConnectivityIndex[iface];
//           int vertexIndex = _polyhedraFaceConnectivityIndex[iface];
//           faceConnectivityIndex[j+1] = faceConnectivityIndex[j] + nVertexFace;
//           for (int k = 0; k < nVertexFace; k++) {
//             if (reorient)
//               faceConnectivity[faceConnectivityIndex[j+1]-1-k] = _polyhedraFaceConnectivity[vertexIndex+k];
//             else
//               faceConnectivity[faceConnectivityIndex[j]+k] = _polyhedraFaceConnectivity[vertexIndex+k];
//           }

//           fvm_triangulate_polygon(3,
//                                   nVertexFace,
//                                   _coords,
//                                   NULL,
//                                   faceConnectivity,
//                                   FVM_TRIANGULATE_MESH_DEF,
//                                   triangulateConnectivity[triangulateConnectivityIndex[j]],
//                                   state);

//         }

//         _computeMeshProperties2D(nFacePolyhedra,
//                                  &faceConnectivityIndex[0],
//                                  &faceConnectivity[0],
//                                  &faceNormal,
//                                  &faceSurface,
//                                  &faceCenter);

//         double cellSurface = 0.;

//         refCellCenterCoords[3*(nStandardElement+i)]   = 0.;
//         refCellCenterCoords[3*(nStandardElement+i)+1] = 0.;
//         refCellCenterCoords[3*(nStandardElement+i)+2] = 0.;
//         refCellVolume[nStandardElement+i] = 0.;

//         for (int j = 0; j < nFace; j++) {
//           cellSurface += faceNormal[i];
//           refCellCenterCoords[3*(nStandardElement+i)]   += faceNormal[i] * faceCenter[3*i];
//           refCellCenterCoords[3*(nStandardElement+i)+1] += faceNormal[i+1] * faceCenter[3*i+1];
//           refCellCenterCoords[3*(nStandardElement+i)+2] += faceNormal[i+2] * faceCenter[3*i+2];
//           refCellVolume[(nStandardElement+i)] += faceNormal[i] * faceCenter[3*i] +
//                                                faceNormal[i+1] * faceCenter[3*i+1] +
//                                                faceNormal[i+2] * faceCenter[3*i+2];
//         }
//         refCellVolume[nStandardElement+i] *= 1./3.;
//         refCellCenterCoords[3*(nStandardElement+i)]   /= cellSurface;
//         refCellCenterCoords[3*(nStandardElement+i)+1] /= cellSurface;
//         refCellCenterCoords[3*(nStandardElement+i)+2] /= cellSurface;
//       }
//       fvm_triangulate_state_destroy(state);
    }
  }

  void Mesh::_computeMeshProperties()
  {
    if (_cellCenterCoords == NULL)
      _cellCenterCoords = new std::vector<double>(3*_nElts);
    else if (_cellCenterCoords->size() < 3*_nElts)
      _cellCenterCoords->resize(3*_nElts);

    if (_cellVolume == NULL)
      _cellVolume = new std::vector<double>(_nElts);
    else if (_cellVolume->size() < _nElts)
      _cellVolume->resize(_nElts);

    if (_nDim == 1)
        _computeMeshProperties1D();
    else if (_nDim == 2) {
      std::vector<double> faceNormal(3*_nElts);
      _computeMeshProperties2D(_nElts,
                               _eltConnectivityIndex,
                               _eltConnectivity,
                               &faceNormal,
                               _cellVolume,
                               _cellCenterCoords);
    }
    else if (_nDim == 3)
      _computeMeshProperties3D();

  }

}
