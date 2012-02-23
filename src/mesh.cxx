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

#include <cassert>
#include <cmath>

#include <iostream>

#include <bftc_error.h>
#include <bftc_printf.h>

#include <fvmc_nodal_append.h>
#include <fvmc_nodal_order.h>
#include <fvmc_parall.h>

#include "mesh.hxx"
#include "quickSort.h"

#define ABS(a)     ((a) <  0  ? -(a) : (a))
//TODO: Remplacer les #define MAX, ABS... par des fonctions inline
namespace cwipi {

  Mesh::Mesh(const MPI_Comm &localComm,
             const int nDim,
             const int nVertex,
             const int nElts,
             double* coords,
             int *eltConnectivityIndex,
             int *eltConnectivity
             )
    : _localComm(localComm),
      _nDim(nDim), _nVertex(nVertex), _nElts(nElts), _nPolyhedra(0), _coords(coords),
      _eltConnectivityIndex(eltConnectivityIndex), _eltConnectivity(eltConnectivity),
      _polyhedraFaceIndex(NULL), _polyhedraCellToFaceConnectivity(NULL),
      _polyhedraFaceConnectivityIndex(NULL), _polyhedraFaceConnectivity(NULL), 
      _polyhedraCellToVertexConnectivity(NULL),
      _polyhedraCellToVertexConnectivityIndex(NULL), _cellCenterCoords(NULL),
      _cellVolume(NULL), _normalFace(NULL),_fvmNodal(NULL), 
      _polygonIndex(NULL), _isNodalFinalized(false),
      _characteristicLength(NULL),
      _isDegenerated(NULL)
      
  {

    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvmc_parall_set_mpi_comm(localComm);

    //
    // Check dim

    if (_nDim > 3 || _nDim < 1)
      bftc_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);

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
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
          
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
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
            ++nbPoly;
          }

          else
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }
      }
    }

    //
    // Sorting

    if (!sorted) {

      switch (_nDim) {

      case 1 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Bug for edges\n");
        break;

      case 2 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : triangle, quadrangle\n");
        break;

      case 3 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : tetraedra, pyramid, prism, hexaedra\n");
        break;

      default :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "unknown dimension : %i\n", _nDim);
        break;
      }
    }

    //
    // fvmc_nodal building

    _fvmNodal = fvmc_nodal_create("Mesh", 3);

    //
    // Sections building

    switch (_nDim) {

    case 1 :


      fvmc_nodal_append_shared(_fvmNodal,
                              _nElts,
                              FVMC_EDGE,
                              NULL,
                              NULL,
                              NULL,
                              _eltConnectivity,
                              NULL);
      break;

    case 2 :

      int nTriangleSum;
      int nQuadrangleSum;
      int nPolySum;

      MPI_Allreduce (&nbTriangle, &nTriangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbQuadrangle, &nQuadrangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPoly, &nPolySum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      if (nbTriangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTriangle,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTriangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbQuadrangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbQuadrangle,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 3*nbTriangle,
                                NULL);


      else if (nQuadrangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPoly != 0) {
        _polygonIndex = new int[nbPoly+1];
        for(int i = 0; i < nbPoly+1; i++) {
          _polygonIndex[i] = _eltConnectivityIndex[nbTriangle+nbQuadrangle+i]-_eltConnectivityIndex[nbTriangle+nbQuadrangle];
        }

        fvmc_nodal_append_shared(_fvmNodal,
                                nbPoly,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                _polygonIndex,
                                _eltConnectivity + 3*nbTriangle + 4*nbQuadrangle,
                                NULL);
      }

      else if (nPolySum != 0) {

        //bftc_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for a empty polygon section\n");

        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);
      }
    break;

    case 3 :
      int nTetraSum;
      int nPyramidSum;
      int nPrismSum;
      int nHexaedraSum;

      MPI_Allreduce (&nbTetra, &nTetraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPyramid, &nPyramidSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPrism, &nPrismSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbHexaedra, &nHexaedraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      if (nbTetra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTetra,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTetraSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPyramid != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPyramid,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra,
                                NULL);

      else if (nPyramidSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPrism != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPrism,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra + 5*nbPyramid,
                                NULL);

      else if (nPrismSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbHexaedra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbHexaedra,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity +  4*nbTetra + 5*nbPyramid + 6*nbPrism,
                                NULL);

      else if (nbHexaedra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);


      break;
    }

     MPI_Barrier(localComm);
     fvmc_parall_set_mpi_comm(oldFVMComm);

  }

  void Mesh::_finalizeNodal()
  {

    if (!_isNodalFinalized) {

      MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
      if (oldFVMComm != MPI_COMM_NULL)
        MPI_Barrier(oldFVMComm);
      fvmc_parall_set_mpi_comm(_localComm);

      _isNodalFinalized = true;

      //
      // Shared vertices
      
      fvmc_nodal_set_shared_vertices(_fvmNodal, _coords);
      
      //
      // Order Fvm_nodal
      
      int localCommSize = 0;
      int localRank = 0;
      
      MPI_Comm_size(_localComm, &localCommSize);
      MPI_Comm_rank(_localComm, &localRank);
      
      int *allNElts     = new int[localCommSize];
      unsigned int *globalEltNum = new unsigned int[_nElts];
      
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
      case 2 :
        fvmc_nodal_order_faces(_fvmNodal, globalEltNum);
        break;
      case 3 :
        fvmc_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
      }

      fvmc_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

      delete [] globalEltNum;

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
      
      
      fvmc_nodal_order_vertices(_fvmNodal, globalVertexNum);
      fvmc_nodal_init_io_num(_fvmNodal, globalVertexNum, 0);
      
      delete[] globalVertexNum;
      delete[] allNElts;
      
#if defined(DEBUG) && 0
      fvmc_nodal_dump(_fvmNodal);
#endif
      
      MPI_Barrier(_localComm);
      fvmc_parall_set_mpi_comm(oldFVMComm);
    }
  }


  /////////
  /////////////
  ///////////
  /////////

  Mesh::Mesh(const MPI_Comm &localComm,
             fvmc_nodal_t* fvmc_nodal)
    : _localComm(localComm),
      _nDim(fvmc_nodal_get_dim(fvmc_nodal)), _nVertex(0),
      _nElts(0), _nPolyhedra(0), _coords(NULL),
      _eltConnectivityIndex(NULL), _eltConnectivity(NULL),
      _polyhedraFaceIndex(NULL), _polyhedraCellToFaceConnectivity(NULL),
      _polyhedraFaceConnectivityIndex(NULL), _polyhedraFaceConnectivity(NULL), _polyhedraCellToVertexConnectivity(NULL), 
      _polyhedraCellToVertexConnectivityIndex(NULL), _cellCenterCoords(NULL),
      _cellVolume(NULL), _normalFace(NULL),_fvmNodal(NULL), _polygonIndex(NULL), _isNodalFinalized(true)


  {
    //
    // Copy
    
    // TODO: Attention dans ce cas on alloue les vecteurs en interne sans les liberer
    //       Restructurer en creant des const * pour les valeurs partagees et détruite les valeurs 
    //       creees dans le destructeur



    fvmc_nodal_get_vertex(fvmc_nodal,
                         &_nElts,
                         &_eltConnectivityIndex,
                         &_eltConnectivity);


    fvmc_nodal_get_coords(fvmc_nodal,
                         &_nVertex,
                         &_coords);


    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvmc_nodal_destroy(fvmc_nodal); 

    fvmc_parall_set_mpi_comm(localComm);

    //
    // Check dim

    if (_nDim > 3 || _nDim < 1)
      bftc_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);

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
          int nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
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
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }
      }

      else if (_nDim == 3) {

        for (int i = 0; i < _nElts; i++) {
          int nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
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
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
            ++nbPoly;
          }

          else
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }

      }
    }

    //
    // Sorting

    if (!sorted) {

      switch (_nDim) {

      case 1 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Bug for edges\n");
        break;

      case 2 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : triangle, quadrangle\n");
        break;

      case 3 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : tetraedra, pyramid, prism, hexaedra\n");
        break;

      default :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "unknown dimension : %i\n", _nDim);
        break;
      }
    }

    //
    // fvmc_nodal building

    _fvmNodal = fvmc_nodal_create("Mesh", 3);

    //
    // Sections building

    switch (_nDim) {

    case 1 :


      fvmc_nodal_append_shared(_fvmNodal,
                              _nElts,
                              FVMC_EDGE,
                              NULL,
                              NULL,
                              NULL,
                              _eltConnectivity,
                              NULL);
      break;

    case 2 :

      int nTriangleSum;
      int nQuadrangleSum;
      int nPolySum;

      MPI_Allreduce (&nbTriangle, &nTriangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbQuadrangle, &nQuadrangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPoly, &nPolySum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      if (nbTriangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTriangle,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTriangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbQuadrangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbQuadrangle,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 3*nbTriangle,
                                NULL);


      else if (nQuadrangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPoly != 0) {
        _polygonIndex = new int[nbPoly+1];
        for(int i = 0; i < nbPoly+1; i++) {
          _polygonIndex[i] = _eltConnectivityIndex[nbTriangle+nbQuadrangle+i]-_eltConnectivityIndex[nbTriangle+nbQuadrangle];
        }

        fvmc_nodal_append_shared(_fvmNodal,
                                nbPoly,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                _polygonIndex,
                                _eltConnectivity + 3*nbTriangle + 4*nbQuadrangle,
                                NULL);
      }

      else if (nPolySum != 0) {

        //bftc_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for a empty polygon section\n");

        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);
      }
    break;

    case 3 :
      int nTetraSum;
      int nPyramidSum;
      int nPrismSum;
      int nHexaedraSum;

      MPI_Allreduce (&nbTetra, &nTetraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPyramid, &nPyramidSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPrism, &nPrismSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbHexaedra, &nHexaedraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      if (nbTetra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTetra,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTetraSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPyramid != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPyramid,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra,
                                NULL);

      else if (nPyramidSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPrism != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPrism,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra + 5*nbPyramid,
                                NULL);

      else if (nPrismSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbHexaedra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbHexaedra,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity +  4*nbTetra + 5*nbPyramid + 6*nbPrism,
                                NULL);

      else if (nbHexaedra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);


      break;
    }

    //
    // Shared vertices

    fvmc_nodal_set_shared_vertices(_fvmNodal, _coords);

    //
    // Order Fvm_nodal

    int localCommSize = 0;
    int localRank = 0;

    MPI_Comm_size(_localComm, &localCommSize);
    MPI_Comm_rank(_localComm, &localRank);

    int *allNElts     = new int[localCommSize];
    unsigned int *globalEltNum = new unsigned int[_nElts];

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
      case 2 :
        fvmc_nodal_order_faces(_fvmNodal, globalEltNum);
        break;
      case 3 :
        fvmc_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
    }

    fvmc_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

    delete [] globalEltNum;

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


    fvmc_nodal_order_vertices(_fvmNodal, globalVertexNum);
    fvmc_nodal_init_io_num(_fvmNodal, globalVertexNum, 0);

    delete[] globalVertexNum;
    delete[] allNElts;

    #if defined(DEBUG) && 0
    fvmc_nodal_dump(_fvmNodal);
    #endif

    MPI_Barrier(localComm);
    fvmc_parall_set_mpi_comm(oldFVMComm);

  }








  Mesh::~Mesh()
  {
    delete _cellCenterCoords;
    delete _cellVolume;
    delete _normalFace;
    delete[] _polygonIndex;
    fvmc_nodal_destroy(_fvmNodal);
  }


  void Mesh::addPolyhedra(const int nElt,
                          int *faceIndex,
                          int *cellToFaceConnectivity,
                          int *faceConnectivityIndex,
                          int *faceConnectivity)
  {
    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvmc_parall_set_mpi_comm(_localComm);

    if (_fvmNodal == NULL)
      bftc_error(__FILE__, __LINE__, 0, "No mesh to add element\n");

    _nPolyhedra += nElt;
    _nElts += nElt;

    _polyhedraFaceIndex              = faceIndex;
    _polyhedraCellToFaceConnectivity = cellToFaceConnectivity;
    _polyhedraFaceConnectivityIndex  = faceConnectivityIndex;
    _polyhedraFaceConnectivity       = faceConnectivity;

    if (nElt > 0)

      fvmc_nodal_append_shared(_fvmNodal,
                              nElt,
                              FVMC_CELL_POLY,
                              faceIndex,
                              cellToFaceConnectivity,
                              faceConnectivityIndex,
                              faceConnectivity,
                              NULL);
    else {
      //bftc_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for an empty polyedron section\n");


      fvmc_nodal_append_shared(_fvmNodal,
                                 0,
                                 FVMC_CELL_POLY,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL);
    }

    if (_cellCenterCoords != NULL || _cellVolume != NULL)
      _computeMeshProperties();
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
    std::vector<double> &refCharacteristicLength = *_characteristicLength;
    std::vector<bool>   &refIsDegenerated        = *_isDegenerated;

    for (int i = 0; i < _nElts ; i++) {
      refIsDegenerated[i] = false;
      nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
      assert(nCurrentEltVertex == 2);
      index = _eltConnectivityIndex[i];
      _computeCellCenterCoordsWithVertex(i, nCurrentEltVertex, index,
                                         _eltConnectivity, _cellCenterCoords);
      int index = _eltConnectivityIndex[i];
      int pt1 = _eltConnectivity[index]   - 1;
      int pt2 = _eltConnectivity[index+1] - 1;
      refCellVolume[i] = 
        sqrt((_coords[3*pt2]-_coords[3*pt1])*(_coords[3*pt2]-_coords[3*pt1])
            +(_coords[3*pt2+1]-_coords[3*pt1+1])*(_coords[3*pt2+1]-_coords[3*pt1+1]) 
            +(_coords[3*pt2+2]-_coords[3*pt1+2])*(_coords[3*pt2+2]-_coords[3*pt1+2]));
      refCharacteristicLength[i] = refCellVolume[i];
      if (refCharacteristicLength[i] < GEOM_EPS_DIST) {
        refIsDegenerated[i] = true;
        bftc_printf("Warning computeMeshProperties : linear element '%i' is degenerated "
                   "(distance between 2 vertices = %12.5e)\n",
                   i,
                   refCharacteristicLength[i]);
      }
    }
  }

  void Mesh::_computeMeshProperties2D(const int  nElts,
                                      const int *faceConnectivityIndex,
                                      const int *faceConnectivity,
                                      std::vector<double> *faceNormal,
                                      std::vector<double> *characteristicLength,
                                      std::vector<bool>   *isDegenerated,
                                      std::vector<double> *faceSurface,
                                      std::vector<double> *faceCenter)

  {
    const double big = 1e32; 
    int nCurrentEltVertex;
    double v1[3];
    double v2[3];
    double v3[3];
    double barycentre[3];
    std::vector <double> triNormal;
    std::vector <double> triBarycentre;
    int index;
    std::vector<double> &refFaceSurface = *faceSurface;
    std::vector<double> &refFaceNormal  = *faceNormal;
    std::vector<double> &refCharacteristicLength = *characteristicLength;
    std::vector<bool>   &refIsDegenerated        = *isDegenerated;

    int maxCurrentEltVertex  = 0;
    double surftot = 0.;

    for (int i = 0; i < nElts ; i++) {
      nCurrentEltVertex = faceConnectivityIndex[i+1] - faceConnectivityIndex[i];
      refCharacteristicLength[i] = big;
      refIsDegenerated[i] = false;
      if (nCurrentEltVertex > maxCurrentEltVertex)
        maxCurrentEltVertex = nCurrentEltVertex;
    }

    triNormal.resize(3*maxCurrentEltVertex, 0.);
    triBarycentre.resize(3*maxCurrentEltVertex, 0.);

    for (int i = 0; i < nElts ; i++) {
      nCurrentEltVertex = faceConnectivityIndex[i+1] - faceConnectivityIndex[i];
      index = faceConnectivityIndex[i];
      _computeCellCenterCoordsWithVertex(i, nCurrentEltVertex, index,
                                         faceConnectivity,  faceCenter);
      if (nCurrentEltVertex == 3) {
        const int pt1 = faceConnectivity[index]   - 1;
        const int pt2 = faceConnectivity[index+1] - 1;
        const int pt3 = faceConnectivity[index+2] - 1;
        double d1;
        double d2;
        double d3;

        v1[0] = _coords[3*pt2]   - _coords[3*pt1];
        v1[1] = _coords[3*pt2+1] - _coords[3*pt1+1];
        v1[2] = _coords[3*pt2+2] - _coords[3*pt1+2];

        v2[0] = _coords[3*pt3]   - _coords[3*pt1];
        v2[1] = _coords[3*pt3+1] - _coords[3*pt1+1];
        v2[2] = _coords[3*pt3+2] - _coords[3*pt1+2];

        v3[0] = _coords[3*pt3]   - _coords[3*pt2];
        v3[1] = _coords[3*pt3+1] - _coords[3*pt2+1];
        v3[2] = _coords[3*pt3+2] - _coords[3*pt2+2];

        d1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
        d2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
        d3 = sqrt(v3[0]*v2[0] + v3[1]*v3[1] + v3[2]*v3[2]);

        refCharacteristicLength[i] = std::min(d1, d2);
        refCharacteristicLength[i] = std::min(d3, refCharacteristicLength[i]);

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

        refCharacteristicLength[i] = big;

        for (int k = 0; k < nCurrentEltVertex ; k++) {
          int pt1 = faceConnectivity[index + k] - 1;
          int pt2 = faceConnectivity[index + (k+1)%nCurrentEltVertex] - 1;
          double xx = _coords[3*pt1    ] - _coords[3*pt2    ];
          double yy = _coords[3*pt1 + 1] - _coords[3*pt2 + 1];
          double zz = _coords[3*pt1 + 2] - _coords[3*pt2 + 2];
          double lEdge = sqrt(xx*xx + yy*yy + zz*zz);

          refCharacteristicLength[i] = std::min(refCharacteristicLength[i], lEdge); 

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

          if (triNormal[3*k]   * refFaceNormal[3*i] +
              triNormal[3*k+1] * refFaceNormal[3*i+1] +
              triNormal[3*k+2] * refFaceNormal[3*i+2] < 0.)

            triSurface *= -1;
          refFaceSurface[i] += triSurface;

          refFaceCenter[3*i]   += triSurface * triBarycentre[3*k];
          refFaceCenter[3*i+1] += triSurface * triBarycentre[3*k+1];
          refFaceCenter[3*i+2] += triSurface * triBarycentre[3*k+2];

        }

        double eps_loc = computeGeometricEpsilon(refCharacteristicLength[i], GEOM_EPS_SURF);

        if (refFaceSurface[i] <= eps_loc) {
          refIsDegenerated[i] = true;
          bftc_printf("Warning computeMeshProperties : surface element '%i' is degenerated (surface = %12.5e)\n", 
                     i, 
                     refFaceSurface[i]);
          refFaceCenter[3*i] = barycentre[0];
          refFaceCenter[3*i+1] = barycentre[1];
          refFaceCenter[3*i+2] = barycentre[2];
        }
        else {
          refFaceCenter[3*i]   /= refFaceSurface[i];
          refFaceCenter[3*i+1] /= refFaceSurface[i];
          refFaceCenter[3*i+2] /= refFaceSurface[i];
        }

#ifdef INFINITY
        if (isinf(ABS(refFaceCenter[3*i]))   ||
            isinf(ABS(refFaceCenter[3*i+1])) ||
            isinf(ABS(refFaceCenter[3*i+2]))) {
          refIsDegenerated[i] = true;
          bftc_printf("Warning computeMeshProperties : surface element '%i' is degenerated "
                     "(center coordinates = (%12.5e, %12.5e, %12.5e))\n", 
                     i, 
                     refFaceCenter[3*i],
                     refFaceCenter[3*i+1],
                     refFaceCenter[3*i+2]);
          refFaceCenter[3*i] = barycentre[0];
          refFaceCenter[3*i+1] = barycentre[1];
          refFaceCenter[3*i+2] = barycentre[2];
        }
#endif

      }
      refFaceSurface[i] = ABS(refFaceSurface[i]);
      surftot += refFaceSurface[i];
    }
  }

  void Mesh::_computeMeshProperties3D()
  {

    // TODO: Refaire completement les calcul des centres cellules et volume (comme dans CEDRE)

    std::vector<double> & refCellCenterCoords = *_cellCenterCoords;
    std::vector<double> & refCellVolume = *_cellVolume;
    std::vector<double> &refCharacteristicLength = *_characteristicLength;
    std::vector<bool>   &refIsDegenerated        = *_isDegenerated;
    std::vector<double> barycentre(3, 0.);

    int nStandardElement = _nElts - _nPolyhedra;
    if (nStandardElement > 0) {
      std::vector<int>  faceConnectivityIndex(5,0);
      std::vector<int>  faceConnectivity(12,0);
      std::vector<int>  tetraConnec(24,0);
      int ntetra =0;
      std::vector<double> faceSurface(4,0.);
      std::vector<double> faceCharacteristicLength(4,0.);
      std::vector<bool>   faceIsDegenerated(4,false);
      std::vector<double> faceCenter(3*4,0.);
      std::vector<double> faceNormal(3*4,0.);

      int nFace;

      for (int i = 0; i < nStandardElement ; i++) {
        int nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
        int index = _eltConnectivityIndex[i];

        refIsDegenerated[i] = false;

        //
        // Element spliting

        for (int j = 0; j < barycentre.size(); j++)
          barycentre[j] = 0.;

        switch(nCurrentEltVertex) {

        case 4 :
          ntetra = 1;
          tetraConnec[0] = _eltConnectivity[index];
          tetraConnec[1] = _eltConnectivity[index+1];
          tetraConnec[2] = _eltConnectivity[index+2];
          tetraConnec[3] = _eltConnectivity[index+3];
          for (int j = 0; j < 4 ; j++) {
            barycentre[0] += _coords[3 * tetraConnec[j]    ]; 
            barycentre[1] += _coords[3 * tetraConnec[j] + 1]; 
            barycentre[2] += _coords[3 * tetraConnec[j] + 2]; 
          }
          barycentre[0] / 4;
          barycentre[1] / 4;
          barycentre[2] / 4;
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
          for (int j = 0; j < 4 ; j++) {
            barycentre[0] += _coords[3 * tetraConnec[j]    ]; 
            barycentre[1] += _coords[3 * tetraConnec[j] + 1]; 
            barycentre[2] += _coords[3 * tetraConnec[j] + 2]; 
          }
          barycentre[0] += _coords[3 * tetraConnec[5]    ]; 
          barycentre[1] += _coords[3 * tetraConnec[5] + 1]; 
          barycentre[2] += _coords[3 * tetraConnec[5] + 2]; 
          barycentre[0] / 5;
          barycentre[1] / 5;
          barycentre[2] / 5;
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
          for (int j = 0; j < 4 ; j++) {
            barycentre[0] += _coords[3 * tetraConnec[j]    ]; 
            barycentre[1] += _coords[3 * tetraConnec[j] + 1]; 
            barycentre[2] += _coords[3 * tetraConnec[j] + 2]; 
          }
          barycentre[0] += _coords[3 * tetraConnec[5]    ]; 
          barycentre[1] += _coords[3 * tetraConnec[5] + 1]; 
          barycentre[2] += _coords[3 * tetraConnec[5] + 2]; 
          barycentre[0] += _coords[3 * tetraConnec[6]    ]; 
          barycentre[1] += _coords[3 * tetraConnec[6] + 1]; 
          barycentre[2] += _coords[3 * tetraConnec[6] + 2]; 
          barycentre[0] / 6;
          barycentre[1] / 6;
          barycentre[2] / 6;
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
          for (int j = 0; j < 4 ; j++) {
            barycentre[0] += _coords[3 * tetraConnec[j]    ]; 
            barycentre[1] += _coords[3 * tetraConnec[j] + 1]; 
            barycentre[2] += _coords[3 * tetraConnec[j] + 2]; 
          }
          barycentre[0] += _coords[3 * tetraConnec[5]    ]; 
          barycentre[1] += _coords[3 * tetraConnec[5] + 1]; 
          barycentre[2] += _coords[3 * tetraConnec[5] + 2]; 

          barycentre[0] += _coords[3 * tetraConnec[13]    ]; 
          barycentre[1] += _coords[3 * tetraConnec[13] + 1]; 
          barycentre[2] += _coords[3 * tetraConnec[13] + 2];
 
          barycentre[0] += _coords[3 * tetraConnec[15]    ]; 
          barycentre[1] += _coords[3 * tetraConnec[15] + 1]; 
          barycentre[2] += _coords[3 * tetraConnec[15] + 2];

          barycentre[0] += _coords[3 * tetraConnec[6]    ]; 
          barycentre[1] += _coords[3 * tetraConnec[6] + 1]; 
          barycentre[2] += _coords[3 * tetraConnec[6] + 2];

          barycentre[0] / 8;
          barycentre[1] / 8;
          barycentre[2] / 8;
          break;
        }

        refCellCenterCoords[3*i] = 0.;
        refCellCenterCoords[3*i+1] = 0.;
        refCellCenterCoords[3*i+2] = 0.;
        refCellVolume[i] = 0.;

        for (int j = 0; j < ntetra; j++) {
          std::vector<double> tetraCenterCoords(3,0.);
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
                                   &faceCharacteristicLength,
                                   &faceIsDegenerated,
                                   &faceSurface,
                                   &faceCenter);

          for (int k = 0; k < 4; k++) {
            tetraCenterCoords[0] += _coords[3*(tetraConnec[4*j+k]-1)];
            tetraCenterCoords[1] += _coords[3*(tetraConnec[4*j+k]-1)+1];
            tetraCenterCoords[2] += _coords[3*(tetraConnec[4*j+k]-1)+2];
          }
          tetraCenterCoords[0] /= 4;
          tetraCenterCoords[1] /= 4;
          tetraCenterCoords[2] /= 4;

          double cellSurface = 0.;
          for (int k = 0; k < nFace; k++) {
            if (faceIsDegenerated[k]) {
              bftc_printf("Warning : a face of an volume element '%i'"
                         " is degenerated (surface face '%i' = %12.5e)\n", 
                         i, faceSurface[k]); 
            }
            cellSurface += faceSurface[k];
            tetraVolume +=  faceNormal[3*k]   * faceCenter[3*k] +
                            faceNormal[3*k+1] * faceCenter[3*k+1] +
                            faceNormal[3*k+2] * faceCenter[3*k+2];
          }

          tetraVolume *= 1./3.;

          refCellCenterCoords[3*i]   += tetraVolume*tetraCenterCoords[0] ;
          refCellCenterCoords[3*i+1] += tetraVolume*tetraCenterCoords[1];
          refCellCenterCoords[3*i+2] += tetraVolume*tetraCenterCoords[2];
          refCellVolume[i] += tetraVolume;
        }

        //
        // Check element
        // 

        double eps_loc = computeGeometricEpsilon(refCharacteristicLength[i], GEOM_EPS_VOL);

        if (refCellVolume[i] <= eps_loc) {
          refIsDegenerated[i] = true;
          bftc_printf("Warning computeMeshProperties : volume element '%i' is degenerated (volume = %12.5e)\n", 
                     i, 
                     refCellVolume[i]);
          refCellCenterCoords[3*i]   = barycentre[1];
          refCellCenterCoords[3*i+1] = barycentre[2];
          refCellCenterCoords[3*i+2] = barycentre[3];
        }
        else {
          refCellCenterCoords[3*i]   /= refCellVolume[i];
          refCellCenterCoords[3*i+1] /= refCellVolume[i];
          refCellCenterCoords[3*i+2] /= refCellVolume[i];
        }
      }
    }

    if (_nPolyhedra > 0) {

      //
      // Polyedra splitting

      // Not yet implemented
      bftc_error(__FILE__, __LINE__, 0, "Not implemented yet\n");

//       fvmc_tesselation_t *fvmc_tesselation = fvmc_tesselation_create(FVMC_CELL_POLY,
//                                                                   _nPolyhedra,
//                                                                   _polyhedraFaceIndex,
//                                                                   _polyhedraCellToFaceConnectivity,
//                                                                   _polyhedraFaceConnectivityIndex,
//                                                                   _polyhedraFaceConnectivity,
//                                                                   null);

//       fvmc_tesselation_init(fvmc_tesselation, 3,_coords, NULL, NULL);


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

//      refIsDegenerated[i] = false;

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


//       fvmc_triangulate_state_t* state = fvmc_triangulate_state_create(maxVertexFace);

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

//           fvmc_triangulate_polygon(3,
//                                   nVertexFace,
//                                   _coords,
//                                   NULL,
//                                   faceConnectivity,
//                                   FVMC_TRIANGULATE_MESH_DEF,
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
//       fvmc_triangulate_state_destroy(state);
    }
  }

  void Mesh::_computeMeshProperties()
  {
    if (_cellCenterCoords == NULL) {
      _cellCenterCoords = new std::vector<double>(3*_nElts);
      _cellVolume = new std::vector<double>(_nElts);
      _characteristicLength = new std::vector<double>(_nElts);
      _isDegenerated = new std::vector<bool>(_nElts);
    }
    else {
      _cellCenterCoords->resize(3*_nElts);
      _cellVolume->resize(_nElts);
      _characteristicLength->resize(_nElts);
      _isDegenerated = new std::vector<bool>(_nElts);
    }

    if (_nDim == 1)
        _computeMeshProperties1D();
    else if (_nDim == 2) {
      if (_cellCenterCoords == NULL)
        _normalFace = new std::vector<double>(3*_nElts);
      else
        _normalFace->resize(3*_nElts);
      _computeMeshProperties2D(_nElts,
                               _eltConnectivityIndex,
                               _eltConnectivity,
                               _normalFace,
                               _characteristicLength,
                               _isDegenerated,
                               _cellVolume,
                               _cellCenterCoords);
    }
    else if (_nDim == 3)
      _computeMeshProperties3D();

  }
  
  void Mesh::_computeMeshPolyhedraProperties(){
    
    int nbrVertex = 0;
    int nbrVertexOld = 0;
    int sizeTabLoc = 12;
    int indFaceVertex;
    int indVertex;
    int nbrVertexFace;
    
    std::vector<int> vertexPolyLoc (sizeTabLoc);
    std::vector<int> vertexPolyBool(_nVertex);

    if(_polyhedraCellToVertexConnectivity == NULL)
       _polyhedraCellToVertexConnectivity = new std::vector<int>(0);
    
    if (_polyhedraCellToVertexConnectivityIndex == NULL)
      _polyhedraCellToVertexConnectivityIndex = new std::vector<int>(0);

    std::vector<int> & refPolyhedraCellToVertexConnectivity = *_polyhedraCellToVertexConnectivity;
    std::vector<int> & refPolyhedraCellToVertexConnectivityIndex = *_polyhedraCellToVertexConnectivityIndex;

    for (int i = 0; i < _nVertex; i++)
      vertexPolyBool[i] = 0;    
    
    /**** Premiere boucle pour connaitre le nombre total de sommets ****/

    for(int iPoly = 0; iPoly < _nPolyhedra ; iPoly++){

       int nFacePolyhedra = _polyhedraFaceIndex[iPoly+1] - _polyhedraFaceIndex[iPoly];

       for(int iFace = 0; iFace < nFacePolyhedra ; iFace++){

         nbrVertexFace = _polyhedraFaceConnectivityIndex[_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly] + iFace]]
                       - _polyhedraFaceConnectivityIndex[_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly] + iFace] - 1];
         
         indFaceVertex = _polyhedraFaceConnectivityIndex[_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly]  + iFace] - 1];

         for (int iVertex = 0 ; iVertex < nbrVertexFace ; iVertex++){

           indVertex = _polyhedraFaceConnectivity[indFaceVertex + iVertex];

           if(vertexPolyBool[indVertex - 1] == 0){

             if(nbrVertex > sizeTabLoc){
               sizeTabLoc *=2;
               vertexPolyLoc.resize(sizeTabLoc);
             }

             vertexPolyLoc[nbrVertex - nbrVertexOld] = indVertex;
             nbrVertex++;
             vertexPolyBool[indVertex - 1] = 1;       
           }           

         }

       }

       for(int iVertexPoly = 0 ; iVertexPoly < nbrVertex - nbrVertexOld ; iVertexPoly++)
         vertexPolyBool[vertexPolyLoc[iVertexPoly] - 1] = 0;
 
       nbrVertexOld = nbrVertex;
    }

    vertexPolyLoc.clear();

    refPolyhedraCellToVertexConnectivity.resize(nbrVertex);
    refPolyhedraCellToVertexConnectivityIndex.resize(_nPolyhedra+1);
    refPolyhedraCellToVertexConnectivityIndex[0] = 0;

    nbrVertex = 0;
    nbrVertexOld = 0;

    /**** Calcul des tableaux d'index et de connectivite entre poly et sommets ****/

    for(int iPoly = 0; iPoly < _nPolyhedra ; iPoly++){

       int nFacePolyhedra = _polyhedraFaceIndex[iPoly+1] - _polyhedraFaceIndex[iPoly];

      for(int iFace = 0; iFace < nFacePolyhedra ; iFace++){

        nbrVertexFace = _polyhedraFaceConnectivityIndex[_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly] + iFace ] ]
                      - _polyhedraFaceConnectivityIndex[_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly] + iFace ] - 1];
         
        indFaceVertex = _polyhedraFaceConnectivityIndex[_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly] + iFace] - 1];        

        for (int iVertex = 0 ; iVertex < nbrVertexFace ; iVertex++){

          indVertex = _polyhedraFaceConnectivity[indFaceVertex + iVertex];

          if(vertexPolyBool[indVertex - 1] == 0){
            refPolyhedraCellToVertexConnectivity[nbrVertex] = indVertex;            
            nbrVertex++;
            vertexPolyBool[indVertex - 1] = 1;                   
          }           
          
        }
        
      }
      
      for(int iVertexPoly = nbrVertexOld ; iVertexPoly < nbrVertex ; iVertexPoly++)
        vertexPolyBool[refPolyhedraCellToVertexConnectivity[iVertexPoly] - 1] = 0;
      
      refPolyhedraCellToVertexConnectivityIndex[iPoly+1] = nbrVertex;
      
      nbrVertexOld = nbrVertex;

    }    
  }      

}
