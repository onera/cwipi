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
//Bug mpich2
//#define MPICH_IGNORE_CXX_SEEK 1

#include <mpi.h>

#include <cassert>
#include <cmath>

#include <iostream>

#include <bft_error.h>
#include <bft_printf.h>

#include <fvm_nodal_append.h>
#include <fvm_nodal_order.h>
#include <fvm_parall.h>

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
      _polyhedraFaceConnectivityIndex(NULL), _polyhedraFaceConnectivity(NULL), _cellCenterCoords(NULL),
      _cellVolume(NULL), _fvmNodal(NULL), _polygonIndex(NULL)

  {

    MPI_Comm oldFVMComm = fvm::fvm_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvm::fvm_parall_set_mpi_comm(localComm);

    //
    // Check dim

    if (_nDim > 3 || _nDim < 1)
      bft::bft_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);

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
            bft::bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
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
            bft::bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
            ++nbPoly;
          }

          else
            bft::bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }
      }
    }

    //
    // Sorting

    if (!sorted) {

      switch (_nDim) {

      case 1 :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Bug for edges\n");
        break;

      case 2 :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : triangle, quadrangle\n");
        break;

      case 3 :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : tetraedra, pyramid, prism, hexaedra\n");
        break;

      default :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "unknown dimension : %i\n", _nDim);
        break;
      }
    }

    //
    // fvm_nodal building

    _fvmNodal = fvm::fvm_nodal_create("Mesh", 3);

    //
    // Sections building

    switch (_nDim) {

    case 1 :


      fvm::fvm_nodal_append_shared(_fvmNodal,
                              _nElts,
                              fvm::FVM_EDGE,
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
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbTriangle,
                                fvm::FVM_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTriangleSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbQuadrangle != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbQuadrangle,
                                fvm::FVM_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 3*nbTriangle,
                                NULL);


      else if (nQuadrangleSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_FACE_QUAD,
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

        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbPoly,
                                fvm::FVM_FACE_POLY,
                                NULL,
                                NULL,
                                _polygonIndex,
                                _eltConnectivity + 3*nbTriangle + 4*nbQuadrangle,
                                NULL);
      }

      else if (nPolySum != 0) {

        //bft::bft_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for a empty polygon section\n");

        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_FACE_POLY,
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
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbTetra,
                                fvm::FVM_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTetraSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPyramid != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbPyramid,
                                fvm::FVM_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra,
                                NULL);

      else if (nPyramidSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPrism != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbPrism,
                                fvm::FVM_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra + 5*nbPyramid,
                                NULL);

      else if (nPrismSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbHexaedra != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbHexaedra,
                                fvm::FVM_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity +  4*nbTetra + 5*nbPyramid + 6*nbPrism,
                                NULL);

      else if (nbHexaedra != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);


      break;
    }

    //
    // Shared vertices

    fvm::fvm_nodal_set_shared_vertices(_fvmNodal, _coords);

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

    for(int i = 0; i < nElts; i++)
      globalEltNum[i] = nGlobal + i + 1;

    switch (_nDim) {
      case 2 :
        fvm::fvm_nodal_order_faces(_fvmNodal, globalEltNum);
        break;
      case 3 :
        fvm::fvm_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
    }

    fvm::fvm_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

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


    fvm::fvm_nodal_order_vertices(_fvmNodal, globalVertexNum);
    fvm::fvm_nodal_init_io_num(_fvmNodal, globalVertexNum, 0);

    delete[] globalVertexNum;
    delete[] allNElts;

    #if defined(DEBUG) && 0
    fvm::fvm_nodal_dump(_fvmNodal);
    #endif

    MPI_Barrier(localComm);
    fvm::fvm_parall_set_mpi_comm(oldFVMComm);

  }



  /////////
  /////////////
  ///////////
  /////////

  Mesh::Mesh(const MPI_Comm &localComm,
             fvm::fvm_nodal_t* fvm_nodal)
    : _localComm(localComm),
      _nDim(fvm::fvm_nodal_get_dim(fvm_nodal)), _nVertex(0),
      _nElts(0), _nPolyhedra(0), _coords(NULL),
      _eltConnectivityIndex(NULL), _eltConnectivity(NULL),
      _polyhedraFaceIndex(NULL), _polyhedraCellToFaceConnectivity(NULL),
      _polyhedraFaceConnectivityIndex(NULL), _polyhedraFaceConnectivity(NULL), _cellCenterCoords(NULL),
      _cellVolume(NULL), _fvmNodal(NULL), _polygonIndex(NULL)

  {
    //
    // Copy
    
    // TODO: Attention dans ce cas on alloue les vecteurs en interne sans les liberer
    //       Restructurer en creant des const * pour les valeurs partagees et détruite les valeurs 
    //       creees dans le destructeur



    fvm::fvm_nodal_get_vertex(fvm_nodal,
                         &_nElts,
                         &_eltConnectivityIndex,
                         &_eltConnectivity);


    fvm::fvm_nodal_get_coords(fvm_nodal,
                         &_nVertex,
                         &_coords);


    MPI_Comm oldFVMComm = fvm::fvm_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvm::fvm_nodal_destroy(fvm_nodal); 

    fvm::fvm_parall_set_mpi_comm(localComm);

    //
    // Check dim

    if (_nDim > 3 || _nDim < 1)
      bft::bft_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);

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
            bft::bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
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
            bft::bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
            ++nbPoly;
          }

          else
            bft::bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }

      }
    }

    //
    // Sorting

    if (!sorted) {

      switch (_nDim) {

      case 1 :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Bug for edges\n");
        break;

      case 2 :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : triangle, quadrangle\n");
        break;

      case 3 :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : tetraedra, pyramid, prism, hexaedra\n");
        break;

      default :
        bft::bft_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "unknown dimension : %i\n", _nDim);
        break;
      }
    }

    //
    // fvm_nodal building

    _fvmNodal = fvm::fvm_nodal_create("Mesh", 3);

    //
    // Sections building

    switch (_nDim) {

    case 1 :


      fvm::fvm_nodal_append_shared(_fvmNodal,
                              _nElts,
                              fvm::FVM_EDGE,
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
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbTriangle,
                                fvm::FVM_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTriangleSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbQuadrangle != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbQuadrangle,
                                fvm::FVM_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 3*nbTriangle,
                                NULL);


      else if (nQuadrangleSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_FACE_QUAD,
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

        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbPoly,
                                fvm::FVM_FACE_POLY,
                                NULL,
                                NULL,
                                _polygonIndex,
                                _eltConnectivity + 3*nbTriangle + 4*nbQuadrangle,
                                NULL);
      }

      else if (nPolySum != 0) {

        //bft::bft_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for a empty polygon section\n");

        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_FACE_POLY,
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
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbTetra,
                                fvm::FVM_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTetraSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPyramid != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbPyramid,
                                fvm::FVM_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra,
                                NULL);

      else if (nPyramidSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPrism != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbPrism,
                                fvm::FVM_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra + 5*nbPyramid,
                                NULL);

      else if (nPrismSum != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbHexaedra != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                nbHexaedra,
                                fvm::FVM_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity +  4*nbTetra + 5*nbPyramid + 6*nbPrism,
                                NULL);

      else if (nbHexaedra != 0)
        fvm::fvm_nodal_append_shared(_fvmNodal,
                                0,
                                fvm::FVM_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);


      break;
    }

    //
    // Shared vertices

    fvm::fvm_nodal_set_shared_vertices(_fvmNodal, _coords);

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
        fvm::fvm_nodal_order_faces(_fvmNodal, globalEltNum);
        break;
      case 3 :
        fvm::fvm_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
    }

    fvm::fvm_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

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


    fvm::fvm_nodal_order_vertices(_fvmNodal, globalVertexNum);
    fvm::fvm_nodal_init_io_num(_fvmNodal, globalVertexNum, 0);

    delete[] globalVertexNum;
    delete[] allNElts;

    #if defined(DEBUG) && 0
    fvm::fvm_nodal_dump(_fvmNodal);
    #endif

    MPI_Barrier(localComm);
    fvm::fvm_parall_set_mpi_comm(oldFVMComm);

  }








  Mesh::~Mesh()
  {
    delete _cellCenterCoords;
    delete _cellVolume;
    delete[] _polygonIndex;
    fvm::fvm_nodal_destroy(_fvmNodal);
  }


  void Mesh::addPolyhedra(const int nElt,
                          int *faceIndex,
                          int *cellToFaceConnectivity,
                          int *faceConnectivityIndex,
                          int *faceConnectivity)
  {
    MPI_Comm oldFVMComm = fvm::fvm_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvm::fvm_parall_set_mpi_comm(_localComm);

    if (_fvmNodal == NULL)
      bft::bft_error(__FILE__, __LINE__, 0, "No mesh to add element\n");

    _nPolyhedra += nElt;
    _nElts += nElt;

    _polyhedraFaceIndex              = faceIndex;
    _polyhedraCellToFaceConnectivity = cellToFaceConnectivity;
    _polyhedraFaceConnectivityIndex  = faceConnectivityIndex;
    _polyhedraFaceConnectivity       = faceConnectivity;

    if (nElt > 0)

      fvm::fvm_nodal_append_shared(_fvmNodal,
                              nElt,
                              fvm::FVM_CELL_POLY,
                              faceIndex,
                              cellToFaceConnectivity,
                              faceConnectivityIndex,
                              faceConnectivity,
                              NULL);
    else {
      //bft::bft_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for an empty polyedron section\n");


      fvm::fvm_nodal_append_shared(_fvmNodal,
                                 0,
                                 fvm::FVM_CELL_POLY,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL);
    }
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
        fvm::fvm_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
      default :
        bft::bft_error(__FILE__, __LINE__,0, "Polyhedra is 3D a element !\n");
        break;
    }

    fvm::fvm_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

    delete [] globalEltNum;

#if defined(DEBUG) && 0
    fvm::fvm_nodal_dump(_fvmNodal);
#endif

    MPI_Barrier(_localComm);
    fvm::fvm_parall_set_mpi_comm(oldFVMComm);
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

    triNormal.resize(3*maxCurrentEltVertex, 0.);
    triBarycentre.resize(3*maxCurrentEltVertex, 0.);

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

          if (triNormal[3*k]   * refFaceNormal[3*i] +
              triNormal[3*k+1] * refFaceNormal[3*i+1] +
              triNormal[3*k+2] * refFaceNormal[3*i+2] < 0.)

            triSurface *= -1;
          refFaceSurface[i] += triSurface;

          refFaceCenter[3*i]   += triSurface * triBarycentre[3*k];
          refFaceCenter[3*i+1] += triSurface * triBarycentre[3*k+1];
          refFaceCenter[3*i+2] += triSurface * triBarycentre[3*k+2];

        }

        refFaceCenter[3*i]   /= refFaceSurface[i];
        refFaceCenter[3*i+1] /= refFaceSurface[i];
        refFaceCenter[3*i+2] /= refFaceSurface[i];
#ifdef NAN
        if (isnan(refFaceCenter[3*i]) ||
            isnan(refFaceCenter[3*i]+1) ||
            isnan(refFaceCenter[3*i]+2)) {
          std::cout << "degenerated element (nan):" << i << std::endl;
          refFaceCenter[3*i] = barycentre[0];
          refFaceCenter[3*i+1] = barycentre[1];
          refFaceCenter[3*i+2] = barycentre[2];
        }
#endif


#ifdef INFINITY
        if (isinf(ABS(refFaceCenter[3*i])) ||
            isinf(ABS(refFaceCenter[3*i]+1)) ||
            isinf(ABS(refFaceCenter[3*i]+2))) {
          std::cout << "degenerated element (inf) :" << i << std::endl;
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
    std::vector<double> & refCellCenterCoords = *_cellCenterCoords;
    std::vector<double> & refCellVolume = *_cellVolume;

    int nStandardElement = _nElts - _nPolyhedra;
    if (nStandardElement > 0) {
      std::vector<int>  faceConnectivityIndex(5,0);
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
        refCellCenterCoords[3*i] /= refCellVolume[i];
        refCellCenterCoords[3*i+1] /= refCellVolume[i];
        refCellCenterCoords[3*i+2] /= refCellVolume[i];
      }
    }

    if (_nPolyhedra > 0) {

      //
      // Polyedra splitting

      // Not yet implemented
      bft::bft_error(__FILE__, __LINE__, 0, "Not implemented yet\n");

//       fvm::fvm_tesselation_t *fvm::fvm_tesselation = fvm::fvm_tesselation_create(fvm::FVM_CELL_POLY,
//                                                                   _nPolyhedra,
//                                                                   _polyhedraFaceIndex,
//                                                                   _polyhedraCellToFaceConnectivity,
//                                                                   _polyhedraFaceConnectivityIndex,
//                                                                   _polyhedraFaceConnectivity,
//                                                                   null);

//       fvm::fvm_tesselation_init(fvm::fvm_tesselation, 3,_coords, NULL, NULL);


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


//       fvm::fvm_triangulate_state_t* state = fvm::fvm_triangulate_state_create(maxVertexFace);

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

//           fvm::fvm_triangulate_polygon(3,
//                                   nVertexFace,
//                                   _coords,
//                                   NULL,
//                                   faceConnectivity,
//                                   fvm::FVM_TRIANGULATE_MESH_DEF,
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
//       fvm::fvm_triangulate_state_destroy(state);
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
