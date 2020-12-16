#ifndef __SURF_MESH_GEN_H__
#define __SURF_MESH_GEN_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2012-2017  ONERA

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
#include <map>
#include <vector>

#include "pdm_printf.h"

#include "cwp.h"


/**
 * \cond
 */

using namespace std;

namespace cwipi {

  class surfMeshGenerator
  {

    public:
      surfMeshGenerator();

      ~surfMeshGenerator();

      void init(int nx, int ny, int nPart, MPI_Comm* comm, double prop, double width, double randomVar);

      void computeMesh();


      double* specialFieldTriGet(int i_part);
      double* specialFieldQuadGet(int i_part);
      double* specialFieldPolyGet(int i_part);

      int nTriGet(int i_part) {
        return _nTri[i_part];
      }

      int nQuadGet(int i_part){
        return _nQuad[i_part];
      }

      int nPolyGet(int i_part){
        return _nPoly[i_part];
      }

      int* connecTriGet (int i_part){
        return _eltsConnecTri[i_part];
      }

      CWP_g_num_t* eltsGnumTriGet(int i_part){
        return _eltsGnumTri[i_part];
      }

      int* connecQuadGet(int i_part){
        return _eltsConnecQuad[i_part];
      }

      CWP_g_num_t* eltsGnumQuadGet(int i_part){
        return _eltsGnumQuad[i_part];
      }

      int* connecPolyGet(int i_part){
        return _eltsConnecPoly[i_part];
      }

      int* connecPolyIndexGet(int i_part){
        return _eltsConnecPolyIndex[i_part];
      }

      CWP_g_num_t* eltsGnumPolyGet(int i_part){
        return _eltsGnumPoly[i_part];
      }

      int nVtxGet(int i_part){
        return _nVtx[i_part];
      }

      double* coordsGet(int i_part){
        return _coords[i_part];
      }

      CWP_g_num_t* vtxGnumGet(int i_part){
        return _vtxGnum[i_part];
      }

      int nEltsGet(int i_part){
        return _nElts[i_part];
      }

      int* connecGet(int i_part){
        return _eltsConnec[i_part];
      }

      int* connecIndexGet(int i_part){
        return _eltsConnecIndex[i_part];
      }

      CWP_g_num_t* eltsGnumGet(int i_part){
        return _eltsGnum[i_part];
      }

      int* faceEdgeIdxGet(int i_part){
        return _faceEdgeIdx[i_part];
      }

      int* faceEdgeGet(int i_part){
        return _faceEdge[i_part];
      }

      int* edgeVtxIdxGet(int i_part){
        return _edgeVtxIdx[i_part];
      }

      int* edgeVtxGet(int i_part){
        return _edgeVtx[i_part];
      }


      int nEdgeGet(int i_part){
        return _nEdge[i_part];
      }

      int nFaceGet(int i_part){
        return _nFace[i_part];
      }


      CWP_g_num_t* faceLNToGNGet(int i_part){
        return (CWP_g_num_t*)_faceLNToGN[i_part];
      }


    private:
      double _motif(double x, double y);
      double _inBox(double x, double y, double x1, double y1 ,double x2, double y2);
      double _inCircle(double x, double y, double R) ;
      double _inBox2(double x, double y) ;

      int      _nx         ; /*!< Number of segment division on the whole domain (global) >*/
      int      _ny         ; /*!< Number of segment division on the whole domain (global) >*/

      PDM_MPI_Comm _comm      ; /*!< MPI Communicator >*/
      int          _rank;
      int          _commSize;
      int          _color;
      double   _prop      ; /*!< Proportion of MPI process owning the interface mesh >*/
      PDM_MPI_Comm _interfComm; /*!< Interface MPI Communicator >*/
      int          _interfCommSize;

      std::vector<double*> _coords;
      std::vector<CWP_g_num_t*> _vtxGnum;

      std::vector<int> _nVtx;
      std::vector<int> _nElts;
      std::vector<int> _nPoly;
      std::vector<int*> _eltsConnecPolyIndex;
      std::vector<int*> _eltsConnecPoly;
      std::vector<CWP_g_num_t*> _eltsGnumPoly;
      std::vector <double*>  _specialFieldPoly;
      std::vector <double*>  _surfVecPoly     ;
      std::vector <double*>  _centerPoly      ;
      std::vector <double*>  _charLengthPoly  ;
      std::vector <int*   >  _isDegPoly       ;

      std::vector<int> _nTri;
      std::vector<int*> _eltsConnecTri;
      std::vector<CWP_g_num_t*> _eltsGnumTri;
      std::vector <double*>  _specialFieldTri;
      std::vector <double*>  _surfVecTri     ;
      std::vector <double*>  _centerTri      ;
      std::vector <double*>  _charLengthTri  ;
      std::vector <int*   >  _isDegTri       ;


      std::vector<int> _nQuad;
      std::vector<int*> _eltsConnecQuad;
      std::vector<CWP_g_num_t*> _eltsGnumQuad;
      std::vector <double*>  _specialFieldQuad;
      std::vector <double*>  _surfVecQuad     ;
      std::vector <double*>  _centerQuad      ;
      std::vector <double*>  _charLengthQuad  ;
      std::vector <int*   >  _isDegQuad       ;

      std::vector<int*> _eltsConnecIndex;
      std::vector<int*> _eltsConnec;
      std::vector<CWP_g_num_t*> _eltsGnum;

      std::vector<int> _nFace   ;
      std::vector<int*> _faceEdgeIdx;
      std::vector<int*> _faceEdge   ;
      std::vector<int> _nEdge   ;
      std::vector<int*> _edgeVtxIdx ;
      std::vector<int*> _edgeVtx    ;
      std::vector<CWP_g_num_t*> _faceLNToGN ;

      double _width;
      double _xmin;
      double _xmax;
      double _ymin;
      double _ymax;

      double _randomVar;
      int    _nPart;
  };//end class surfMeshGenerator
}


/**
 * \endcond
 */


#endif //__SURF_MESH_GEN_H__
