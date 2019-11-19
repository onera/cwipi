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

using namespace std;

namespace cwipi {

  class surfMeshGenerator
  {
  
    public:
      surfMeshGenerator();

      ~surfMeshGenerator();

      void init(int nx, int ny, int nPart, MPI_Comm* comm, double prop, double width, double randomVar);
  
      void computeMesh();
      
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

      
    private:
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
      
      std::vector<int> _nTri;
      std::vector<int*> _eltsConnecTri;
      std::vector<CWP_g_num_t*> _eltsGnumTri;
      
      std::vector<int> _nQuad;
      std::vector<int*> _eltsConnecQuad;
      std::vector<CWP_g_num_t*> _eltsGnumQuad;

      std::vector<int*> _eltsConnecIndex;
      std::vector<int*> _eltsConnec;      
      std::vector<CWP_g_num_t*> _eltsGnum;
                  
      
      double _width;
      double _xmin;
      double _xmax;
      double _ymin;
      double _ymax;
   
      double _randomVar;
      int    _nPart;
  };//end class surfMeshGenerator      
}



#endif //__SURF_MESH_GEN_H__
