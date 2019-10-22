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
#include "singleton.hpp"

using namespace std;

namespace cwipi {

  class surfMeshGenerator
     : public Singleton <surfMeshGenerator>
  {
    friend class Singleton <surfMeshGenerator>;
  
    public:
      surfMeshGenerator();

      ~surfMeshGenerator();

      void init(int nx, int ny, MPI_Comm* comm, double prop, double width, double randomVar);
  
      void computeMesh();
      
      int nTriGet() {
        return _nTri;
      }
      
      int nQuadGet(){
        return _nQuad;
      }
      
      int nPolyGet(){
        return _nPoly;
      }

      int* connecTriGet (){
        return _eltsConnecTri;
      }
      
      CWP_g_num_t* eltsGnumTriGet(){
        return _eltsGnumTri;
      }      
      
      int* connecQuadGet(){
        return _eltsConnecQuad;
      }

      CWP_g_num_t* eltsGnumQuadGet(){
        return _eltsGnumQuad;
      }
      
      int* connecPolyGet(){
        return _eltsConnecPoly;
      }
      
      int* connecPolyIndexGet(){
        return _eltsConnecPolyIndex;
      }

      CWP_g_num_t* eltsGnumPolyGet(){
        return _eltsGnumPoly;
      }
      
      int nVtxGet(){
        return _nVtx;
      }

      double* coordsGet(){
        return _coords;
      }

      CWP_g_num_t* vtxGnumGet(){
        return _vtxGnum;
      }
      
      int nEltsGet(){
        return _nElts;
      }    

      int* connecGet(){
        return _eltsConnec;
      }
      
      int* connecIndexGet(){
        return _eltsConnecIndex;
      }

      CWP_g_num_t* eltsGnumGet(){
        return _eltsGnum;
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
      
      double* _coords;
      CWP_g_num_t* _vtxGnum;
      
      int _nVtx;
      int _nElts;
      int _nPoly;
      int* _eltsConnecPolyIndex;
      int* _eltsConnecPoly;
      CWP_g_num_t* _eltsGnumPoly;
      
      int _nTri;
      int* _eltsConnecTri;
      CWP_g_num_t* _eltsGnumTri;
      
      int _nQuad;
      int* _eltsConnecQuad;
      CWP_g_num_t* _eltsGnumQuad;

      int* _eltsConnecIndex;
      int* _eltsConnec;      
      CWP_g_num_t* _eltsGnum;
                  
      
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
