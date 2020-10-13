#ifndef __SURF_MESH_GEN_DB_H__
#define __SURF_MESH_GEN_DB_H__
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
#include <string>


#include "pdm_printf.h"
#include "pdm_mesh_nodal.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_part.h"

#include "cwp.h"
#include "singleton.hpp"
#include "surfMeshGenerator.hxx"


/**
 * \cond
 */


using namespace std;

namespace cwipi {

  class surfMeshGeneratorDB
     : public Singleton <surfMeshGeneratorDB>
  {
    friend class Singleton <surfMeshGeneratorDB>;

    public:
      surfMeshGeneratorDB();

      ~surfMeshGeneratorDB();


      void createMember(string genName);

      void destroyMember(string genName);

      surfMeshGenerator* memberGet(string genName);

    private:

      std::map<string,surfMeshGenerator*> _dataBase;

  };//end class surfMeshGeneratorDB
}


/**
 * \endcond
 */


#endif //__SURF_MESH_GEN_DB_H__
