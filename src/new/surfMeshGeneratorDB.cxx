/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-2017  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more detailstr_options.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/license_width/>.
*/



#include "pdm_part.h"
#include "pdm_timer.h"
#include "pdm.h"
#include "pdm_config.h"
#include "time.h"
#include "pdm_mesh_nodal.h"
#include "pdm_poly_surf_gen.h"
#include "surfMeshGeneratorDB.hxx"
#include <assert.h>
#include <stdio.h>
#include <string.h>

/**
 * \cond
 */

namespace cwipi {

surfMeshGeneratorDB::surfMeshGeneratorDB()
{

}

surfMeshGeneratorDB::~surfMeshGeneratorDB()
{

}

void surfMeshGeneratorDB::createMember(string genName)
{

  surfMeshGenerator* surfMesh = new surfMeshGenerator();
  _dataBase.insert( std::pair<string,surfMeshGenerator*>( genName, surfMesh ));

}

void surfMeshGeneratorDB::destroyMember(string genName)
{

  std::map<string,surfMeshGenerator*>::iterator itr = _dataBase.begin();
  while (itr != _dataBase.end()) {
    if ( itr->first == genName) {
       std::map<string,surfMeshGenerator*>::iterator toErase = itr;
       ++itr;
       _dataBase.erase(toErase);
    } else {
       ++itr;
    }
  }

}



surfMeshGenerator* surfMeshGeneratorDB::memberGet(string genName)
{
  std::map<string,surfMeshGenerator*>::iterator itr = _dataBase.begin();
  while (itr != _dataBase.end()) {
    if ( itr->first == genName) {
       return itr-> second;
       ++itr;
    } else {
       ++itr;
    }
  }
  return nullptr;
}


}//end namespace cwipi


/**
 * \endcond
 */
