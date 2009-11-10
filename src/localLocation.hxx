/*
 * localLocation.hxx
 *
 *  Created on: Oct 16, 2009
 *      Author: equemera
 */

#ifndef LOCATIONTOLOCALMESH_HXX_
#define LOCATIONTOLOCALMESH_HXX_

#include <string>
#include <vector>

namespace cwipi
{

class LocationToLocalMesh
{
public:
  LocationToLocalMesh();
  virtual ~LocationToLocalMesh();

private :
  int    *_barycentricCoordinatesIndex;
  double *_barycentricCoordinates;
  int     _nDistantpoint;
  int    *_location;
};

}

#endif /* LOCALLOCATION_HXX_ */
