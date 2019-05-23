#ifndef __GEOMLOCATION_H__
#define __GEOMLOCATION_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2012  ONERA

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

#include "mesh.hxx"
#include "geometry.hxx"
#include "field.hpp"

namespace cwipi {

  class GeomLocation:
    public Geometry {
  
    public:
      GeomLocation();
      
      virtual ~GeomLocation();
     
     void locate_cell_point_setting_request(int id_dist);
     
     void locate_cell_point_setting_surface(int id_dist);     
      void locate_cell_point_setting_code_1_local() ;
      
      void locate_cell_point_compute(int id_dist)             ;
      void locate_cell_point_get(int id_dist)      ;
      void locate_cell_point_get_cpl(int id_dist)      ;
      
      void redistribution_cell_point_request(int id_gnum_location) ;
      void redistribution_cell_point_set    (int id_gnum_location) ;
      void location_compute                 (int id_gnum_location) ;
                  
      void location_get(int id_gnum_location)      ;
      void location_get_cpl(int id_gnum_location)      ;
      void redistribution_cell_point_filling_of_redistribution_array ();
      void redistribution_index_communication(int tag)    ;
      void redistribution_communication()          ;
      void redistribution_communication2()          ;
      void redistribution_wait_and_targets_array_filling();

                 
      virtual void issend(Field <double>* sendingField);

      virtual double* interpolate (Field <double>* referenceField);  

    private:

      
   


  }; //end GeomLocation
  

  
}
#endif //__GEOMLOCATION_H__


