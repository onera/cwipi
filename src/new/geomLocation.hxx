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
     
     void locate_setting_request(int id_dist);
     
     void locate_setting_surface(int id_dist);     
      void locate_setting_code_1_local() ;
      
      void locate_compute(int id_dist)             ;
      void locate_get(int id_dist)      ;
      void locate_get_cpl(int id_dist)      ;
      
      void broadcasting_request(int id_gnum_location) ;
      void broadcasting_set    (int id_gnum_location) ;
      void location_compute                 (int id_gnum_location) ;
                  
      void location_get(int id_gnum_location)      ;
      void location_get_cpl(int id_gnum_location)  ;
      void broadcasting_filling_of_broadcasting_array ();
      void broadcasting_index_communication()    ;
      void broadcasting_communication()          ;
      void broadcasting_communication2()          ;
      void broadcasting_wait_and_targets_array_filling();

                 
      virtual void issend(Field <double>* sendingField);

      virtual double* interpolate (Field <double>* referenceField);  

    private:

      
   


  }; //end GeomLocation
  

  
}
#endif //__GEOMLOCATION_H__


