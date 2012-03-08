#ifndef __CWIPI_VECTORUTILITIES_H__
#define __CWIPI_VECTORUTILITIES_H__
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

#include <cmath>
#include <algorithm>

#include <bftc_error.h>

namespace cwipi {

  /// Minimum value allowed for geometric computation
  
  const double GEOMETRY_EPS_MIN  = 1e-30; 

  /// Constant value used to compute geomtric epsilon for volume

  const double GEOMETRY_EPS_VOL  = 1e-15;
 
  /// Constant value used to compute geomtric epsilon for surface

  const double GEOMETRY_EPS_SURF = 1e-12;

  /// Minimum distance between two vertices 

  const double GEOMETRY_EPS_DIST = 1e-9; 

  ///
  /// \brief Compute a dynamic geometric epsilon from a characteristic length
  ///
  ///   @param [in]  characteristicLength  Characteristic length
  ///   @param [in]  consEpsilon           Constant part
  ///   @return                            Geometric epsilon
  ///

  inline double geometricEpsilon(const double characteristicLength,
                                 const double constEpsilon)
    
  {
    return std::max(constEpsilon * characteristicLength, GEOMETRY_EPS_MIN);
  }

  ///
  /// \brief Compute cross product
  ///
  ///   @param [in]      vect1  first vector
  ///   @param [in]      vect2  second vector
  ///   @param [out]     cross  vect1 X vect2
  ///

  inline void crossProduct(const double *vect1, 
                           const double *vect2,
                           double *cross)
  { 
    cross[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
    cross[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
    cross[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];
  }


  ///
  /// \brief Compute dot product
  ///
  ///   @param [in]      vect1  first vector
  ///   @param [in]      vect2  second vector
  ///   @return                 vect1 . vect2
  ///

  inline double dotProduct(const double* vect1, 
                           const double* vect2)

  {
    return vect1[0]*vect2[0] + vect1[1]*vect2[1] + vect1[2]*vect2[2];
  }

  ///
  /// \brief Compute the dot product between two vectors
  /// 
  /// @param [in]  p1V1  Coordinates of the first point of the first vector
  /// @param [in]  p2V1  Coordinates of the second point of the first vector
  /// @param [in]  p1V2  Coordinates of the first point of the second vector
  /// @param [in]  p2V2  Coordinates of the second point of the second vector
  /// @return            value of the dot product
  ///

  inline double dotProduct(const double* p1V1, 
                           const double* p2V1,
                           const double* p1V2, 
                           const double* p2V2)
  {

    return (p2V1[0] - p1V1[0])*(p2V2[0] - p1V2[0])
         + (p2V1[1] - p1V1[1])*(p2V2[1] - p1V2[1])
         + (p2V1[2] - p1V1[2])*(p2V2[2] - p1V2[2]);
  }

  ///
  /// \brief Compute Norm
  ///
  ///   @param [in]      vect1  first vector
  ///   @param [in]      vect2  second vector
  ///   @param [inout]   prod_vect vect1 X vect2
  ///
  
  inline double norm(const double* vect)
  {
    return sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
  }

  ///
  /// \brief Compute the norm of a vector
  /// 
  /// @param [in]  p1V1  Coordinates of the first point of the first vector
  /// @param [in]  p2V1  Coordinates of the second point of the first vector
  /// @return             value of the norm
  ///

  inline double norm(const double* p1V1,
                     const double* p2V1)
  {
    return sqrt(dotProduct(p1V1,p2V1,p1V1,p2V1));
  }
  
  ///
  /// \brief Normalized cross product
  ///
  ///   @param [in]      vect1  first vector
  ///   @param [in]      vect2  second vector
  ///   @param [out]     cross  vect1 X vect2
  ///

  inline void normalizedCrossProduct(const double *vect1, 
                                     const double *vect2,
                                     const double tolerance,
                                     double *cross)
  { 
    cross[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
    cross[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
    cross[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];

    double normCP = norm(cross);

    if (normCP < tolerance)
      bftc_error(__FILE__, __LINE__, 0,
                 "Error normalizedCrossProduct : cross product norm is 0 (%12.5e)\n",
                 normCP);

    cross[0] /= normCP;
    cross[1] /= normCP;
    cross[2] /= normCP;

  }

  ///
  /// \brief Compute determinant
  ///
  ///   @param [in]      vect1  First vector
  ///   @param [in]      vect2  Second vector
  ///   @param [in]      vect3  Third vector
  ///   @return                 Determinant   
  ///

  inline double determinant(const double* vect1,
                            const double* vect2,
                            const double* vect3)
  {
    return ((vect1[1] * vect2[2] - vect2[1] * vect1[2]) * vect3[0])
         + ((vect2[0] * vect1[2] - vect1[0] * vect2[2]) * vect3[1])
         + ((vect1[0] * vect2[1] - vect2[0] * vect1[1]) * vect3[2]);
  }

  ///
  /// \brief Compute the cosinus between two vectors
  /// 
  /// @param [in]  vect1      Coordinates of the first vector
  /// @param [in]  vect2      Coordinates of the second vector
  /// @param [in]  tolerance  Tolerance for (||vect|| == 0)                       
  /// @return                 Value of the cosinus
  ///
  
  inline double cosinus(const double* vect1, 
                        const double* vect2,
                        const double tolerance)
  {

    double norm1 = norm(vect1);
    double norm2 = norm(vect2);

    if (norm1 < tolerance || norm2 < tolerance)
      bftc_error(__FILE__, __LINE__, 0,"cosinus : vector norm are 0 (vect1 %f,vect2 %f)\n",
                 norm1, norm2);
    return (dotProduct(vect1,vect2)/(norm1*norm2));
  }


  ///
  /// \brief Compute the cosinus between two vectors
  /// 
  /// @param [in]  p1V1  Coordinates of the first point of the first vector
  /// @param [in]  p2V1  Coordinates of the second point of the first vector
  /// @param [in]  p1V2  Coordinates of the first point of the second vector
  /// @param [in]  p2V2  Coordinates of the second point of the second vector
  /// @return            value of the cosinus
  ///
  
  inline double cosinus (const double* p1V1, 
                         const double* p2V1,
                         const double* p1V2, 
                         const double* p2V2,
                         const double tolerance)
  {

    double norm1 = norm(p1V1, p2V1);
    double norm2 = norm(p1V2, p2V2);

   if (norm1 < tolerance || norm2 < tolerance)
      bftc_error(__FILE__, __LINE__, 0,"cosinus : vector norm are 0 (vect1 %f,vect2 %f)\n",
                     norm1, norm2);
    
    return (dotProduct(p1V1,p2V1,p1V2,p2V2)
            /(norm1*norm2));
  }

  ///
  /// \brief Sign of a dot product
  /// 
  /// @param [in]  vect1   Coordinates of the first vector
  /// @param [in]  vect2   Coordinates of the second vector
  /// @return              Return 1 if the two vectors are in the same direction, else -1
  ///
  
  inline int signDotProduct (const double* vect1,
                             const double* vect2)
  {
    return ((dotProduct(vect1,vect2) >= 0 ) ? 1 : -1 );
  }

}

#endif //__CWIPI_VECTORUTILITIES_H__
