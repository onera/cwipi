#ifndef __SUPPORT_H__
#define __SUPPORT_H__
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

namespace cwipi {

  /** 
   * \class Support support.hxx "support.hxx"
   * \brief Geometry support
   *
   *  This class computes defines th geometry support (mesh)
   * 
   */

  class Support {
    
  public:

    /**
     * \brief Constructor
     *
     */

    Support();

    /**
     * \brief Destructor
     *
     */

    virtual ~Support();


  private:
    
  //   Support &operator=(const Support &other);  /*!< Assigment operator not available */
  //   Support (const Support& other);            /*!< Copy constructor not available */

  // private:
  //   std::vector<double>         *_tmpVertexField;  // Evite une allocation a chaque extrapolation
  //   std::vector<double>         *_tmpDistantField; //TODO: Fusionner _tmpDistantField utiliser pour exchange
  //                                          // et les comm asynchrones
  //   std::map<int, std::vector<double> * > &_tmpDistantFieldsIssend; //TODO: temporaire : A revoir lors
  //                                                                   // de la restructuration
  //   std::map<int, const double * >        &_tmpLocalFieldsIrecv;
  //   std::map<int, const char * >          &_tmpExchangeNameIrecv;
  //   std::map<int, int >                   &_tmpStrideIrecv;
  //   std::map<int, int >                   &_tmpTimeStepIrecv;
  //   std::map<int, double >                &_tmpTimeValueIrecv;
  //   std::map<int, const char * >          &_tmpFieldNameIrecv;
  //   Support                                *support;
  };

}

#endif //__SUPPORT_H__
