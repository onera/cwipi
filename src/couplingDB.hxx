#ifndef __COUPLING_DB_H__
#define __COUPLING_DB_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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

#include <map>
#include <string>

#include "singleton.hpp"
#include "cwp.h"

using namespace std;

namespace cwipi {

  class Coupling;
  class CodeProperties;
  
  /** 
   * \class Coupling couplingDB.hxx "couplingDB.hxx"
   * \brief Coupling instances storage
   *
   *  This class stores coupling instances
   * 
   */

  class CouplingDB : public Singleton <CouplingDB>
  {
    friend class Singleton <CouplingDB>;
    friend class Coupling;
    
  public:

    /**
     * \brief Building and storage a coupling object 
     *
     * This function creates a coupling object and defines its properties.
     *
     * \param [in]  localCodeProperties  Source code
     * \param [in]  cplId                Coupling identifier
     * \param [in]  cplCodeProperties    Coupled code properties
     * \param [in]  commType             Communication type
     * \param [in]  geomAlgo             Geometric algorithm
     * \param [in]  supportType          Support type
     * \param [in]  nPart                Number of interface partition 
     * \param [in]  movingStatus         Support moving status
     * \param [in]  recvFreqType         Type of receiving frequency
     *
     */

    void 
    couplingCreate
    (
     const CodeProperties        &localCodeProperties,
     const string                &cplId,
     const CodeProperties        &coupledCodeProperties,
     const CWP_Comm_t           commType,
     const CWP_Geom_t           geomAlgo,
     const CWP_Support_t        supportType,
     const int                    nPart,
     const CWP_Displacement_t  movingStatus,
     const CWP_Freq_t           recvFreqType
    );
    
    /**
     * \brief Deletion a coupling object int the database.
     *
     * \param [in]  localCodeProperties  Source code
     * \param [in]  cplId                Coupling identifier
     *
     */

    void 
    couplingDel
    (
     const CodeProperties &localCodeProperties,
     const string         &cplId
    );
    
    /**
     * \brief Return a coupling object from it identifier
     *
     * \param [in]  localCodeProperties  Source code
     * \param [in]  cplId                Coupling identifier
     *
     */

    inline Coupling& 
    couplingGet
    (
     const CodeProperties &localCodeProperties,
     const string &cplId
    );
    
    /**
     * \brief Return if a coupling identifier exists  
     *
     * \param [in]  localCodeProperties  Source code
     * \param [in]  cplId                Coupling identifier
     *
     * \return status
     */

    inline bool 
    couplingIs
    (
     const CodeProperties &localCodeProperties,
     const string &cplId
    );
    
  private:

    /**
     * \brief Constructor
     *
     */

    CouplingDB();

    /**
     * \brief Copy constructor
     *
     */

    CouplingDB
    (
     const CouplingDB &other
    );

    /**
     * \brief Default assignment (unavailable)
     *
     */

    CouplingDB & 
    operator=
    (
     const CouplingDB &other
    );

    /**
     * \brief Destructor
     *
     */

    virtual ~CouplingDB();
    
  private:
    map < const CodeProperties *, map < string, Coupling * > > & _couplingDB; /*!< Couplings storage for each local Code */

  };
}

#endif
