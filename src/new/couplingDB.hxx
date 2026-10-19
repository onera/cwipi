#ifndef __COUPLING_DB_H__
#define __COUPLING_DB_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2026  ONERA

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

#include "singleton.hpp"

#include "cwp.h"

#include <map>
#include <string>

/**
 * \cond
 */

namespace cwipi {

  class Coupling;
  class CodeProperties;

  /**
   * \class CouplingDB couplingDB.hxx "couplingDB.hxx"
   * \brief Coupling instances storage.
   *
   * This class stores coupling instances
   *
   */

  class CouplingDB :
    public Singleton <CouplingDB>
  {

    friend class Singleton <CouplingDB>;
    friend class Coupling;

  public:

    /**
     * \brief Constructor.
     *
     */
    CouplingDB();


    /**
     * \brief Destructor.
     *
     */
    virtual ~CouplingDB();


    /**
     * \brief Building and storage a coupling object.
     *
     * This function creates a coupling object and defines its properties.
     *
     * \param [in]  localCodeProperties   Source code
     * \param [in]  cplId                 Coupling identifier
     * \param [in]  coupledCodeProperties Coupled code properties
     * \param [in]  entities_dim          Interface dimension
     * \param [in]  commType              Communication type
     * \param [in]  spatialInterpAlgo     Spatial interpolation algorithm
     * \param [in]  nPart                 Number of interface partition
     * \param [in]  movingStatus          Support moving status
     * \param [in]  recvFreqType          Type of receiving frequency
     *
     */
    void
    couplingCreate
    (
            CodeProperties       &localCodeProperties,
      const std::string          &cplId,
            CodeProperties       &coupledCodeProperties,
      const CWP_Interface_t       entities_dim,
      const CWP_Comm_t            commType,
      const CWP_Spatial_interp_t  spatialInterpAlgo,
      const int                   nPart,
      const CWP_Dynamic_mesh_t    movingStatus,
      const CWP_Time_exch_t       recvFreqType
    );


    /**
     * \brief Delete a coupling object in the database.
     *
     * \param [in]  localCodeProperties Source code
     * \param [in]  cplId               Coupling identifier
     *
     */
    void
    couplingDel
    (
      const CodeProperties &localCodeProperties,
      const std::string    &cplId
    );


    /**
     * \brief Return a coupling object from its identifier.
     *
     * \param [in]  localCodeProperties Source code
     * \param [in]  cplId               Coupling identifier
     *
     * \return Reference to the coupling object
     */
    inline Coupling&
    couplingGet
    (
      const CodeProperties &localCodeProperties,
      const std::string    &cplId
    );


    /**
     * \brief Return the list of coupling objects for a code.
     *
     * \param [in]  localCodeProperties Source code
     *
     * \return Reference to the map of coupling objects
     */
    inline std::map<std::string, Coupling *>&
    couplingListGet
    (
      const CodeProperties &localCodeProperties
    );


    /**
     * \brief Return if a coupling identifier exists.
     *
     * \param [in]  localCodeProperties Source code
     * \param [in]  cplId               Coupling identifier
     *
     * \return status
     */
    inline bool
    couplingIs
    (
      const CodeProperties &localCodeProperties,
      const std::string    &cplId
    );


    /**
     * \brief Begin time step.
     *
     * \param [in]  localCodeProperties Source code
     * \param [in]  current_time        Current time
     *
     */
    void
    timeStepBegin
    (
      const CodeProperties &localCodeProperties,
      const double          current_time
    );


    /**
     * \brief End time step.
     *
     * \param [in]  localCodeProperties Source code
     *
     */
    void
    timeStepEnd
    (
      const CodeProperties &localCodeProperties
    );

  private:

    std::map<const CodeProperties *, std::map<std::string, Coupling *>> &_couplingDB; /*!< Couplings storage for each local code */

  };

}

/**
 * \endcond
 */

#endif // __COUPLING_DB_H__