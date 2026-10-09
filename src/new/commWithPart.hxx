#ifndef __COMMWITHPART_H__
#define __COMMWITHPART_H__
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

#include "communication.hxx"

#include "cwp.h"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class CommWithPart commWithPart.hxx "commWithPart.hxx"
   * \brief Communication with partitioned support.
   *
   * This class manages communication when the support is partitioned
   *
   */

  class CommWithPart :
    public Communication
  {

  public:

    /**
     * \brief Constructor.
     *
     */
    CommWithPart();


    /**
     * \brief Destructor.
     *
     */
    virtual ~CommWithPart() override;


    /**
     * \brief Return the communication type.
     *
     * \return Communication type
     */
    inline CWP_Comm_t
    commTypeGet() override;

  private:

    /**
     * \brief Build the coupling communicator.
     *
     * \param [in]  cplCodeCommType Coupling communicator type
     */
    void
    cplCommCreate
    (
      CWP_Comm_t cplCodeCommType
    ) override;

    CWP_Comm_t _commType; /*!< Communication type */

  };

  inline CWP_Comm_t
  CommWithPart::commTypeGet()
  {
    return _commType;
  }

}

/**
 * \endcond
 */

#endif // __COMMWITHPART_H__