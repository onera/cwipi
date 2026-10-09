#ifndef __COMMWITHOUTPART_H__
#define __COMMWITHOUTPART_H__
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
   * \class CommWithoutPart commWithoutPart.hxx "commWithoutPart.hxx"
   * \brief Communication without partitioned support.
   *
   * This class manages communication when the support is not partitioned.
   *
   */

  class CommWithoutPart :
    public Communication {

  public:

    /**
     * \brief Constructor.
     *
     */
    CommWithoutPart();


    /**
     * \brief Destructor.
     *
     */
    virtual ~CommWithoutPart();


    /**
     * \brief Return the communication type.
     *
     * \return Communication type
     */
    inline CWP_Comm_t
    commTypeGet();

  private:

    /**
     * \brief Build the coupling communicator.
     *
     * \param [in]  cplCodeCommType Coupling communicator type
     */
    void
    _cplCommCreate
    (
      CWP_Comm_t cplCodeCommType
    );

    CWP_Comm_t _commType; /*!< Communication type */

  };

  inline CWP_Comm_t
  CommWithoutPart::commTypeGet()
  {
    return _commType;
  }

}

/**
 * \endcond
 */

#endif // __COMMWITHOUTPART_H__