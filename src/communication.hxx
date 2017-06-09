#ifndef __COMMUNICATION_H__
#define __COMMUNICATION_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2013  ONERA

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

#include "codeProperties.hxx"

namespace cwipi {

  /** 
   * \class Communication communication.hxx "communication.hxx"
   * \brief Communication abstract interface
   *
   *  This class is the abstract interface of communication instances.
   *  A communication instance manage MPI communication between coupled codes
   */

  class Communication {
    
  public :

    /**
     *
     * \brief Constructor.
     *
     */

    Communication();

    /**
     * \brief Destructor.
     *
     */

    virtual ~Communication();

    /**
     *
     * \brief Initialize coupling communicators.
     *
     * \param [in]  localCodeProperties   Local code properties
     * \param [in]  cplCodeProperties     Coupled code properties
     *
     */

    void 
    init
    (
     CodeProperties &localCodeProperties, 
     CodeProperties &cplCodeProperties 
     );

    /**
     *
     * \brief Return the communicator type
     *
     */

    virtual CWP_Comm_t 
    commTypeGet() = 0;

    /**
     *
     * \brief Synchronise
     *
     */

    virtual void
    sync
    (
     void *tab, 
     MPI_Datatype mpiType, 
     int tabSize
    ) = 0;

  protected :

    /**
     *
     * \brief Building coupling communicator from \ref CWP_Comm_t.
     *
     */

    virtual void _cplCommCreate(CWP_Comm_t cplCodeCommType) = 0;


  private :

    /**
     *
     * \brief Assigment operator
     *
     */
    
    Communication 
    &operator=
    (const Communication &other);

    /**
     *
     * \brief Copy constructor
     *
     */

    Communication (const Communication& other); 

  protected:

    bool      _isCplRank;               /*!< Is the current rank coupled */
    MPI_Comm  _unionComm;               /*!< Union communicator between coupled codes */
    MPI_Comm  _fvmComm;                 /*!< FVM communicator 
                                          (part of local communicator) */
    MPI_Comm _cplGroup;                  /*!< Coupling group 
                                          (part of merger inter communicator */
    MPI_Comm _cplComm;                  /*!< Coupling communicator 
                                          (part of merger inter communicator */
    int      _locCodeRootRankCplComm;  /*!< Root rank associated to the coupled code
                                          into the coupling communicator */
    int      _cplCodeRootRankCplComm;  /*!< Root rank associated to the coupled code
                                          into the coupling communicator */
    CodeProperties *_localCodeProperties; /*!< Pointer to the local code properties */
    CodeProperties *_cplCodeProperties;   /*!< Pointer to the coupled code properties */
  };
}
#endif //__COMMUNICATION_H__
