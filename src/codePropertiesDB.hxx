#ifndef __CODE_PROPERTIES_DB_H__
#define __CODE_PROPERTIES_DB_H__
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

#include <mpi.h>
#include <map>
#include <vector>

#include <bftc_printf.h>

#include "cwp.h"
#include "singleton.hpp"

using namespace std;

namespace cwipi {

  class CodeProperties;

  /** 
   * \class CodePropertiesDB 
   *        codePropertiesDB.hxx 
   *        "codePropertiesDB.hxx"
   *
   * \brief Codes properties management.
   *
   *  This class manages codes properties :
   *  - Local control parameters,
   *  - Distant control parameters,
   *  - MPI communicators
   *  - . 
   * 
   */

  class CodePropertiesDB 
    : public Singleton <CodePropertiesDB>
  {

    friend class Singleton <CodePropertiesDB>;

  public:

    /**
     * \brief MPI Communicator Initialization.
     *
     * This function builds the current code intra-communicator from
     * the current name and the MPI communicator containing all processes of
     * all codes.
     *
     * \param [in]  name         Current code name
     * \param [in]  globalComm   MPI communicator containing all processes 
     *                           of all codes
     *
     * \return                   Current code intra-communicator
     *
     */

    MPI_Comm 
    init
    (
     const char* name, 
     const MPI_Comm globalComm
    );

    /**
     * \brief Taking into account of a proxy function to redirect output
     *
     * \param [in]  proxyFunction  Function
     *
     */

    inline void 
    printfProxySet
    (
     bftc_printf_proxy_t *const proxyFunction
    );

    /**
     * \brief Return local code MPI intra communicator.
     *
     * \return  MPI Intra communicator
     *
     */

    inline const MPI_Comm &
    intraCommGet() const;

    /**
     * \brief Return MPI communicator containing all processes of all codes.
     *
     * \return  Global MPI Communicator
     *
     */

    inline const MPI_Comm &
    globalCommGet() const;

    /**
     * \brief Return the current code first rank into the global communicator.
     *
     * \return MPI rank
     *
     */

    inline const int &
    locFirstRankGet() const;

    /**
     * \brief Return the current code last rank into the global communicator.
     *
     * \return MPI rank
     *
     */

    inline const int &
    locLastRankGet() const;

    /**
     * \brief Return the distant code first rank into the global communicator.
     *
     * \param [in]  codeName  Code name
     *
     * \return      MPI rank
     *
     */

    inline const int &
    distFirstRankGet
    (
     const string &codeName
    ) const;

    /**
     * \brief Return the distant code last rank into the global communicator.
     *
     * \param [in]  codeName  Code name
     *
     * \return      MPI rank
     *
     */

    inline const int &
    distLastRankGet
    (
     const string &codeName
    ) const;

    /**
     * \brief Return the distant code properties.
     *
     * \param [in]  codeName  Code name
     *
     * \return      Properties
     *
     */

    inline const CodeProperties &
    distCodePropertiesGet
    (
     const string &codeName
    ) const;

    /**
     * \brief Return the current code properties.
     *
     * \param [in]  codeName  Code name
     *
     * \return      Properties
     *
     */

    inline const CodeProperties &
    locCodePropertiesGet() const ;

    /**
     * \brief Return number of local parameters
     *
     * \return  Number of parameters
     *
     */

    template < typename T > 
    int 
    locCtrlParamNGet() const;

    /**
     * \brief Return of local parameters
     *
     * \return  List of parameters
     *
     */

    template < typename T > 
    char **
    locCtrlParamListGet() const;

    /**
     * \brief Chek name parameter
     *
     * \param [in]  name  Parameter name to check
     *
     * \return  1 : true / 0 : false
     *
     */

    template < typename T > 
    int
    locCtrlParamIs
    (
     const string &name
    ) const;

    /**
     * \brief Add a control paramater.
     *
     * \param [in]  name   Parameter name
     * \param [in]  value  Initial value 
     *
     */

    template < typename T > 
    void 
    locCtrlParamAdd
    (
     const string &name, 
     const T       value
    );

    /**
     * \brief Set a control paramater.
     *
     * \param [in]  name   Parameter name
     * \param [in]  value  Initial value 
     *
     */

    template < typename T > 
    void 
    locCtrlParamSet
    (
     const string &name, 
     const T       value
    );

    /**
     * \brief Get a control paramater.
     *
     * \param [in]  name   Parameter name
     *
     * \return             Value           
     *
     */

    template < typename T > 
    const T &
    locCtrlParamGet
    (
     const string &name
    );

    /**
     * \brief Cancel a control paramater.
     *
     * \param [in]  name   Parameter name
     *
     */

    template < typename T > 
    void 
    locCtrlParamCancel
    (
     const string &name
    );

    /**
     * \brief Return number of local parameters
     *
     * \param [in]  codeName  Code name
     *
     * \return  Number of parameters
     *
     */

    template < typename T > 
    int 
    distCtrlParamNGet
    (
     const string &codeName
    ) const;


    /**
     * \brief Return of local parameters
     *
     * \param [in]  codeName  Code name
     *
     * \return  List of parameters
     *
     */

    template < typename T > 
    char **
    distCtrlParamListGet
    (
     const string &codeName
    ) const;


    /**
     * \brief Chek name parameter
     *
     * \param [in]  codeName  Code name
     * \param [in]  name  Parameter name to check
     *
     * \return  1 : true / 0 : false
     *
     */

    template < typename T > 
    int
    distCtrlParamIs
    (
     const string &codeName,
     const string &name
    ) const;


    /**
     * \brief Get a distant control paramater.
     *
     * \param [in]  codeName  Code name
     * \param [in]  name      Parameter name
     *
     * \return             Value           
     *
     */

    template < typename T > 
    const T &
    distCtrlParamGet
    (
      const string &codeName,
      const string &name
    );

    /**
     * \brief Reduce a parameter.
     *
     * \param [in]  op     Operator from \ref CWP_Op_t
     * \param [in]  name   Parameter name
     * \param [in]  nCode  Number of code
     * \param       ...    Code names
     *
     * \return             Operation result
     *
     */

    template < typename T > 
    void
    ctrlParamReduce
    (
     const CWP_Op_t op, 
     const string    &name,
     T               *res,
     const int        nCode,
     va_list         *pa                       
    );

    /**
     * \brief Dump properties.  
     *
     */

    void 
    dump();

    /**
     * \brief Lock access to local parameters from a distant code  
     *
     */

    inline void 
    lock();

    /**
     * \brief unlock access to local parameters from a distant code  
     *
     */

    inline void 
    unLock();

  private:

    /**
     * \brief Default constructor.
     *
     */

    CodePropertiesDB();

    /**
     * \brief Copy constructor (unavailable)
     *
     */

    CodePropertiesDB
    (
     const CodePropertiesDB &other
     );

    /**
     * \brief Default assignment (unavailable)
     *
     */

    CodePropertiesDB & 
    operator=
    (
     const CodePropertiesDB &other
     );

    /**
     * \brief Destructor.
     *
     */

    virtual ~CodePropertiesDB();

    /**
     * \brief Parameters non blocking send
     *
     */

    template < typename T > 
    void 
    _issendParameters();

    /**
     * \brief Compute the buffer length to send int values
     *
     * \param[in] locCtrlParam  Local control parameters  
     *
     * \return   Buffer length
     *
     */

    inline int 
    _issendLBuffGet
    (
     map <string, int> &locCtrlParam
     );

    /**
     * \brief Compute the buffer length to send double values
     *
     * \param[in] locCtrlParam  Local control parameters  
     *
     * \return   Buffer length
     *
     */

    inline int
    _issendLBuffGet
    (
     map <string, double> &locCtrlParam
    );

    /**
     * \brief Compute the buffer length to send string values
     *
     * \param[in] locCtrlParam  Local control parameters  
     *
     * \return   Buffer length
     *
     */

    inline int
    _issendLBuffGet
    (
     map <string, string> &locCtrlParam
    );

    /**
     * \brief Copy int values into the buffer
     *
     * \param[in] buff            Buffer  
     * \param[in] locCtrlParam  Local control parameters  
     *
     */

    inline void
    _issendBuffCopy
    (
     unsigned char      *buff,
     map <string, int> &locCtrlParam
     );

    /**
     * \brief Copy double values into the buffer
     *
     * \param[in] buff            Buffer  
     * \param[in] locCtrlParam  Local control parameters  
     *
     */

    inline void
    _issendBuffCopy
    (
     unsigned char        *buff,
     map <string, double> &locCtrlParam
     );
    
    /**
     * \brief Copy double values into the buffer
     *
     * \param[in] buff            Buffer  
     * \param[in] locCtrlParam  Local control parameters  
     *
     */

    inline void
    _issendBuffCopy
    (
     unsigned char        *buff,
     map <string, string> &locCtrlParam
    );


    /**
     * \brief Define int values from the received buffer
     *
     * \param[in] lNames          Total length of parameters names 
     * \param[in] buff            Buffer  
     * \param[in] distCtrlParam   Distant control parameters  
     *
     */

    inline void
    _irecvBuffCopy
    (
     const int            lNames,
     const unsigned char *buff,
     map <string, int>   &distCtrlParam
     );


    /**
     * \brief Define double values from the received buffer
     *
     * \param[in] lNames          Total length of parameters names 
     * \param[in] buff            Buffer  
     * \param[in] distCtrlParam   Distant control parameters  
     *
     */

    inline void
    _irecvBuffCopy
    (
     const int             lNames,
     const unsigned char  *buff,
     map <string, double> &distCtrlParam
     );

    /**
     * \brief Define string values from the received buffer
     *
     * \param[in] lNames          Total length of parameters names 
     * \param[in] buff            Buffer  
     * \param[in] distCtrlParam   Distant control parameters  
     *
     */

    inline void
    _irecvBuffCopy
    (
     const int             lNames,
     const unsigned char  *buff,
     map <string, string> &distCtrlParam
     );

    /**
     * \brief Lock status non blocking sending
     *
     */

    void 
    _issendLock();

    /**
     * \brief Non blocking sending of parameters issend cancellation
     *
     */

    template < typename T > 
    void 
    _issendParameterCancel();

    /**
     * \brief Get parameter lock status of a distant code
     *
     * \param [in]  codeName   Distant code name
     *
     */

    void 
    _lockStatusGet
    (
     const string &codeName
    );

    /**
     * \brief Receive parameter from a distant code
     *
     * \param [in]  codeName   Distant code name
     *
     */

    template < typename T > 
    void 
    _irecvParameters
    (
     const string &codeName
    );

  private:
    map <string, CodeProperties * > & _distCodePropertiesDB;       /*!< Distant code 
                                                                           properties data base */
    CodeProperties                  * _locCodeProperties;           /*!< Local code properties */

    map < string, map < string, vector<MPI_Request> * > > _issendMPIrequest; /*!< MPI Request for 
                                                                                  parameter sending */
    map<string, string>               _issendNameBuffs;               /*!< Issend buffers to storage 
                                                                           parameter names */
    map<string, vector<unsigned char > > _issendValBuffs;             /*!< Issend buffers to storage 
                                                                           parameter values */
    string                            _recvNameBuff;                  /*!< Receive buffer to storage 
                                                                           distant parameter names */
    vector<unsigned char>             _recvValBuff;                   /*!< Receive buffer to storage 
                                                                           distant parameter values */
    map <string, int>                 _distLockStatus;             /*!< Parameters lock status of
                                                                           distant applications */ 
    map <string, MPI_Request >        _issendLockMPIrequest;          /*!< MPI request for parameter 
                                                                           lock status sending */ 
    int                               _issendLockStatus;              /*!< Parameter lock status */
  private:
    static const int _nIssend;                                        /*!< Number of issend 
                                                                           to send parameters */
  };
}

#endif /* __CODE_PROPERTIES_DB_H__ */
