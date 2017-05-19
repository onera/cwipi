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
     * \param [in]  globalComm      MPI communicator containing all processes 
     *                              of all codes
     * \param [in]  n_codes         Number of codes on the current rank
     * \param [in]  code_names      Codes names on the current rank
     * \param [in]  is_coupled_rank Current rank is it a coupled rank
     * \param [out] 
     *
     * \return                   Current code intra-communicator
     *
     */

    void 
    init
    (
     const MPI_Comm     globalComm,
     const int          n_codes,
     const char**       code_names, 
     const CWP_Status_t is_coupled_rank,
     MPI_Comm           *intra_comms       
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
      * \parm[in]   localCodeName  Local code name
      * 
      * \return  MPI Intra communicator
      *
      */

    inline const MPI_Comm &
    intraCommGet
    (
     const string & localCodeName
    ) const;

     
    /**
     * \brief Return MPI communicator containing all processes of all codes.
     *
     * \return  Global MPI Communicator
     *
     */

    inline const MPI_Comm &
    globalCommGet() const;

    
    /**
     * \brief Return the code properties.
     *
     * \param [in]  codeName  Code name
     *
     * \return      Properties
     *
     */

    inline const CodeProperties &
    codePropertiesGet
    (
     const string &codeName
    ) const;


    /**
     * \brief Set a control paramater.
     *
     * \param [in]  localCodeName   Local code name
     * \param [in]  name            Parameter name
     * \param [in]  value           Initial value 
     *
     */
    
    template < typename T > 
    void 
    ctrlParamAdd
    (
     const string &localCodeName, 
     const string &name, 
     const T       value
    );


    /**
     * \brief Set a control paramater.
     *
     * \param [in]  localCodeName   Local code name
     * \param [in]  name            Parameter name
     * \param [in]  value           Initial value 
     *
     */
    
    template < typename T > 
    void 
    ctrlParamSet
    (
     const string &localCodeName, 
     const string &name, 
     const T       value
    );

    
    /**
     * \brief Cancel a control paramater.
     *
     * \param [in]  localCodeName   Local code name
     * \param [in]  name            Parameter name
     *
     */

    template < typename T > 
    void 
    ctrlParamCancel
    (
     const string &localCodeName, 
     const string &name
    );

    /**
     * \brief Return the number of parameters
     *
     * \param [in]  codeName  Code name
     *
     * \return  Number of parameters
     *
     */

    template < typename T > 
    int 
    ctrlParamNGet
    (
     const string &codeName
    ) const;


    /**
     * \brief Return of the list of parameters
     *
     * \param [in]  codeName  Code name
     *
     * \return  List of parameters
     *
     */

    template < typename T > 
    char **
    ctrlParamListGet
    (
     const string &codeName
    ) const;


    /**
     * \brief Chek name parameter existence
     *
     * \param [in]  codeName   Code name
     * \param [in]  name       Parameter to check
     *
     * \return  1 : true / 0 : false
     *
     */

    template < typename T > 
    int
    ctrlParamIs
    (
     const string &codeName,
     const string &name
    ) const;


    /**
     * \brief Get the value of a control paramater.
     *
     * \param [in]  codeName  Code name
     * \param [in]  name      Parameter name
     *
     * \return             Value           
     *
     */

    template < typename T > 
    const T &
    ctrlParamGet
    (
      const string &codeName,
      const string &name
    );

    /**
     * \brief Reduce a parameter through a list of codes. The available processes
     *        are sum, max and min. 
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
     const CWP_Op_t  op, 
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
     * \param [in]  codeName  Code name to lock
     *
     */

    inline void 
    lock
    (
    const string &codeName
    );

    /**
     * \brief unlock access to local parameters from a distant code  
     *
     * \param [in]  codeName  Code name to unlock
     *
     */

    inline void 
    unLock
    (
    const string &codeName
    );

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
     * \param [in]  codeName  Local code name
     *
     */

    template < typename T > 
    void 
    _issendParameters
    (
     const string &codeName
    );

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
     * \param [in]  codeName  Local code name
     *
     */

    void 
    _issendLock
    (
    const string &codeName
    );

    /**
     * \brief Non blocking sending of parameters issend cancellation
     *
     * \param [in]  codeName  Local code name
     *
     */

    template < typename T > 
    void 
    _issendParameterCancel
    (
     const string &codeName
    );

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
    MPI_Comm                          _globalComm;             /*!< Global communicator */  
    map <string, CodeProperties * > & _codePropertiesDB;       /*!< Distant code 
                                                                           properties data base */
    map <string, CodeProperties * > & _locCodePropertiesDB;       /*!< Local code properties */

    map < string, map < string, map < string, vector<MPI_Request> * > > > _issendMPIrequest; /*!< MPI Request for 
                                                                                  parameter sending */
    map < string, map<string, string> >               _issendNameBuffs;               /*!< Issend buffers to storage 
                                                                           parameter names */
    map < string, map<string, vector<unsigned char > > >_issendValBuffs;             /*!< Issend buffers to storage 
                                                                           parameter values */
    string                            _recvNameBuff;                  /*!< Receive buffer to storage 
                                                                           distant parameter names */
    vector<unsigned char>             _recvValBuff;                   /*!< Receive buffer to storage 
                                                                           distant parameter values */
    map <string, int>                 _lockStatus;                    /*!< Parameters lock status of
                                                                           distant applications */ 
    map <string, map <string, MPI_Request > > _issendLockMPIrequest;          /*!< MPI request for parameter 
                                                                           lock status sending */ 

    int                               _tagLockStatusBase;             /*!< The MPI tag base for lock status exchanges */          

    int                               _tagIntParameterBase;              /*!< The MPI tag base for int parameter exchanges */
    int                               _tagDoubleParameterBase;              /*!< The MPI tag base for double parameter exchanges */
    int                               _tagStringParameterBase;              /*!< The MPI tag base for string parameter exchanges */
    
    bool                               _isLocalCodeRootrank;          /*!< Current is it a local root rank 
                                                                      *   in the global communicator */ 
  private:
    static const int _nIssend;                                        /*!< Number of issend 
                                                                           to send parameters */
  };
}

#endif /* __CODE_PROPERTIES_DB_H__ */
