#ifndef __CODE_PROPERTIES_DB_I_H__
#define __CODE_PROPERTIES_DB_I_H__
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

#include <cassert>
#include <vector>
#include <typeinfo>

#include "codeProperties.hxx"
#include "codePropertiesDB.hxx"

using namespace std;

namespace cwipi {

  /**
   * \brief Taking into account of a proxy function to redirect output
   *
   * \param [in]  proxyFunction  Function
   *
   */

  void 
  CodePropertiesDB::printfProxySet
  (
   bftc_printf_proxy_t *const proxyFunction
  )
  {
    if (proxyFunction != NULL)
      bftc_printf_proxy_set(proxyFunction);
  }

  /**
   * \brief Return local code MPI intra communicator.
   *
   * \parm[in]   localCodeName  Local code name
   * 
   * \return  MPI Intra communicator
   *
   */

  const MPI_Comm &
  CodePropertiesDB::intraCommGet(const string & localCodeName) const
  {
    return _locCodeProperties[localCodeName]->intraCommGet();
  }

  /**
   * \brief Return MPI communicator containing all processes of all codes.
   *
   * \return  Global MPI Communicator
   *
   */

  const MPI_Comm &
  CodePropertiesDB::globalCommGet() const
  {    
    return _globalComm;
  }

    
  /**
   * \brief Return the code properties.
   *
   * \param [in]  codeName  Code name
   *
   * \return      Properties
   *
   */

  inline const CodeProperties &
  CodePropertiesDB::codePropertiesGet
  (
   const string &codeName
  ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);
    if (p == _codePropertiesDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", codeName.c_str());
    assert(p->second != NULL);
    return *(p->second);
  }


  /**
   * \brief Add a control paramater.
   *
   * \param [in]  localCodeName   Local code name
   * \param [in]  name            Parameter name
   * \param [in]  value           Initial value 
   *
   */

  template < typename T > 
  void 
  CodePropertiesDB::ctrlParamAdd
  (
   const string &localCodeName, 
   const string &name, 
   const T       value
  )
  {
    const map <string, CodeProperties * >::iterator p = 
      _locCodeProperties.find(localCodeName);
    if (p == _locCodeProperties.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", localCodeName.c_str());
    p->second->ctrlParamAdd(name, value);
  }


  /**
   * \brief set a control paramater.
   *
   * \param [in]  localCodeName   Local code name
   * \param [in]  name            Parameter name
   * \param [in]  value           Initial value 
   *
   */

  template < typename T > 
  void 
  CodePropertiesDB::ctrlParamSet
  (
   const string &localCodeName, 
   const string &name, 
   const T       value
  )
  {
    const map <string, CodeProperties * >::iterator p = 
      _locCodeProperties.find(localCodeName);
    if (p == _locCodeProperties.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", localCodeName.c_str());
    p->second->ctrlParamSet(name, value);
  }

    
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
   * \brief Cancel a control paramater.
   *
   * \param [in]  name   Parameter name
   *
   */

  template < typename T > 
  void 
  CodePropertiesDB::locCtrlParamCancel
  (
   const string &name
  )
  {
    _locCodeProperties->ctrlParamCancel<T>(name);
  }

  /**
   * \brief Get a distant control paramater.
   *
   * \param [in]  name   Parameter name
   *
   * \return             Value           
   *
   */
  
  template < typename T > 
  const T &
  CodePropertiesDB::distCtrlParamGet
  (
   const string &codeName,
   const string &name
  ) 
  {
    _irecvParameters<T>(codeName);

    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);

    if (p == _codePropertiesDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' code not found \n", codeName.c_str());
    T *value;
    p->second->ctrlParamGet(name, &value);
    return *value;
  }

  /**
   * \brief Reduce a parameter.
   *
   * \param [in]  op     Operator from \ref CWP_Op_t
   * \param [in]  name   Parameter name
   * \param [out] res    Result          
   * \param [in]  nCode  Number of code
   * \param       ...    Code names
   *
   * \return             Operation result
   *
   */

  template < typename T > 
  void
  CodePropertiesDB::ctrlParamReduce
  (
   const CWP_Op_t  op, 
   const string     &name,
   T                *res,
   const int         nCode,
   va_list          *pa
  )
  {

    T *valLoc;
    _locCodeProperties->ctrlParamGet(name, &valLoc);
    *res =  *valLoc;

    for (int k = 0; k < nCode; k++) {

      string codeName = string((char*)va_arg(*pa, char *));

      cout << "reduce " << typeid(T).name() << " " << codeName.c_str() << endl;

      const map <string, CodeProperties * >::iterator p = 
        _codePropertiesDB.find(codeName);

      if (p == _codePropertiesDB.end())
        bftc_error(__FILE__, __LINE__, 0,
                   "'%s' code not found \n", codeName.c_str());

      T *distParam;
      p->second->ctrlParamGet(name, &distParam);

      switch (op) {
      case CWP_OP_MIN:
        *res = max(*distParam, *res);
        
        break;
      
      case CWP_OP_MAX:
        *res = min(*distParam, *res);
        break;

      case CWP_OP_SUM:
        *res += *distParam;
        break;
      }
    }
  }

  /**
   * \brief Lock access to local parameters from a distant code  
   *
   */
  
  void 
  CodePropertiesDB::lock()
  {
    _issendParameterCancel<int>();
    _issendParameterCancel<double>();
    _issendParameterCancel<string>();
    _issendLockStatus = 1;
    _issendLock();
  }
  
  /**
   * \brief unlock access to local parameters from a distant code  
   *
   */

  void 
  CodePropertiesDB::unLock()
  {
    _issendLockStatus = 0;
    _issendParameters<int>();
    _issendParameters<double>();
    _issendParameters<string>();
    _issendLock();
  }

  /**
   * \brief Parameters non blocking send
   *
   */
  
  template < typename T > 
  void 
  CodePropertiesDB::_issendParameters()
  {
    if (!_issendLockStatus) {

      typedef typename map <string, T>::iterator IteratorMapT;
      typedef map <string, CodeProperties * >::iterator IteratorMapAppli;
  
      string nameType = typeid(T).name();

      int tagInit = 0;
      for (int k = 0; k < nameType.size(); k++)
        tagInit += nameType[k];

      map <string, T> *locCtrlParamPt; 
      _locCodeProperties->ctrlParamGet(&locCtrlParamPt);

      map <string, T> &locCtrlParam = *locCtrlParamPt;

      int NLocCtrlParam = locCtrlParamPt->size();

      int locFirstRank = _locCodeProperties->firstRankGet();

      const MPI_Comm& intraComm  = _locCodeProperties->intraCommGet();
      const MPI_Comm& globalComm = _locCodeProperties->globalCommGet();

      int intraCommSize = -1;
      int currentRank   = -1;

      MPI_Comm_rank(intraComm, &currentRank);
      MPI_Comm_size(intraComm, &intraCommSize);

      if (currentRank == 0) {

        int nAppli = _codePropertiesDB.size() + 1;

        //
        // Kill Existing messages

        _issendParameterCancel<T>();

        //
        // Define name parameters buffer

        int lStrings = locCtrlParam.size(); // To take into account '\0'
        for (IteratorMapT p1 = locCtrlParam.begin(); 
                          p1 != locCtrlParam.end(); 
                              p1++) {
          lStrings += p1->first.size();
        }
        
        if (_issendNameBuffs[typeid(T).name()].size() < lStrings) 
          _issendNameBuffs[typeid(T).name()].resize(lStrings);

        lStrings = 0;
        for (IteratorMapT p1 = locCtrlParam.begin(); 
                          p1 != locCtrlParam.end(); 
                          p1++) {
          for (int k = 0; k <  p1->first.size(); k++) 
            _issendNameBuffs[typeid(T).name()][lStrings + k] =  p1->first[k];
            
          lStrings += p1->first.size();
 
          _issendNameBuffs[typeid(T).name()][lStrings] = '\0';
            
          lStrings += 1;
        }

        //
        // Define values parameters buffer

        size_t sizeT = sizeof(T);
        int lBuff    = _issendLBuffGet(locCtrlParam);

        if (_issendValBuffs[typeid(T).name()].size() < lBuff) 
          _issendValBuffs[typeid(T).name()].resize(lBuff);

        int ival = 0 ;
        unsigned char *issendValBuff = &(_issendValBuffs[typeid(T).name()][0]);

        ival=0;

// FIXME: Il doit y avoir une double boucle IteratorMapT
        for (IteratorMapT p1 = locCtrlParam.begin(); 
                          p1 != locCtrlParam.end(); 
                          p1++) {

          _issendBuffCopy(issendValBuff, locCtrlParam);
        }

        //
        // Issend

        for (IteratorMapAppli p = _codePropertiesDB.begin(); 
                              p != _codePropertiesDB.end(); 
                              p++) {

          int distFirstRank = p->second->firstRankGet();

          vector<MPI_Request> & vectRequest = 
            *(_issendMPIrequest[typeid(T).name()][p->first]);

          //
          // Parameters name


          int tag = tagInit+'_'+'n'+'a'+'m'+'e'+'_'+'p'+'a'+'r'+'a'+'m'+ 
            nAppli * locFirstRank + distFirstRank;

          MPI_Issend(&(_issendNameBuffs[typeid(T).name()][0]), 
                     lStrings, MPI_CHAR, distFirstRank, tag,
                     globalComm, &(vectRequest[0]));

          //
          // Values

          tag = tagInit+'_'+'v'+'a'+'l'+'_'+'p'+'a'+'r'+'a'+'m'+ 
            nAppli * locFirstRank + distFirstRank;

          MPI_Issend(issendValBuff, 
                     lBuff, 
                     MPI_UNSIGNED_CHAR, 
                     distFirstRank,
                     tag,
                     globalComm,
                     &(vectRequest[1]));

        }
      }
    }
  }

  /**
   * \brief Cancel current parameters non blocking send
   *
   */
  
  template < typename T > 
  void 
  CodePropertiesDB::_issendParameterCancel()
  {
    typedef map <string, CodeProperties * >::iterator IteratorMapAppli;

    int nAppli         = _codePropertiesDB.size() + 1;
    int locFirstRank = _locCodeProperties->firstRankGet();
    
    const MPI_Comm& intraComm  = _locCodeProperties->intraCommGet();
    const MPI_Comm& globalComm = _locCodeProperties->globalCommGet();
    
    int locCommSize = -1;
    int currentRank = -1;
    
    MPI_Comm_rank(intraComm, &currentRank);
    MPI_Comm_size(intraComm, &locCommSize);

    if (currentRank == 0) {

      for (IteratorMapAppli p  = _codePropertiesDB.begin(); 
                            p != _codePropertiesDB.end(); 
                            p++) {

        //
        // Kill Existing messages

        int flag;
        vector<MPI_Request> & vectRequest = 
          *(_issendMPIrequest[typeid(T).name()][p->first]);

        for(int k = 0; k < _nIssend; k++) {

          MPI_Test(&(vectRequest[k]), &flag, MPI_STATUS_IGNORE);

          if (!flag) {
            MPI_Cancel(&(vectRequest[k]));
            MPI_Request_free(&(vectRequest[k]));
          }
        }
      }
    }
  }

  /**
   * \brief Compute the buffer length to send int values
   *
   * \param[in] locCtrlParam  Local control parameters  
   *
   * \return   Buffer length
   *
   */

  int 
  CodePropertiesDB::_issendLBuffGet
  (
   map <string, int> &locCtrlParam
  )
  {
    return sizeof(int) * locCtrlParam.size();
  }

  /**
   * \brief Compute the buffer length to send double values
   *
   * \param[in] locCtrlParam  Local control parameters  
   *
   * \return   Buffer length
   *
   */

  int
  CodePropertiesDB::_issendLBuffGet
  (
   map <string, double> &locCtrlParam
  )
  {
    return sizeof(double) * locCtrlParam.size();
  }

  /**
   * \brief Compute the buffer length to send string values
   *
   * \param[in] locCtrlParam  Local control parameters  
   *
   * \return   Buffer length
   *
   */

  int
  CodePropertiesDB::_issendLBuffGet
  (
   map <string, string> &locCtrlParam
  )
  {
    typedef map <string, string>::iterator IC;
    int lBuff = 0;

    for (IC p  = locCtrlParam.begin(); 
            p != locCtrlParam.end(); 
            p++) {
      lBuff += p->second.size() + 1;
    }

    return lBuff;
  }

  /**
   * \brief Copy int values into the buffer
   *
   * \param[in] buff            Buffer  
   * \param[in] locCtrlParam  Local control parameters  
   *
   */

  void
  CodePropertiesDB::_issendBuffCopy
  (
   unsigned char      *buff,
   map <string, int> &locCtrlParam
   )
  {
    typedef typename map <string, int>::iterator IteratorMapT;
    int ival = 0;

    for (IteratorMapT p1  = locCtrlParam.begin(); 
                      p1 != locCtrlParam.end(); 
                      p1++) {

      unsigned char * valChar =  (unsigned char *) (&(p1->second));

      for (int k = 0; k < sizeof(int); k++)
        buff[ival++] = valChar[k];
    }
  }


  /**
   * \brief Copy double values into the buffer
   *
   * \param[in] buff            Buffer  
   * \param[in] locCtrlParam  Local control parameters  
   *
   */

  void
  CodePropertiesDB::_issendBuffCopy
  (
   unsigned char        *buff,
   map <string, double> &locCtrlParam
   )
  {
    typedef typename map <string, double>::iterator IteratorMapT;
    int ival = 0;

    for (IteratorMapT p1  = locCtrlParam.begin(); 
                      p1 != locCtrlParam.end(); 
                      p1++) {

      unsigned char * valChar =  (unsigned char *) (&(p1->second));

      for (int k = 0; k < sizeof(double); k++)
        buff[ival++] = valChar[k];
    }
  }

    
  /**
   * \brief Copy double values into the buffer
   *
   * \param[in] buff            Buffer  
   * \param[in] locCtrlParam  Local control parameters  
   *
   */

  void
  CodePropertiesDB::_issendBuffCopy
  (
   unsigned char        *buff,
   map <string, string> &locCtrlParam
   )
  {
    typedef map <string, string>::iterator IC;
    int lBuff = 0;

    for (IC p  = locCtrlParam.begin(); 
            p != locCtrlParam.end(); 
            p++) {
      unsigned char * valChar =  (unsigned char *) (&(p->second[0]));
      for (int k = 0; k < p->second.size(); k++) {
        buff[lBuff++] = valChar[k];
      }
      buff[lBuff++] = '\0';
    }
  }

  /**
   * \brief Define int values from the received buffer
   *
   * \param[in] lNames          Total length of parameters names 
   * \param[in] buff            Buffer  
   * \param[in] distCtrlParam   Distant control parameters  
   *
   */

  void
  CodePropertiesDB::_irecvBuffCopy
  (
   const int            lNames,
   const unsigned char *buff,
   map <string, int>   &distCtrlParam
  )
  {
    int next = 0;
    int k    = 0;
    const int *valParamPt = (int *) buff; 
    
    while (next < lNames) {
      string key = &(_recvNameBuff[0]) + next;
      distCtrlParam[key] = valParamPt[k++];
      next += key.size()+1;
    }
  }


  /**
   * \brief Define double values from the received buffer
   *
   * \param[in] lNames          Total length of parameters names 
   * \param[in] buff            Buffer  
   * \param[in] distCtrlParam   Distant control parameters  
   *
   */
  
  void
  CodePropertiesDB::_irecvBuffCopy
  (
   const int            lNames,
   const unsigned char  *buff,
   map <string, double> &distCtrlParam
  )
  {
    int next = 0;
    int k    = 0;
    const double *valParamPt = (double *) buff; 
    
    while (next < lNames) {
      string key = &(_recvNameBuff[0]) + next;
      distCtrlParam[key] = valParamPt[k++];
      next += key.size()+1;
    }
  }
  
  /**
   * \brief Define string values from the received buffer
   *
   * \param[in] lNames          Total length of parameters names 
   * \param[in] buff            Buffer  
   * \param[in] distCtrlParam   Distant control parameters  
   *
   */
  
  void
  CodePropertiesDB::_irecvBuffCopy
  (
   const int            lNames,
   const unsigned char  *buff,
   map <string, string> &distCtrlParam
  )
  {
    int next  = 0;
    int next1 = 0;
    const char *valParamPt = (char *) buff; 
    
    while (next < lNames) {
      string key = (char *) (&(_recvNameBuff[0]) + next);
      string val = (char *) (valParamPt + next1);
 
      distCtrlParam[key] = val;
      next  += key.size()+1;
      next1 += val.size()+1;
    }
  }

  /**
   * \brief Get the parameter lock status of a distant code
   *
   * \param [in]  Code name 
   *
   */

  template < typename T > 
  void 
  CodePropertiesDB::_irecvParameters
  (
   const string &codeName
  )
  {
    MPI_Status status;

    string nameType = typeid(T).name();
    
    int tagInit = 0;
    for (int k = 0; k < nameType.size(); k++)
      tagInit += nameType[k];

    //
    // Parameters exchange between first ranks

    typedef typename map <string, T>::iterator IteratorMap;
    typedef map <string, CodeProperties * >::iterator IteratorMapAppli;

    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);
    if (p == _codePropertiesDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", codeName.c_str());

    map <string, T> *distCtrlParamPt; 
    p->second->ctrlParamGet(&distCtrlParamPt);

    map <string, T> &distCtrlParam = *distCtrlParamPt;
    
    int locFirstRank = _locCodeProperties->firstRankGet();
    int locLastRank  = _locCodeProperties->lastRankGet();
    
    const MPI_Comm& intraComm  = _locCodeProperties->intraCommGet();
    const MPI_Comm& globalComm = _locCodeProperties->globalCommGet();

    int locCommSize = -1;
    int currentRank = -1;
    
    MPI_Comm_rank(intraComm, &currentRank);
    MPI_Comm_size(intraComm, &locCommSize);

    int nAppli = _codePropertiesDB.size() + 1;
    
    int flag;
    int NDistControlParameters;
    int lStrings;

    int lBuff;

    //
    // Wait lock status is unlocked

    _lockStatusGet(p->first);
    while(_distLockStatus[p->first])
      _lockStatusGet(p->first);

    //
    // Get value

    if (currentRank == 0) {
      
      int distFirstRank = p->second->_firstRank;

      int tag = tagInit+'_'+'n'+'a'+'m'+'e'+'_'+'p'+'a'+'r'+'a'+'m'+ 
        nAppli * distFirstRank + locFirstRank;

      int mpiError = MPI_ERR_COMM; // Initialize to an error go into "while" loop

      while(mpiError != MPI_SUCCESS) {
        
        flag = 0;
        int num = 0;
        const int numMax = 1000; 
        while(!flag && (num < numMax)) {
          mpiError = MPI_Iprobe(distFirstRank, 
                                tag, 
                                globalComm, 
                                &flag, 
                                &status);
          num++;
        }
        if (mpiError != MPI_SUCCESS) {
          cout << "erreur mpi_iprobe " << locFirstRank << 
            " " << distFirstRank << " " << endl << flush; 
          continue;
        }

        if ((num > 1) && (num < numMax))
          cout << "Warning : pb iprobe latency " << locFirstRank << 
            " " << distFirstRank << " "  << num << endl << flush; 

        if (flag) {
          
          distCtrlParam.clear();
            
          mpiError = MPI_Get_count(&status, 
                                   MPI_CHAR, 
                                   &lStrings);


          if (mpiError != MPI_SUCCESS)  {
            cout << "erreur mpi_get_count" << locFirstRank << 
              " " << distFirstRank << " " << endl << flush; 
            continue;
          }

          if (lStrings > _recvNameBuff.size())
            _recvNameBuff.resize(lStrings);

          mpiError = MPI_Recv(&(_recvNameBuff[0]), 
                              lStrings, 
                              MPI_CHAR, 
                              distFirstRank,
                              tag,
                              globalComm, 
                              MPI_STATUS_IGNORE);

          if (mpiError != MPI_SUCCESS) {
            cout << "erreur mpi_recv1" << locFirstRank << 
              " " << distFirstRank << " " << endl << flush; 
            continue;
          }

          //
          // Receive values
          
          tag = tagInit+'_'+'v'+'a'+'l'+'_'+'p'+'a'+'r'+'a'+'m'+ 
            nAppli * distFirstRank + locFirstRank;

          mpiError = MPI_Probe(distFirstRank,
                               tag,
                               globalComm, 
                               &status);

          if (mpiError != MPI_SUCCESS) {
            cout << "erreur mpi_probe" << locFirstRank << " " << distFirstRank << " " << endl << flush; 
            continue;
          }

          mpiError = MPI_Get_count(&status, 
                                   MPI_CHAR, 
                                   &lBuff);

          if (mpiError != MPI_SUCCESS) {
            cout << "erreur mpi_gget_count2" << locFirstRank << " " << distFirstRank << " " << endl << flush; 
            continue;
          }

          if (lBuff > _recvValBuff.size())
            _recvValBuff.resize(lBuff);

          unsigned char *recvValBuff = (unsigned char *) &(_recvValBuff[0]);
          

          mpiError = MPI_Recv(recvValBuff, 
                              lBuff, 
                              MPI_UNSIGNED_CHAR, 
                              distFirstRank, 
                              tag,
                              globalComm, 
                              MPI_STATUS_IGNORE);

          if (mpiError != MPI_SUCCESS)  {
            cout << "erreur mpi_recv2" << locFirstRank << " " << distFirstRank << " " << endl << flush; 
            continue;
          }
        }
      }
    }

    if (locCommSize > 1) {
      MPI_Bcast(&flag, 
                1, 
                MPI_INT, 
                0, 
                intraComm);

      if (flag) {
        MPI_Bcast(&NDistControlParameters, 
                  1, 
                  MPI_INT, 
                  0, 
                  intraComm);
          
        MPI_Bcast(&lStrings, 
                  1, 
                  MPI_INT, 
                  0, 
                  intraComm);
          
        if (lStrings > _recvNameBuff.size())
          _recvNameBuff.resize(lStrings);
          
        MPI_Bcast(&(_recvNameBuff[0]), 
                  lStrings, 
                  MPI_CHAR, 
                  0, 
                  intraComm);
          
        MPI_Bcast(&lBuff, 
                  1, 
                  MPI_INT, 
                  0, 
                  intraComm);

        if (lBuff > _recvValBuff.size())
          _recvValBuff.resize(lBuff);

        unsigned char *recvValBuff = (unsigned char *) &(_recvValBuff[0]);

         MPI_Bcast(recvValBuff, 
                  lBuff, 
                  MPI_UNSIGNED_CHAR, 
                  0, 
                  intraComm);
      }
    }

    if (flag) {

      unsigned char *recvValBuff = (unsigned char *) &(_recvValBuff[0]);

      _irecvBuffCopy(lStrings,
                     recvValBuff, 
                     distCtrlParam);

    }
  }

  /**
   * \brief Return number of local parameters
   *
   * \return  Number of parameters
   *
   */
  
  template < typename T > 
  int 
  CodePropertiesDB::locCtrlParamNGet() const
  {
    return _locCodeProperties->ctrlParamNGet<T>();
 }

  /**
   * \brief Return of local parameters
   *
   * \return  List of parameters
   *
   */
  
  template < typename T > 
  char **
  CodePropertiesDB::locCtrlParamListGet() const
  {
    return _locCodeProperties->ctrlParamListGet<T>();
  }
  
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
  CodePropertiesDB::locCtrlParamIs
  (
     const string &name
   ) const
  {
    return _locCodeProperties->ctrlParamIs<T>(name);
  }


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
  CodePropertiesDB::distCtrlParamNGet
  (
   const string &codeName
   ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);

    if (p == _codePropertiesDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' code not found \n", codeName.c_str());

    return p->second->ctrlParamNGet<T>();
  }
  

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
  CodePropertiesDB::distCtrlParamListGet
  (
   const string &codeName
   ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);

    if (p == _codePropertiesDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' code not found \n", codeName.c_str());

    return p->second->ctrlParamListGet<T>();
  }

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
  CodePropertiesDB::distCtrlParamIs
  (
   const string &codeName,
   const string &name
   ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);

    if (p == _codePropertiesDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                 "'%s' code not found \n", codeName.c_str());

    return p->second->ctrlParamIs<T>(name);
  }
  

} // namespace cwipi

#endif /* __CODE_PROPERTIES_DB_H__ */
