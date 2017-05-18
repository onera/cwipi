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

#include <cstring>
#include <ctime>
#include <cstdlib>
#include <list>
#include <map>

#include <bftc_mem.h>

#include <fvmc_parall.h>

#include "cwipi_config.h"
#include "codePropertiesDB.hxx"
#include "codePropertiesDB_i.hxx"

using namespace std;

namespace cwipi {

  const int CodePropertiesDB::_nIssend = 2;

  /**
   * \brief Default constructor.
   *
   */

  CodePropertiesDB::CodePropertiesDB()
    : _codePropertiesDB(*(new map <string, CodeProperties * > ())),
      _locCodePropertiesDB(*(new map <string, CodeProperties * > ()))
  {
  }

  /**
   * \brief Destructor.
   *
   */

  CodePropertiesDB::~CodePropertiesDB()
  {
#if defined(DEBUG) && 0
    cout << "destroying CodePropertiesDB." << endl;
#endif

    if (!_codePropertiesDB.empty()) {
      typedef map <string, CodeProperties * >::iterator CI;
      for (CI p = _codePropertiesDB.begin();
           p != _codePropertiesDB.end(); p++) {
        if (p->second != NULL)
          delete p->second;
      }
      _codePropertiesDB.clear();

    }
    delete &_codePropertiesDB;
    delete &_locCodePropertiesDB;

    if (!_issendMPIrequest.empty()) {
      typedef map <string, vector<MPI_Request> * >::iterator CI;
      typedef map <string, map <string, vector<MPI_Request> * > >::iterator CI2;
      typedef map <string, map <string, map <string, vector<MPI_Request> * > > >::iterator CI1;

      for (CI1 p1  = _issendMPIrequest.begin();
               p1 != _issendMPIrequest.end(); 
               p1++) {
        for (CI2 p2  = p1->second.begin();
                 p2 != p1->second.end(); 
                 p2++) {
          for (CI p  = p2->second.begin();
                  p != p2->second.end(); 
                  p++) {
            if (p->second != NULL) {
              for (int k = 0; k < p->second->size(); k++)
                if ((*(p->second))[k] != MPI_REQUEST_NULL)
                  MPI_Request_free(&(*(p->second))[k]);
              delete p->second;
            }
          }
          p2->second.clear();
        }
        p1->second.clear();
      }
      _issendMPIrequest.clear();
    }


    typedef map <string, map <string, MPI_Request > >::iterator CI33;
    typedef map <string, MPI_Request >::iterator CI3;
    for (CI33 p33  = _issendLockMPIrequest.begin();
              p33 != _issendLockMPIrequest.end(); 
              p33++) {
      
      for (CI3 p3  = p33->second.begin();
               p3 != p33->second.end(); 
               p3++) {
      
        if (p3->second != MPI_REQUEST_NULL) {
          MPI_Request_free(&(p3->second));
        }
      }
      p33->second.clear();
    }
    _issendLockMPIrequest.clear();

    _lockStatus.clear();

    typedef map <string, string >::iterator CI4;
    typedef map <string, map <string, string > >::iterator CI41;
    for (CI41 p41  = _issendNameBuffs.begin();
              p41 != _issendNameBuffs.end(); 
              p41++) {
      for (CI4 p4  = p41->second.begin();
               p4 != p41->second.end(); 
               p4++) {
        p4->second.clear();
      }
      p41->second.clear();
    }
    _issendNameBuffs.clear();

    typedef map <string, vector<unsigned char > >::iterator CI5;
    typedef map <string, map <string, vector<unsigned char > > >::iterator CI51;
    for (CI51 p51  = _issendValBuffs.begin();
              p51 != _issendValBuffs.end(); 
              p51++) {
      for (CI5 p5  = p51->second.begin();
               p5 != p51->second.end(); 
               p5++) {
        p5->second.clear();
      }
      p51->second.clear();
    }
    _issendValBuffs.clear();

    _recvValBuff.clear();
    _recvNameBuff.clear();
    
  }

  
  /**
   * \brief MPI Communicator Initialization.
   *
   * This function builds the current code intra-communicator from
   * the current name and the MPI communicator containing all processes of
   * all codes.
   *
   * \param [in]  n_codes      Number of codes on the current rank
   * \param [in]  code_names   Codes names on the current rank
   * \param [in]  globalComm   MPI communicator containing all processes 
   *                           of all codes
   *
   * \return                   Current code intra-communicator
   *
   */

  void 
  CodePropertiesDB::init
  (
  const MPI_Comm     globalComm,
  const int          n_codes,
  const char**       code_names, 
  const CWP_Status_t is_coupled_rank,
  MPI_Comm           *intra_comms       
 )
  {

    // Initialize MPI
    // --------------
    
    int currentRank;
    int globalCommSize;
    
    _globalComm = globalComm;

    int flag;

    MPI_Initialized(&flag);
    if (!flag)
      bftc_error(__FILE__, __LINE__, 0, "MPI is not initialized\n");

    MPI_Comm_rank(globalComm, &currentRank);
    MPI_Comm_size(globalComm, &globalCommSize);

    // Search codes
    // ------------

    int index = 0;
    int properties[3]; /* The first property is the number of rank codes
                        * The second is the length of concatenated code names 
                        * The third is coupled rank state */

    properties[0] = n_codes;
    properties[1] = 0;
    for (int i = 0; i < n_codes; i++) {
      properties[1] += strlen(code_names[i]) + 1;
    }
    properties[2] = is_coupled_rank;

    int *allProperties = new int[3*globalCommSize];

    MPI_Allgather(&properties,
                  3,
                  MPI_INT,
                  allProperties,
                  3,
                  MPI_INT,
                  globalComm);

    char *concatenateNames = new char[properties[1]];
    char *_concatenateNames = concatenateNames;
    for (int i = 0; i < n_codes; i++) {
      strcpy(_concatenateNames, code_names[i]);
      _concatenateNames += strlen(code_names[i]) + 1;
    }

    int *iproc = new int[globalCommSize + 1];
    int *codesLengthName = new int[globalCommSize];
    iproc[0] = 0;
    for(int i = 0; i < globalCommSize; i++) {
      codesLengthName[i] = allProperties[3*i+1];
      iproc[i+1] = allProperties[3*i+1] + iproc[i];
    }

    int totalLength = iproc[globalCommSize];
    char *mergeNames = new char[totalLength];
    
    MPI_Allgatherv((void*) concatenateNames,
                   properties[1],
                   MPI_CHAR,
                   mergeNames,
                   codesLengthName,
                   iproc,
                   MPI_CHAR,
                   globalComm);

    delete[] concatenateNames;
    delete[] codesLengthName;

    map < string, vector < int > * >  coupledRankCode;
    map < string, vector < int > * >  rankCode;
    
    int id = 1;
     
    _isLocalCodeRootrank = false;
    
    for (int irank = 0; irank < globalCommSize; irank++) {

      for (int k = 0; k < allProperties[3*irank]; k++) {
        assert(index <= totalLength);

        const char *ptCurrentName = mergeNames + index;

        string currentName = string(ptCurrentName);

        typedef map <string, CodeProperties *>::iterator Iterator;

        Iterator p = _codePropertiesDB.find(currentName);
        
        if (p == _codePropertiesDB.end()) {

          if (!_isLocalCodeRootrank) {
            _isLocalCodeRootrank = currentRank == irank;
          }
          CodeProperties *currentCodeProperties =
            new CodeProperties(currentName, id++, irank, 
                               currentRank == irank, globalComm);

          coupledRankCode[currentName] = new vector <int> ();
          coupledRankCode[currentName]->reserve(globalCommSize);
          rankCode[currentName] = new vector <int> ();
          rankCode[currentName]->reserve(globalCommSize);

          pair<string, CodeProperties *>
            newPair(currentName, currentCodeProperties);

          _codePropertiesDB.insert(newPair);

          if (currentRank == irank) {
            _locCodePropertiesDB.insert(newPair);
          }

          _lockStatus[string(currentName)] = 0;
        }

        else {
          
          if (currentRank == irank) {
            if (_locCodePropertiesDB.find(currentName) == _locCodePropertiesDB.end()) {
               _locCodePropertiesDB[currentName] = p->second;
            }
          }
        
        }
        
        if (allProperties[3*irank+2]) {
          coupledRankCode[currentName]->push_back(irank);
        }
        rankCode[currentName]->push_back(irank);

        if (currentRank == irank) {
          _codePropertiesDB[currentName]->isCoupledRankset(is_coupled_rank);
        }

        index += currentName.size() + 1;
        assert(index <= totalLength);
      }
    }

    typedef map <string, CodeProperties *>::iterator IteratorCP;
     
    for (IteratorCP p1  = _locCodePropertiesDB.begin(); 
                    p1 != _locCodePropertiesDB.end(); 
                    p1++) {
      for (IteratorCP p2  = _codePropertiesDB.begin(); 
                      p2 != _codePropertiesDB.end(); 
                      p2++) {
        if ((p1->second != p2->second) && 
            (p1->second->_rootRankInGlobalComm != p2->second->_rootRankInGlobalComm)) {
        
          _issendLockMPIrequest[p1->first][p2->first] = MPI_REQUEST_NULL;
          _issendMPIrequest[p1->first][typeid(int).name()][p2->first] =  
            new vector<MPI_Request>(_nIssend, MPI_REQUEST_NULL);

          _issendMPIrequest[p1->first][typeid(double).name()][p2->first] =  
            new vector<MPI_Request>(_nIssend, MPI_REQUEST_NULL);

          _issendMPIrequest[p1->first][typeid(string).name()][p2->first] =  
            new vector<MPI_Request>(_nIssend, MPI_REQUEST_NULL);
        }
      }
    }
         
    delete [] allProperties;
    delete [] iproc;
    delete [] mergeNames;
    
    // Create intra code communicators
    // -------------------------------

    MPI_Group globalGroup;
    MPI_Comm_group (globalComm, &globalGroup);
    
    typedef map <string, vector < int > * >::iterator Iterator;

    for (Iterator p = rankCode.begin(); 
                  p != rankCode.end(); p++){
      const int *_ranks = &((*p->second)[0]);
      int _n_ranks = p->second->size();
      MPI_Group_incl (globalGroup, _n_ranks, _ranks, 
                     &(_codePropertiesDB[p->first]->_intraGroup));
      MPI_Comm_create (globalComm, 
                       _codePropertiesDB[p->first]->_intraGroup,
                       &(_codePropertiesDB[p->first]->_intraComm));
      _codePropertiesDB[p->first]->_rootRankInGlobalComm = _ranks[0];
      delete p->second;
    }
    
    rankCode.clear();

    // Create intra coupled code group
    // -------------------------------

    for (Iterator p = coupledRankCode.begin(); 
                  p != coupledRankCode.end(); p++){
      const int *_ranks = &((*p->second)[0]);
      int _n_ranks = p->second->size();
      MPI_Group_incl (globalGroup, _n_ranks, _ranks, 
                     &(_codePropertiesDB[p->first]->_intraCoupledGroup));
      delete p->second;
    }

    coupledRankCode.clear();

    _tagLockStatusBase = 1000;
    _tagParameterBase = _tagLockStatusBase + _codePropertiesDB.size();
        
    for (int i = 0; i < n_codes; i++) {
      const string &nameStr = code_names[i];
      intra_comms[i] = _codePropertiesDB[nameStr]->_intraComm;
      _issendLock(nameStr);
    }

  }

  /**
   * \brief Dump properties.  
   *
   */

  void 
  CodePropertiesDB::dump()
  {
    typedef map <string, CodeProperties *>::iterator Iterator;

    bftc_printf("\nCode properties\n\n");
    for (Iterator p = _codePropertiesDB.begin(); 
                  p != _codePropertiesDB.end(); p++){
      _irecvParameters<int>(p->first);
      _irecvParameters<double>(p->first);
      _irecvParameters<string>(p->first);
      p->second->dump();
      bftc_printf("\n");
    }
    bftc_printf_flush();
  }

  /**
   * \brief Get the parameter lock status of a distant code
   *
   * \param [in]  Code name 
   *
   */

  void 
  CodePropertiesDB::_lockStatusGet
  (
   const string &codeName
  )
  {

    //
    // Find codeName in local codes
    
    const map <string, CodeProperties * >::iterator p1 = 
      _codePropertiesDB.find(codeName);
    
    if (p1 == _codePropertiesDB.end()) {
      bftc_error(__FILE__, __LINE__, 0,
              "'%s' code not found \n", codeName.c_str());
    }

    const map <string, CodeProperties * >::iterator p = 
      _locCodePropertiesDB.find(codeName);

    if (p == _locCodePropertiesDB.end()) {
    
      const MPI_Comm& intraComm  = p->second->intraCommGet();
      const MPI_Comm& globalComm = p->second->globalCommGet();

      int locCommSize = -1;
      int flag;

      MPI_Comm_size(intraComm, &locCommSize);

      if (_isLocalCodeRootrank) {
        int distFirstRank = p->second->_rootRankInGlobalComm;

        int tag = _tagLockStatusBase + p->second->_id; 

        MPI_Iprobe(distFirstRank, 
                   tag, 
                   globalComm, 
                   &flag, 
                   MPI_STATUS_IGNORE);

        if (flag) {

          MPI_Recv(&(_lockStatus[p->first]),
                   1, 
                   MPI_INT, 
                   distFirstRank, 
                   tag,
                   globalComm, 
                   MPI_STATUS_IGNORE);
        }     
      }

      if (locCommSize > 1) {
        MPI_Bcast(&flag, 
                  1, 
                  MPI_INT, 
                  0, 
                  intraComm);

        if (flag) {
          MPI_Bcast(&(_lockStatus[p->first]), 
                    1, 
                    MPI_INT, 
                    0, 
                    intraComm);
        }
      }
    }
  }
   
  /**
   * \brief Lock status non blocking sending
   *
   * \param [in]  codeName  Local code name
   *
   */
  
  void 
  CodePropertiesDB::_issendLock
  (
   const string &codeName
  )
  {
    typedef map <string, CodeProperties * >::iterator IteratorMapCP;

    const IteratorMapCP p = _locCodePropertiesDB.find(codeName); 
    if (p == _locCodePropertiesDB.end()) {
      bftc_error(__FILE__, __LINE__, 0,
              "'%s' code not found \n", codeName.c_str());
    }
    
    CodeProperties * locCodeProperties = p->second;

    int globRootRank = locCodeProperties->_rootRankInGlobalComm;
    
    const MPI_Comm& intraComm  = locCodeProperties->intraCommGet();
    const MPI_Comm& globalComm = locCodeProperties->globalCommGet();
    
    int intraCommSize = -1;
    int intraRank = -1;
    
    MPI_Comm_rank(intraComm, &intraRank);
    MPI_Comm_size(intraComm, &intraCommSize);

    int globalCommSize = -1;
    int globalRank = -1;
    
    MPI_Comm_rank(globalComm, &globalRank);
    MPI_Comm_size(globalComm, &globalCommSize);

    //
    // Kill Existing messages

    if (globRootRank == globalRank) {

      for (IteratorMapCP p1  = _codePropertiesDB.begin(); 
                         p1 != _codePropertiesDB.end(); 
                         p1++) {

        if (locCodeProperties != p1->second) {
          int distRootRank = p1->second->_rootRankInGlobalComm;
          if (distRootRank != globRootRank) {
            if (_issendLockMPIrequest[p->first][p1->first] != MPI_REQUEST_NULL) {
              
              int flag;
              MPI_Test(&(_issendLockMPIrequest[p->first][p1->first]), 
                       &flag, 
                       MPI_STATUS_IGNORE);

              if (!flag) {
                MPI_Cancel(&(_issendLockMPIrequest[p->first][p1->first]));
                MPI_Request_free(&(_issendLockMPIrequest[p->first][p1->first]));
              }
            }
          }
        }
      }

      for (IteratorMapCP p1  = _codePropertiesDB.begin(); 
                         p1 != _codePropertiesDB.end(); 
                         p1++) {

        if (locCodeProperties != p1->second) {
          int distRootRank = p1->second->_rootRankInGlobalComm;
          int tag = _tagLockStatusBase + p->second->_id; 

          MPI_Issend(&_lockStatus[codeName], 
                     1, 
                     MPI_INT, 
                     distRootRank, 
                     tag,
                     globalComm, 
                     &(_issendLockMPIrequest[p->first][p1->first]));
        }
      }
    }
  }

}
