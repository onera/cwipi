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
    : _distCodePropertiesDB(*(new map <string, CodeProperties * > ())),
      _locCodeProperties(NULL)
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
    if (_locCodeProperties != NULL)
      delete _locCodeProperties;

    if (!_distCodePropertiesDB.empty()) {
      typedef map <string, CodeProperties * >::iterator CI;
      for (CI p = _distCodePropertiesDB.begin();
           p != _distCodePropertiesDB.end(); p++) {
        if (p->second != NULL)
          delete p->second;
      }
      _distCodePropertiesDB.clear();

    }
    delete &_distCodePropertiesDB;

    if (!_issendMPIrequest.empty()) {
      typedef map <string, vector<MPI_Request> * >::iterator CI;
      typedef map <string, map <string, vector<MPI_Request> * > >::iterator CI2;

      for (CI2 p2  = _issendMPIrequest.begin();
               p2 != _issendMPIrequest.end(); 
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
      _issendMPIrequest.clear();
    }

    typedef map <string, MPI_Request >::iterator CI3;
    for (CI3 p3  = _issendLockMPIrequest.begin();
             p3 != _issendLockMPIrequest.end(); 
             p3++) {
      if (p3->second != MPI_REQUEST_NULL)
        MPI_Request_free(&(p3->second));
    }
    _issendLockMPIrequest.clear();

    _distLockStatus.clear();

    typedef map <string, string >::iterator CI4;
    for (CI4 p4  = _issendNameBuffs.begin();
             p4 != _issendNameBuffs.end(); 
             p4++) {
      p4->second.clear();
    }
    _issendNameBuffs.clear();

    typedef map <string, vector<unsigned char > >::iterator CI5;
    for (CI5 p5  = _issendValBuffs.begin();
             p5 != _issendValBuffs.end(); 
             p5++) {
      p5->second.clear();
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
   * \param [in]  name         Current code name
   * \param [in]  globalComm   MPI communicator containing all processes 
   *                           of all codes
   *
   * \return                   Current code intra-communicator
   *
   */

  MPI_Comm 
  CodePropertiesDB::init
  (
   const char* name, 
   const MPI_Comm globalComm
  )
  {

    // Initialize MPI
    // --------------

    int currentRank;
    int globalCommSize;
    int color = 0;

    {
      int flag;

      MPI_Initialized(&flag);
      if (!flag)
        bftc_error(__FILE__, __LINE__, 0, "MPI is not initialized\n");

      MPI_Comm_rank(globalComm, &currentRank);
      MPI_Comm_size(globalComm, &globalCommSize);

    }

    // Search codes
    // ------------

    {
      int j = 0;
      int index = 0;
      int totalLength = 0;
      int nameLength = 0;

      string currentString = "";

      nameLength = strlen(name) + 1;

      MPI_Allreduce (&nameLength, &totalLength, 1, MPI_INT, MPI_SUM,
                     globalComm);

      char *mergeNames = new char[totalLength];

      int *namesLength = new int[globalCommSize];
      int *iproc = new int[globalCommSize];

      MPI_Allgather(&nameLength,
                    1,
                    MPI_INT,
                    namesLength,
                    1,
                    MPI_INT,
                    globalComm);

      iproc[0] = 0;
      for(int i = 1; i < globalCommSize; i++)
        iproc[i] = namesLength[i-1] + iproc[i-1];

      MPI_Allgatherv((void*) const_cast <char*> (name),
                     nameLength,
                     MPI_CHAR,
                     mergeNames,
                     namesLength,
                     iproc,
                     MPI_CHAR,
                     globalComm);

      delete[] iproc;
      delete[] namesLength;

      for (int irank = 0; irank < globalCommSize; irank++) {

        assert(index <= totalLength);

        const char *ptCurrentName = mergeNames + index;
        string currentName(ptCurrentName);

        if (currentString != currentName) {

          CodeProperties *currentCodeProperties =
            new CodeProperties(currentName, globalComm);

          if (!strcmp(currentName.c_str(), name)) {
            _locCodeProperties = currentCodeProperties;
            color = j+1;
          }
          else {
            pair<string, CodeProperties *>
              newPair(string(currentName), currentCodeProperties);

            pair<map<string, CodeProperties *>::iterator, bool>
              p = _distCodePropertiesDB.insert(newPair);

            if (!p.second)
              bftc_error(__FILE__, __LINE__, 0,
                        "The MPI ranks range is not continous or\n"
                        "There are two codes with the same name '%s'  \n", 
                         currentName.c_str());

            //
            // issend int param MPI request storage

            _issendMPIrequest[typeid(int).name()][string(currentName)] =  
              new vector<MPI_Request>(_nIssend, MPI_REQUEST_NULL);

            _issendMPIrequest[typeid(double).name()][string(currentName)] =  
              new vector<MPI_Request>(_nIssend, MPI_REQUEST_NULL);

            _issendMPIrequest[typeid(string).name()][string(currentName)] =  
              new vector<MPI_Request>(_nIssend, MPI_REQUEST_NULL);

            //
            // issend lock MPI request storage

            _issendLockMPIrequest[string(currentName)] = MPI_REQUEST_NULL;
            _distLockStatus[string(currentName)] = 0;

          }

          currentCodeProperties->firstRankSet(irank);

          if (currentString != "") {
            if (!strcmp(currentString.c_str(), name))
              _locCodeProperties->lastRankSet(irank-1);
            else
              _distCodePropertiesDB[currentString]->lastRankSet(irank-1);
          }

          currentString = currentName;

          j += 1;
        }

        if (currentString != "") {
          if (!strcmp(currentString.c_str(), name))
            _locCodeProperties->lastRankSet(globalCommSize-1);
          else
            _distCodePropertiesDB[currentString]->lastRankSet(globalCommSize-1);
        }

        index += currentName.size() + 1;
        assert(index <= totalLength);
      }

      delete [] mergeNames;

      // Create current code communicator
      // ---------------------------------------

      MPI_Comm intraComm = MPI_COMM_NULL ;
      MPI_Comm_split(globalComm, color, currentRank, &intraComm);
      _locCodeProperties->intraCommSet(intraComm);

      _issendLockStatus = 0;
      _issendLock();

      return intraComm;
    }
  }

  /**
   * \brief Dump properties.  
   *
   */

  void 
  CodePropertiesDB::dump()
  {
    bftc_printf("\nLoc code properties\n\n");
    _locCodeProperties->dump();

    typedef map <string, CodeProperties *>::iterator Iterator;

    bftc_printf("\nDistant code properties\n\n");
    for (Iterator p = _distCodePropertiesDB.begin(); 
                  p != _distCodePropertiesDB.end(); p++){
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
    int locFirstRank = _locCodeProperties->firstRankGet();
    int nAppli         = _distCodePropertiesDB.size() + 1;
    
    const MPI_Comm& intraComm  = _locCodeProperties->intraCommGet();
    const MPI_Comm& globalComm = _locCodeProperties->globalCommGet();

    int locCommSize = -1;
    int currentRank   = -1;
    
    MPI_Comm_rank(intraComm, &currentRank);
    MPI_Comm_size(intraComm, &locCommSize);

    const map <string, CodeProperties * >::iterator p = 
      _distCodePropertiesDB.find(codeName);
 
   if (p == _distCodePropertiesDB.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", codeName.c_str());

    int flag;
    _distLockStatus[p->first];

    if (currentRank == 0) {
      
      int distFirstRank = p->second->_firstRank;

      int tag = 'l'+'o'+'c'+'k'+'_'+'s'+'t'+'a'+'t'+'u'+'s'+ 
        nAppli * distFirstRank + locFirstRank;
          
      MPI_Iprobe(distFirstRank, 
                 tag, 
                 globalComm, 
                 &flag, 
                 MPI_STATUS_IGNORE);

      if (flag) {

        //
        // Receive lock status
          
        MPI_Recv(&(_distLockStatus[p->first]),
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
        MPI_Bcast(&(_distLockStatus[p->first]), 
                  1, 
                  MPI_INT, 
                  0, 
                  intraComm);
      }
    }
  }
   
  /**
   * \brief Lock status non blocking sending
   *
   */
  
  void 
   CodePropertiesDB::_issendLock()
  {
    typedef map <string, CodeProperties * >::iterator IteratorMapAppli;

    int nAppli         = _distCodePropertiesDB.size() + 1;
    int locFirstRank = _locCodeProperties->firstRankGet();
    
    const MPI_Comm& intraComm  = _locCodeProperties->intraCommGet();
    const MPI_Comm& globalComm = _locCodeProperties->globalCommGet();
    
    int locCommSize = -1;
    int currentRank = -1;
    
    MPI_Comm_rank(intraComm, &currentRank);
    MPI_Comm_size(intraComm, &locCommSize);

    //
    // Kill Existing messages

    if (currentRank == 0) {

      for (IteratorMapAppli p = _distCodePropertiesDB.begin(); 
                            p != _distCodePropertiesDB.end(); 
                            p++) {

        int distFirstRank = p->second->_firstRank;

        if (_issendLockMPIrequest[p->first] != MPI_REQUEST_NULL) {
          
          int flag;
          MPI_Test(&(_issendLockMPIrequest[p->first]), 
                   &flag, 
                   MPI_STATUS_IGNORE);

          if (!flag) {
            MPI_Cancel(&(_issendLockMPIrequest[p->first]));
            MPI_Request_free(&(_issendLockMPIrequest[p->first]));
          }
        }
      }

      for (IteratorMapAppli p = _distCodePropertiesDB.begin(); 
                            p != _distCodePropertiesDB.end(); 
                            p++) {

        int distFirstRank = p->second->_firstRank;
        int tag = 'l'+'o'+'c'+'k'+'_'+'s'+'t'+'a'+'t'+'u'+'s'+ 
          nAppli * locFirstRank + distFirstRank;

        MPI_Issend(&_issendLockStatus, 
                   1, 
                   MPI_INT, 
                   distFirstRank, 
                   tag,
                   globalComm, 
                   &(_issendLockMPIrequest[p->first]));
      }
    }
   }

}
