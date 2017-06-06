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
    
  }

  
  /**
   * \brief MPI Communicator Initialization.
   *
   * This function builds the current code intra-communicator from
   * the current name and the MPI communicator containing all processes of
   * all codes.
   *
   * \param [in]  n_codes         Number of codes on the current rank
   * \param [in]  code_names      Codes names on the current rank
   * \param [in]  is_coupled_rank Current rank is it a coupled rank
   * \param [in]  n_param_max     Maximum number of parameters
   * \param [in]  str_size_max    Maximum size for a string
   * \param [out] intra_coms      Current code intra-communicators
   *
   */

  void 
  CodePropertiesDB::init
  (
  const MPI_Comm     globalComm,
  const int          n_codes,
  const char**       code_names, 
  const CWP_Status_t is_coupled_rank,
  const int          n_param_max,
  const int          str_size_max,      
  MPI_Comm           *intra_comms       
 )
  {

    // Initialize MPI
    // --------------
    
    int currentRank;
    int globalCommSize;
    
    _globalComm = globalComm;
    
    _n_param_max = n_param_max;
    _str_size_max = str_size_max;

    int flag;

    MPI_Initialized(&flag);
    if (!flag) {
      bftc_error(__FILE__, __LINE__, 0, "MPI is not initialized\n");
    }
      
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
                               currentRank == irank, globalComm,
                               _n_param_max, _str_size_max);

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
        }

        else {
          
          if (currentRank == irank) {
            if (_locCodePropertiesDB.find(currentName) == _locCodePropertiesDB.end()) {
               _locCodePropertiesDB[currentName] = p->second;
               p->second->isLocalSet(true);
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
      
    for (IteratorCP p2  = _codePropertiesDB.begin(); 
                    p2 != _codePropertiesDB.end(); 
                    p2++) {
        
      if (p2->second->_rootRankInGlobalComm == currentRank) {
        MPI_Win_create(p2->second->_winGlobData, 
                       4 * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winGlob);
        
        p2->second->_winIntParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(p2->second->_winIntParamIdxNameData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamIdxName);
        
        p2->second->_winIntParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winIntParamNameData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamName);
        
        p2->second->_winIntParamValueData = 
            (int *) malloc (sizeof(int) * n_param_max);
        MPI_Win_create(p2->second->_winIntParamValueData, 
                       n_param_max * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamValue);
        
        p2->second->_winDoubleParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(p2->second->_winDoubleParamIdxNameData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamIdxName);
        
        p2->second->_winDoubleParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winDoubleParamNameData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamName);
        
        p2->second->_winDoubleParamValueData = 
            (double *) malloc (sizeof(double) * n_param_max);
        MPI_Win_create(p2->second->_winDoubleParamValueData, 
                       n_param_max * sizeof(double),
                       sizeof(double), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamValue);
        
        p2->second->_winStrParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(p2->second->_winStrParamIdxNameData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxName);
        
        p2->second->_winStrParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winStrParamNameData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamName);
        
        p2->second->_winStrParamIdxValueData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(p2->second->_winStrParamIdxValueData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxValue);
        
        p2->second->_winStrParamValueData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winStrParamValueData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamValue);
      }
      else {
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winGlob);

        p2->second->_winIntParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamIdxName);

        p2->second->_winIntParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max *str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamName);

        p2->second->_winIntParamValueData = 
            (int *) malloc (sizeof(int) * n_param_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamValue);

        p2->second->_winDoubleParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamIdxName);

        p2->second->_winDoubleParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max *str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamName);

        p2->second->_winDoubleParamValueData = 
            (double *) malloc (sizeof(double) * n_param_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(double), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamValue);

        p2->second->_winStrParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxName);
        
        p2->second->_winStrParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamName);
        
        p2->second->_winStrParamIdxValueData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxValue);
        
        p2->second->_winStrParamValueData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamValue);
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

    //TODO
    
    for (int i = 0; i < n_codes; i++) {
      const string &nameStr = code_names[i];
      intra_comms[i] = _codePropertiesDB[nameStr]->_intraComm;
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
    bftc_printf_flush();
  }

}
