
#include <cstring>

#include <bft_mem.h>

#include <fvm_parall.h>

#include "applicationPropertiesDataBase.hxx"
#include "applicationPropertiesDataBase_i.hxx"

namespace couplings {

  ApplicationPropertiesDataBase::ApplicationPropertiesDataBase()
    : _distantApplicationPropertiesDataBase(*(new std::map <std::string, ApplicationProperties * > ())),
      _localApplicationProperties(NULL)
  {
  }

  ApplicationPropertiesDataBase::~ApplicationPropertiesDataBase()
  {
    // BUG SB: pas bug : j'ai change le std::cout en bft_printf
    bft_printf( "destroying ApplicationPropertiesDataBase.\n" );
    if (_localApplicationProperties != NULL)
      delete _localApplicationProperties;

    if (!_distantApplicationPropertiesDataBase.empty()) {
      typedef std::map <std::string, ApplicationProperties * >::iterator CI;
      for (CI p = _distantApplicationPropertiesDataBase.begin();
           p != _distantApplicationPropertiesDataBase.end(); p++) {
        if (p->second != NULL)
          delete p->second;
      }
      _distantApplicationPropertiesDataBase.clear();

    }
    delete &_distantApplicationPropertiesDataBase;
  }

  MPI_Comm  ApplicationPropertiesDataBase::init(const char* applicationName,
                                                const couplings_mpi_ranks_for_coupling_t mpi_rank,
                                                const MPI_Comm globalComm)
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
        bft_error(__FILE__, __LINE__, 0, "MPI is not initialized\n");

      MPI_Comm_rank(globalComm, &currentRank);
      MPI_Comm_size(globalComm, &globalCommSize);
    }

    // Search applications
    // -------------------

    {
      int j = 0;
      int index = 0;
      int totalLength = 0;
      int nameLength = 0;

      std::string currentString = "";
      char *mergeNames = NULL;
      int currentMpiRanks = 0;

      nameLength = strlen(applicationName) + 1;
      MPI_Allreduce (&nameLength, &totalLength, 1, MPI_INT, MPI_SUM,
                     globalComm);

      BFT_MALLOC(mergeNames, totalLength, char) ;


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

      MPI_Allgatherv((void*) const_cast <char*> (applicationName),
                     nameLength,
                     MPI_CHAR,
                     mergeNames,
                     namesLength,
                     iproc,
                     MPI_CHAR,
                     globalComm);

      // Exchange of the status "one master"
      int *mpi_ranks = new int[globalCommSize];
      MPI_Allgather( (void*)&mpi_rank,
                     1,
                     MPI_INT,
                     mpi_ranks,
                     1,
                     MPI_INT,
                     globalComm );
      std::cout << "rangs echanges" << std::flush;
      for( int i=0; i<globalCommSize; i++) std::cout << mpi_ranks[i]<<" " << std::flush; std::cout << std::endl << std::flush;

      delete[] iproc;
      delete[] namesLength;

      for (int irank = 0; irank < globalCommSize; irank++) {

        assert(index <= totalLength);

        const char *ptCurrentName = mergeNames + index;
        std::string currentName(ptCurrentName);
        std::cout << irank << currentName.c_str() << " " << currentString.c_str() << std::endl << std::flush;

        if (currentString != currentName) {

          ApplicationProperties *currentApplicationProperties =
            new ApplicationProperties(currentName, globalComm);

          if (!strcmp(currentName.c_str(), applicationName)) {
            _localApplicationProperties = currentApplicationProperties;
            color = j+1;
          }
          else {
            std::pair<std::string, ApplicationProperties *>
              newPair(std::string(currentName), currentApplicationProperties);

            std::pair<std::map<std::string, ApplicationProperties *>::iterator, bool>
              p = _distantApplicationPropertiesDataBase.insert(newPair);

            if (!p.second)
              bft_error(__FILE__, __LINE__, 0,
                        "The MPI ranks range is not continous or\n"
                        "There are two applications with the same name '%s'  \n", currentName.c_str());
          }

          currentApplicationProperties->setBeginningRank(irank);
          std::cout << currentName.c_str() << "added" << irank << mpi_ranks[irank] << std::endl << std::flush;
          if( mpi_ranks[irank] == COUPLINGS_MPI_RANKS_ONLY_MASTER ) {
            if( ! strcmp( currentName.c_str(), applicationName ) ){
              _localApplicationProperties->setEndRank(irank);
              std::cout << currentName.c_str() << ":master:local:1:" <<irank <<std::endl<< std::flush;}
            else{
              _distantApplicationPropertiesDataBase[currentName]->setEndRank(irank);
              std::cout << currentName.c_str() << ":master:distant:1:" <<irank <<std::endl<< std::flush;}

          }
          if (currentString != "" && currentMpiRanks == COUPLINGS_MPI_RANKS_ALL_RANKS) {
            if (!strcmp(currentString.c_str(), applicationName)){
              _localApplicationProperties->setEndRank(irank-1);
              std::cout << currentString.c_str() << ":all:local:1:" <<irank-1 <<std::endl<< std::flush;}
            else{
              _distantApplicationPropertiesDataBase[currentString]->setEndRank(irank-1);
              std::cout << currentString.c_str() << ":all:distant:1:" <<irank-1 <<std::endl<< std::flush;}
          }
          currentString = currentName;
          currentMpiRanks = mpi_ranks[irank];

          j += 1;
        }

        if (currentString != "") { // BUG: toujours vrai... n'est-il-pas ? pourrait aussi etre sorti de la boucle ?
        if(currentMpiRanks == COUPLINGS_MPI_RANKS_ALL_RANKS) {
          if (!strcmp(currentString.c_str(), applicationName)) {
            _localApplicationProperties->setEndRank(globalCommSize-1);
            std::cout << currentString.c_str() << ":all:local:2:" <<globalCommSize-1 <<std::endl<< std::flush;
        }
          else{
            _distantApplicationPropertiesDataBase[currentString]->setEndRank(globalCommSize-1);
            std::cout << currentString.c_str() << ":all:distant:2:" <<globalCommSize-1 <<std::endl<< std::flush;
          }
        }
        }

        index += currentName.size() + 1;
        assert(index <= totalLength);
      }
std::cout << "fin" << std::endl << std::flush;
//      for (int irank = 0; irank < globalCommSize; irank++) {
//        const char *ptCurrentName = mergeNames + index;
//                std::string currentName(ptCurrentName);
//int a;
//           if (currentString != currentName) {
//           if (!strcmp(currentString.c_str(), applicationName))
//                    a = _localApplicationProperties->getEndRank();
//                  else
//                    a = _distantApplicationPropertiesDataBase[currentString]->getEndRank();
//           std::cout << a << std::flush;
//           currentString = currentName;
//         }
//      }

      if (mergeNames != NULL)
        BFT_FREE (mergeNames);

      // Create current application communicator
      // ---------------------------------------

      MPI_Comm localComm = MPI_COMM_NULL;
      MPI_Comm_split(globalComm, color, currentRank, &localComm); // BUG: was MPI_COMM_WORLD
      _localApplicationProperties->setLocalComm(localComm);
      std::cout << "commsplit effectue." << std::endl << std::flush;

      MPI_Comm restrictedLocalComm = localComm;
      if( mpi_rank == COUPLINGS_MPI_RANKS_ONLY_MASTER ) {
        int color = MPI_UNDEFINED, localCurrentRank;
        MPI_Comm_rank( localComm, &localCurrentRank );
        if( localCurrentRank == 0 ) color = 0;
        MPI_Comm_split(localComm, color, localCurrentRank, &restrictedLocalComm);
      }
      fvm_parall_set_mpi_comm(restrictedLocalComm);
      std::cout << "fvm set mpi comm effectue : " << fvm_parall_get_rank() << " " << fvm_parall_get_size() << std::endl << std::flush;
      return localComm;
    }
  }


  void ApplicationPropertiesDataBase::mergeParameters(const std::string &applicationName)
  {
    _mergeIntParameters(applicationName);
    _mergeDoubleParameters(applicationName);
  }


  void ApplicationPropertiesDataBase::_mergeIntParameters(const std::string &applicationName)
  {

    typedef std::map <std::string, int>::iterator Iterator;

    MPI_Status status;

    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bft_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());

    std::map <std::string, int> &localControlParameters   = _localApplicationProperties->_intControlParameters;
    std::map <std::string, int> &distantControlParameters = p->second->_intControlParameters;

    int NLocalControlParameters = localControlParameters.size();
    int NDistantControlParameters = distantControlParameters.size();

    int localBeginningRank = _localApplicationProperties->getBeginningRank();
    int localEndRank       = _localApplicationProperties->getEndRank();
    const MPI_Comm& localComm    = _localApplicationProperties->getLocalComm();
    const MPI_Comm& globalComm   = _localApplicationProperties->getGlobalComm();

    int distantBeginningRank = p->second->_beginningRank;
    int distantEndRank       = p->second->_endRank;

    //
    // Clear distant parameters copy
    // BUG SB: en double ??
    distantControlParameters.clear();
    distantControlParameters.clear();

    //
    // Check that all local pocesses are the same parameter values

    int localCommSize = -1;
    int currentRank = -1;

    MPI_Comm_rank(localComm, &currentRank);
    MPI_Comm_size(localComm, &localCommSize);

    if (localCommSize > 1) {
      typedef std::map <std::string, int>::iterator Iterator;

      for (Iterator p = localControlParameters.begin(); p != localControlParameters.end(); p++) {
        int value             = p->second;
        const char *paramName = p->first.c_str();
        int nameSize    = p->first.size();

        if (currentRank == 0) {
          for (int irank = 1; irank < localCommSize; irank++){

            int   distantNameSize = 0;
            char *distantParamName = NULL;
            int   distantValue;

            MPI_Recv(&distantNameSize, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantNameSize != nameSize)
              bft_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");

            distantParamName = new char[distantNameSize+1];
            MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, irank, 0, localComm, &status);
            if (strcmp(distantParamName, paramName))
              bft_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");
            delete[] distantParamName;

            MPI_Recv(&distantValue, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantValue != value)
              bft_error(__FILE__, __LINE__, 0,
                        "Different values for '%s' parameter for the rank 0 and the rank '%i'\n", paramName, irank);
          }
        }

        else {
          MPI_Send(&nameSize, 1, MPI_INT, 0, 0, localComm);
          MPI_Send(const_cast <char* > (paramName), nameSize+1, MPI_CHAR, 0, 0, localComm);
          MPI_Send(&value, 1, MPI_INT, 0, 0, localComm);
        }
      }
    }

    //
    // parameters exchange between beginning ranks

    if (currentRank == 0 ) {

      MPI_Sendrecv(&NLocalControlParameters,   1, MPI_INT, distantBeginningRank, 100,
                   &NDistantControlParameters, 1, MPI_INT, distantBeginningRank, 100,
                   globalComm, &status);

      int NParameterMax = MAX(NLocalControlParameters, NDistantControlParameters);

      Iterator p = localControlParameters.begin();

      for (int i = 0; i < NParameterMax; i++ ) {

        if (i >= NDistantControlParameters) {

          int value             = p->second;
          const char *paramName = p->first.c_str();
          int nameSize          = p->first.size();

          MPI_Send(&nameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm);
          MPI_Send(const_cast <char *> (paramName), nameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm);
          MPI_Send(&value, 1, MPI_INT, distantBeginningRank, 0, globalComm);

        }

        else if (p == localControlParameters.end()){

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          int   distantValue;

          MPI_Recv(&distantNameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);

          distantParamName = new char[distantNameSize+1];
          MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm, &status);

          MPI_Recv(&distantValue, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;

          delete[] distantParamName;

        }

        else {

          int value             = p->second;
          const char *paramName = p->first.c_str();
          int nameSize          = p->first.size();

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          int   distantValue;

          MPI_Sendrecv(&nameSize       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantNameSize, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantParamName = new char[distantNameSize+1];

          MPI_Sendrecv(const_cast <char*> (paramName), nameSize+1,        MPI_CHAR, distantBeginningRank, 0,
                       distantParamName              , distantNameSize+1, MPI_CHAR, distantBeginningRank, 0,
                       globalComm, &status);

          MPI_Sendrecv(&value       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantValue, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;

          delete[] distantParamName;
        }

        if (p != localControlParameters.end())
          p++;
      }
    }

    //
    // Local beginning rank send parameters to other local ranks

    if (localCommSize > 1) {

      int size = 0;
      if (currentRank == 0)
        size =  distantControlParameters.size();

      MPI_Bcast(&size, 1, MPI_INT, 0, localComm);

      Iterator p;
      if (currentRank == 0)
        p = distantControlParameters.begin();

      for (int i = 0 ; i < size; i++) {

        int value;
        char *paramName = NULL;
        int nameSize;

        if (currentRank == 0) {
          value     = p->second;
          paramName = const_cast <char *> (p->first.c_str());
          nameSize  = p->first.size();
          p++;
        }

        MPI_Bcast(&nameSize, 1, MPI_INT, 0, localComm);


        if (currentRank != 0)
          paramName = new char[nameSize+1];

        MPI_Bcast(paramName, nameSize+1, MPI_CHAR, 0, localComm);

        MPI_Bcast(&value, 1, MPI_INT, 0, localComm);

        if (currentRank != 0) {
          distantControlParameters[std::string(paramName)] = value;
          delete[] paramName;
        }
      }
    }
  }


  void ApplicationPropertiesDataBase::_mergeDoubleParameters(const std::string &applicationName)
  {

    typedef std::map <std::string, double>::iterator Iterator;

    MPI_Status status;
    //std::string applicationNameStr = applicationName;

    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bft_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());

    std::map <std::string, double> &localControlParameters   = _localApplicationProperties->_doubleControlParameters;
    std::map <std::string, double> &distantControlParameters = p->second->_doubleControlParameters;

    int NLocalControlParameters = localControlParameters.size();
    int NDistantControlParameters = distantControlParameters.size();

    int localBeginningRank = _localApplicationProperties->getBeginningRank();
    int localEndRank       = _localApplicationProperties->getEndRank();
    const MPI_Comm& localComm    = _localApplicationProperties->getLocalComm();
    const MPI_Comm& globalComm   = _localApplicationProperties->getGlobalComm();

    int distantBeginningRank = p->second->_beginningRank;
    int distantEndRank       = p->second->_endRank;

    //
    // Clear distant parameters copy

    distantControlParameters.clear();
    distantControlParameters.clear();

    //
    // Check that all local pocesses are the same parameter values

    int localCommSize = -1;
    int currentRank = -1;

    MPI_Comm_rank(localComm, &currentRank);
    MPI_Comm_size(localComm, &localCommSize);

    if (localCommSize > 1) {
      typedef std::map <std::string, double>::iterator Iterator;

      for (Iterator p = localControlParameters.begin(); p != localControlParameters.end(); p++) {
        double value             = p->second;
        const char *paramName = p->first.c_str();
        int nameSize    = p->first.size();

        if (currentRank == 0) {
          for (int irank = 1; irank < localCommSize; irank++){

            int   distantNameSize = 0;
            char *distantParamName = NULL;
            double   distantValue;

            MPI_Recv(&distantNameSize, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantNameSize != nameSize)
              bft_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");

            distantParamName = new char[distantNameSize+1];
            MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, irank, 0, localComm, &status);
            if (strcmp(distantParamName, paramName))
              bft_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");
            delete[] distantParamName;

            MPI_Recv(&distantValue, 1, MPI_DOUBLE, irank, 0, localComm, &status);
            if (distantValue != value)
              bft_error(__FILE__, __LINE__, 0,
                        "Different values for '%s' parameter for the rank 0 and the rank '%i'\n", paramName, irank);
          }
        }

        else {
          MPI_Send(&nameSize, 1, MPI_INT, 0, 0, localComm);
          MPI_Send(const_cast <char* > (paramName), nameSize+1, MPI_CHAR, 0, 0, localComm);
          MPI_Send(&value, 1, MPI_DOUBLE, 0, 0, localComm);
        }
      }
    }

    //
    // parameters exchange between beginning ranks

    if (currentRank == 0 ) {

      MPI_Sendrecv(&NLocalControlParameters,   1, MPI_INT, distantBeginningRank, 0,
                   &NDistantControlParameters, 1, MPI_INT, distantBeginningRank, 0,
                   globalComm, &status);

      int NParameterMax = MAX(NLocalControlParameters, NDistantControlParameters);

      Iterator p = localControlParameters.begin();

      for (int i = 0; i < NParameterMax; i++ ) {

        if (i >= NDistantControlParameters) {

          double value             = p->second;
          const char *paramName = p->first.c_str();
          int nameSize          = p->first.size();

          MPI_Send(&nameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm);
          MPI_Send(const_cast <char *> (paramName), nameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm);
          MPI_Send(&value, 1, MPI_DOUBLE, distantBeginningRank, 0, globalComm);

        }

        else if (p == localControlParameters.end()){

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          double   distantValue;

          MPI_Recv(&distantNameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);

          distantParamName = new char[distantNameSize+1];
          MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm, &status);

          MPI_Recv(&distantValue, 1, MPI_DOUBLE, distantBeginningRank, 0, globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;
          delete[] distantParamName;

        }

        else {

          double value             = p->second;
          const char *paramName = p->first.c_str();
          int nameSize          = p->first.size();

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          double   distantValue;

          MPI_Sendrecv(&nameSize       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantNameSize, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantParamName = new char[distantNameSize+1];
          MPI_Sendrecv(const_cast <char*> (paramName), nameSize+1,        MPI_CHAR, distantBeginningRank, 0,
                       distantParamName              , distantNameSize+1, MPI_CHAR, distantBeginningRank, 0,
                       globalComm, &status);

          MPI_Sendrecv(&value       , 1, MPI_DOUBLE, distantBeginningRank, 0,
                       &distantValue, 1, MPI_DOUBLE, distantBeginningRank, 0,
                       globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;
          delete[] distantParamName;

        }

        if (p != localControlParameters.end())
          p++;
      }
    }

    //
    // Local beginning rank send parameters to other local ranks

    if (localCommSize > 1) {

      int size = 0;
      if (currentRank == 0)
        size =  distantControlParameters.size();

      MPI_Bcast(&size, 1, MPI_INT, 0, localComm);

      Iterator p;
      if (currentRank == 0)
        p = distantControlParameters.begin();

      for (int i = 0 ; i < size; i++) {
        double value;
        char *paramName = NULL;
        int nameSize;

        if (currentRank == 0) {
          value     = p->second;
          paramName = const_cast <char *> (p->first.c_str());
          nameSize  = p->first.size();
          p++;
        }

        MPI_Bcast(&nameSize, 1, MPI_INT, 0, localComm);


        if (currentRank != 0)
          paramName = new char[nameSize+1];

        MPI_Bcast(paramName, nameSize+1, MPI_CHAR, 0, localComm);

        MPI_Bcast(&value, 1, MPI_DOUBLE, 0, localComm);

        if (currentRank != 0) {
          distantControlParameters[std::string(paramName)] = value;
          delete[] paramName;
        }
      }
    }
  }

  void ApplicationPropertiesDataBase::dump()
  {
    bft_printf("\nLocal application properties\n\n");
    _localApplicationProperties->dump();

    typedef std::map <std::string, ApplicationProperties *>::iterator Iterator;

    bft_printf("\nDistant application properties\n\n");
    for (Iterator p = _distantApplicationPropertiesDataBase.begin(); p != _distantApplicationPropertiesDataBase.end(); p++){
      p->second->dump();
      bft_printf("\n");
    }
    bft_printf_flush();
  }

}
