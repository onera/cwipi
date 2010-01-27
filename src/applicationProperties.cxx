///////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009  ONERA
//
// Description :
//    This class describes application properties :
//        - intra mpi communicator (localComm)
//        - beginning and end rank in the inter mpi communicator (globalComm)
//        - control parameters storage
//
///////////////////////////////////////////////////////////////////////////////
#include "applicationProperties.hxx"

namespace cwipi
{

  ApplicationProperties::ApplicationProperties(std::string &name,
                                               const MPI_Comm globalComm)

    : _name(name), _globalComm(globalComm),
    _intControlParameters(*(new std::map <std::string, int>())),
    _doubleControlParameters(*(new std::map <std::string, double>())),
    _stringControlParameters(*(new std::map <std::string, std::string>()))
  {
    _localComm = NULL;
    _beginningRank = -999;
    _endRank = -999;
  }

  ApplicationProperties::ApplicationProperties(const ApplicationProperties& other)
    : _name(other._name), _globalComm(other._globalComm), _localComm(other._localComm),
      _beginningRank(other._beginningRank), _endRank(other._endRank),
      _intControlParameters(other._intControlParameters),
      _doubleControlParameters(other._doubleControlParameters),
      _stringControlParameters(other._stringControlParameters)
  {
  }

  void ApplicationProperties::dump()
  {
    bft_printf("'%s' properties\n",_name.c_str());
    bft_printf("  - Ranks in global MPI_comm : %i <= ranks <= %i \n",
               _beginningRank,
               _endRank);
    bft_printf("  - Int Control Parameter :\n");

    typedef std::map <std::string, int>::iterator Iterator1;
    for (Iterator1 p = _intControlParameters.begin(); p != _intControlParameters.end(); p++)
      bft_printf("   * '%s' : %i\n", p->first.c_str(), p->second);
    bft_printf("  - Double Control Parameter :\n");

    typedef std::map <std::string, double>::iterator Iterator2;
    for (Iterator2 p = _doubleControlParameters.begin(); p != _doubleControlParameters.end(); p++)
      bft_printf("   * '%s' : %12.5e\n", p->first.c_str(), p->second);

    bft_printf("  - String Control Parameter :\n");

    typedef std::map <std::string, std::string>::iterator Iterator3;
    for (Iterator3 p = _stringControlParameters.begin(); p != _stringControlParameters.end(); p++)
      bft_printf("   * '%s' : '%s'\n", p->first.c_str(), p->second.c_str());

    bft_printf_flush();
  }

  ApplicationProperties::~ApplicationProperties()
  {

    if (!_intControlParameters.empty())
      _intControlParameters.clear();
    delete &_intControlParameters;
    if (!_doubleControlParameters.empty())
      _doubleControlParameters.clear();
    delete &_doubleControlParameters;
    if (!_stringControlParameters.empty())
      _stringControlParameters.clear();
    delete &_stringControlParameters;
    if (_localComm != NULL)
      MPI_Comm_free(&_localComm);

  }
}
