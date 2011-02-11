#include <fvm_parall.h>

#include "singleton.hpp"
#include "couplingDataBase.hxx"
#include "couplingDataBase_i.hxx"
#include "applicationProperties.hxx"
#include "coupling.hxx"

namespace cwipi {

  CouplingDataBase::CouplingDataBase()
    :  _couplingDataBase(*new std::map <std::string, Coupling * > ())
  {
    _fvmComm = MPI_COMM_NULL;
  }

  CouplingDataBase::~CouplingDataBase()
  {
    typedef std::map <std::string, Coupling * >::iterator Iterator;
    for (Iterator p = _couplingDataBase.begin();
         p != _couplingDataBase.end(); p++) {
      if (p->second != NULL)
        delete p->second;
    }
    _couplingDataBase.clear();

    delete &_couplingDataBase;
  }

  void  CouplingDataBase::createCoupling(const std::string &name,
                                         const cwipi_coupling_type_t couplingType,
                                         const ApplicationProperties& localApplicationProperties,
                                         const ApplicationProperties& coupledApplicationProperties,
                                         const int entitiesDim,
                                         const double tolerance,
                                         const cwipi_solver_type_t solverType,
                                         const int    outputFrequency,
                                         const char  *outputFormat,
                                         const char  *outputFormatOption)
  {

    //
    // deactivate postProcessing if postprocessing is activated
    // for an other coupling with a different coupling type !!!
    // (fvm restriction)

    int newOutputFrequency = outputFrequency;
    int localCommSize;

    MPI_Comm_size(localApplicationProperties.getLocalComm(), &localCommSize);

    //
    // Set fvm parall comm

//     if (_fvmComm != MPI_COMM_NULL) {

//       int fvmCommSize;
//       MPI_Comm_size(_fvmComm, &fvmCommSize);

//       if ( localCommSize != 1) {

//         if ((couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING && fvmCommSize > 1) ||
//             (couplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING && fvmCommSize == 1)) {

//           bft_printf("Warning : Post-processing deactivation for '%s' coupling", name.c_str());
//           bft_printf("          Post-processing is activated for an other coupling type\n");
//           bft_printf("          To activate this Post-processing, deactivate the other\n");

//           newOutputFrequency = -1;
//         }
//       }
//     }

//     else {

//       if (couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {

//         int list = 0;
//         MPI_Group localGroup;
//         MPI_Group fvmGroup;

//         MPI_Comm_group(localApplicationProperties.getLocalComm(), &localGroup);
//         MPI_Group_incl(localGroup, 1, &list, &fvmGroup);
//         MPI_Comm_create(localApplicationProperties.getLocalComm(),
//                         fvmGroup,
//                         &_fvmComm);
//       }

//       else

//         MPI_Comm_dup(localApplicationProperties.getLocalComm(), &_fvmComm);

//       int titi;
//       int tutu;
//       MPI_Comm_size(_fvmComm, &tutu);
//       MPI_Comm_rank(_fvmComm, &titi);
//       std::cout << "fvmComm : " << tutu << " " << titi << std::endl;
//       fvm_parall_set_mpi_comm(_fvmComm);

//     }

    //
    // Create the new coupling

    Coupling *newCoupling = new Coupling(name,
                                         couplingType,
                                         localApplicationProperties,
                                         coupledApplicationProperties,
                                         entitiesDim,
                                         tolerance,
                                         solverType,
                                         newOutputFrequency,
                                         outputFormat,
                                         outputFormatOption);

    std::pair<std::string, Coupling* >
      newPair(std::string(name), newCoupling);

    std::pair<std::map<std::string, Coupling* >::iterator, bool>
      p = _couplingDataBase.insert(newPair);

    if (!p.second)
      bft_error(__FILE__, __LINE__, 0,
                "'%s' existing coupling\n", name.c_str());

  }

  void  CouplingDataBase::deleteCoupling(const std::string &name)
  {
    const std::map <std::string, Coupling * >::iterator p = _couplingDataBase.find(name);
    if (p == _couplingDataBase.end())
      bft_error(__FILE__, __LINE__, 0,
                "'%s' coupling not found \n", name.c_str());

    if (p->second != NULL)
      delete p->second;

    _couplingDataBase.erase(p);
  }
}
