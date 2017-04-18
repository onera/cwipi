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

#include <cassert>
#include <cmath>
#include <cstring>
#include <sstream>

#include <bftc_error.h>
#include <bftc_file.h>

#include <fvmc_parall.h>

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "mesh.hxx"
#include "codeProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "cwp.h"

#include "factory.hpp"
#include "geometry.hxx"
#include "communication.hxx"


/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) &&  !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

using namespace std;

namespace cwipi {

  typedef Factory<Communication, CWP_Comm_t> FC;
  typedef Factory<Geometry, CWP_Geom_t> FG;
  typedef Factory<Support, CWP_Support_t> FS;

  
//TODO: Nouveau Coupling a faire  
//  Coupling::Coupling
//  (
//   const string                &cplId,
//   const CWP_Comm_t           cplType,
//   const CodeProperties        &localCodeProperties,
//   const CodeProperties        &coupledCodeProperties,
//   const CWP_Geom_t           geomAlgo,
//   const CWP_Support_t        supportType,
//   const int                    nPart,
//   const CWP_Displacement_t  movingStatus,
//   const CWP_Freq_t           recvFreqType
//   )
//  :_cplId(cplId),
//   _communication(*(FC::getInstance().CreateObject(cplType))),
//   _localCodeProperties(localCodeProperties),
//   _coupledCodeProperties(coupledCodeProperties),
//   _geometry(*(FG::getInstance().CreateObject(geomAlgo))),
//
////!!!!!!!!!
////!!!!!   A poursuivre ici
////!!!!!!!
//   
// _entitiesDim(entitiesDim),_tolerance(tolerance), _solverType(solverType),
// _outputFormat(outputFormat), _outputFormatOption(outputFormatOption),
// _fvmWriter(NULL), _outputFrequency(outputFrequency), _name(name),
// _couplingType(couplingType),
// _tmpDistantFieldsIssend(*new map<int, vector<double> * > ()),
// _tmpLocalFieldsIrecv(*new map<int, const double * > ()),
// _tmpExchangeNameIrecv(*new  map<int, const char * > ()),
// _tmpStrideIrecv(*new  map<int, int > ()),
// _tmpTimeStepIrecv(*new  map<int, int > ()),
// _tmpTimeValueIrecv(*new  map<int, double > ()),
// _tmpFieldNameIrecv(*new  map<int, const char * > ())
//
//{
//  _tmpVertexField = NULL;
//  _tmpDistantField = NULL;
//  _supportMesh = NULL;
//  _couplingComm = MPI_COMM_NULL;
//  _coupledCodeNRankCouplingComm = -1;
//  _coupledCodeBeginningRankCouplingComm = -1;
//  _mergeComm = MPI_COMM_NULL;
//  _fvmComm = MPI_COMM_NULL;
//  _interpolationFct = NULL;
//  _interpolationFct_f = NULL;
//  _toLocate = true;
//  _isCoupledRank = false;
//  
//
//  //
//  // Create coupling comm
//
//  _createCouplingComm();
//
//  //
//  //
//
//  _locationToDistantMesh = new LocationToDistantMesh(_isCoupledRank,
//                                                     _couplingType,
//                                                     _localCodeProperties);
//
//  _locationToLocalMesh = new LocationToLocalMesh(_solverType,
//                                                 _tolerance,
//                                                 _couplingComm,
//                                                 _coupledCodeNRankCouplingComm,
//                                                 _coupledCodeBeginningRankCouplingComm,
//                                                 _isCoupledRank,
//                                                 _entitiesDim,
//                                                 _localCodeProperties,
//                                                 *_locationToDistantMesh);
//
//
//#ifndef NAN
//  bftc_printf("Warning : NAN macro is undefined -> receiving checking deactivation\n");
//#endif
//
//}
//
//
//Coupling::~Coupling()
//{
//#if defined(DEBUG) && 0
//  cout << "destroying '" << _name << "' coupling" << endl;
//#endif
//
//  typedef map <int, vector<double> * >::iterator Iterator1;
//  for (Iterator1 p = _tmpDistantFieldsIssend.begin();
//       p != _tmpDistantFieldsIssend.end(); p++) {
//    if (p->second != NULL)
//      delete p->second;
//  }
//
//  _tmpDistantFieldsIssend.clear();
//  delete &_tmpDistantFieldsIssend;
//
//  _tmpLocalFieldsIrecv.clear();
//  delete &_tmpLocalFieldsIrecv;
//
//  _tmpExchangeNameIrecv.clear();
//  delete &_tmpExchangeNameIrecv;
//
//  _tmpFieldNameIrecv.clear();
//  delete &_tmpFieldNameIrecv;
//
//  _tmpStrideIrecv.clear();
//  delete &_tmpStrideIrecv;
//
//  _tmpTimeStepIrecv.clear();
//  delete &_tmpTimeStepIrecv;
//
//  _tmpTimeValueIrecv.clear();
//  delete &_tmpTimeValueIrecv;
//
//  if (_isCoupledRank ) {
//
//    delete _tmpVertexField;
//
//    delete _tmpDistantField;
//
//    delete _supportMesh;
//
//    delete _locationToDistantMesh;
//
//    delete _locationToLocalMesh;
//
//    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
//    if (oldFVMComm != MPI_COMM_NULL)
//      MPI_Barrier(oldFVMComm);
//    fvmc_parall_set_mpi_comm(_fvmComm);
//
//    if (_fvmWriter != NULL)
//      fvmc_writer_finalize(_fvmWriter);
//
//    if (_fvmComm != MPI_COMM_NULL)
//      MPI_Barrier(_fvmComm);
//    fvmc_parall_set_mpi_comm(oldFVMComm);
//
//  }
//
//  if (_mergeComm != MPI_COMM_NULL)
//    MPI_Comm_free(&_mergeComm);
//
//  if (_couplingComm != MPI_COMM_NULL)
//    MPI_Comm_free(&_couplingComm);
//
//  if (_fvmComm != MPI_COMM_NULL)
//    MPI_Comm_free(&_fvmComm);
//
//}


} // namespace cwipi

