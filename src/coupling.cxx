
#include <cassert>
#include <cmath>
#include <sstream>

#include <bft_error.h>
#include <bft_file.h>

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "mesh.hxx"
#include "applicationProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "coo_baryc.h"
#include "cwipi.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

extern "C" {
  void PROCF(callfortinterpfct, CALLFORTINTERPFCT)
  ( int *entities_dim,
      int *n_local_vertex,
      int *n_local_element,
      int *n_local_polhyedra,
      int *n_distant_point,
      double *local_coordinates,
      int *local_connectivity_index,
      int *local_connectivity,
      int *local_polyhedra_face_index,
      int *local_polyhedra_cell_to_face_connectivity,
      int *local_polyhedra_face_connectivity_index,
      int *local_polyhedra_face_connectivity,
      double *distant_points_coordinates,
      int *distant_points_location,
      int *distant_points_barycentric_coordinates_index,
      double *distant_points_barycentric_coordinates,
      int *data_dimension,
      int *solver_type,
      double *local_field,
      double *distant_field,
      void *ptFortranInterpolationFct
  );
}
namespace cwipi {

  //TODO: Faire une factory sur le type de couplage

Coupling::Coupling(const std::string& name,
                   const cwipi_coupling_type_t couplingType,
                   const ApplicationProperties& localApplicationProperties,
                   const ApplicationProperties& coupledApplicationProperties,
                   const int entitiesDim,
                   const double tolerance,
                   const cwipi_solver_type_t solverType,
                   const int    outputFrequency,
                   const char  *outputFormat,
                   const char  *outputFormatOption)
:_localApplicationProperties(localApplicationProperties),
 _coupledApplicationProperties(coupledApplicationProperties),
 _entitiesDim(entitiesDim),_tolerance(tolerance), _solverType(solverType),
 _outputFormat(outputFormat), _outputFormatOption(outputFormatOption),
 _fvmWriter(NULL), _outputFrequency(outputFrequency), _name(name),
 _nDistantpoint(0), _couplingType(couplingType)

{
  _tmpVertexField = NULL;
  _tmpDistantField = NULL;
  _supportMesh = NULL;
  _coordsPointsToLocate = NULL;
  _fvmLocator = NULL;
  _couplingComm = MPI_COMM_NULL;
  _coupledApplicationNRankCouplingComm = -1;
  _coupledApplicationBeginningRankCouplingComm = -1;
  _mergeComm = MPI_COMM_NULL;
  _fvmComm = MPI_COMM_NULL;
  _interpolationFct = NULL;
  _toLocate = true;
  _barycentricCoordinatesIndex = NULL;
  _barycentricCoordinates = NULL;
  _nNotLocatedPoint = 0;
  _nPointsToLocate = 0;
  _location = NULL; 
  _notLocatedPoint = NULL;
  _locatedPoint = NULL;
  _isCoupledRank = false;

  //
  // Create coupling comm
  
  _createCouplingComm ();


#ifndef NAN
  bft_printf("Warning : NAN macro is undefined -> receiving checking deactivation\n");
#endif
  
}

std::vector<double> &  Coupling::_extrapolate(double *cellCenterField)
{

  if (_tmpVertexField == NULL)
    _tmpVertexField = new  std::vector<double>(_supportMesh->getNVertex(),0.);

    std::vector<double> &vertexField = *_tmpVertexField;

    for (int i = 0; i < vertexField.size(); i++)
      vertexField[i] = 0.;

    // TODO: Faire l'allocation qu'une fois comme _tmpVertexField

    std::vector<double> volumeVertex(_supportMesh->getNVertex(),0.);

    assert(_supportMesh != NULL);

    const int nElts = _supportMesh->getNElts();
    const int nPoly = _supportMesh->getNPolyhedra();
    const int nStandardElt = nElts - nPoly;

    const int *eltConnectivityIndex = _supportMesh->getEltConnectivityIndex();
    const int *eltConnectivity      = _supportMesh->getEltConnectivity();
    const std::vector<double>& cellVolume       = _supportMesh->getVolume();

    for (int i = 0; i < nStandardElt; i++) {
      const int nEltVertex = eltConnectivityIndex[i+1] - eltConnectivityIndex[i];
      const int index = eltConnectivityIndex[i];
      for (int j = 0; j < nEltVertex; j++) {
        int vertex = eltConnectivity[index+j] - 1;
        volumeVertex[vertex] += cellVolume[i];
        vertexField[vertex]  += cellCenterField[i] * cellVolume[i];
      }
    }

    if (nPoly > 0) {
      std::vector<int> vertexPoly;
      vertexPoly.reserve(30);

      const int *polyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex();
      const int *polyhedraCellToFaceConnectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();
      const int *polyhedraFaceConnectivityIndex = _supportMesh->getPolyhedraFaceConnectivityIndex() ;
      const int *polyhedraFaceConnectivity = _supportMesh->getPolyhedraFaceConnectivity();

      for (int i = 0; i < nPoly; i++) {
        int nFacePolyhedra = polyhedraFaceIndex[i+1] - polyhedraFaceIndex[i];
        int faceIndex = polyhedraCellToFaceConnectivity[i];
        int nVertexFace = 0;
        for (int j = 0; j < nFacePolyhedra; j++) {
          int iface = polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
          int nVertexLocFace = polyhedraFaceConnectivityIndex[iface+1] - polyhedraFaceConnectivityIndex[iface];
          nVertexFace += nVertexLocFace;
          int vertexIndex = polyhedraFaceConnectivityIndex[iface];
          for (int k = 0; k < nVertexLocFace; k++) {
            if (vertexPoly.capacity() <= vertexPoly.size()) {
              int capacity = 2*vertexPoly.capacity();
              vertexPoly.reserve(capacity);
            }
            vertexPoly.push_back(polyhedraFaceConnectivity[vertexIndex+k]);
          }

          quickSort(&vertexPoly[0], 0, vertexPoly.size()-1, NULL);

          int ivertex = -1;

          for (int j = 0; j < vertexPoly.size(); j++) {
            if (ivertex < vertexPoly[j]) {
              ivertex = vertexPoly[j];
              volumeVertex[ivertex - 1] += cellVolume[i];
              vertexField[ivertex - 1]  += cellCenterField[i] * cellVolume[i];
            }
          }
        }
      }
    }

    for (int i = 0; i < _supportMesh->getNVertex(); i++) {
      assert(volumeVertex[i] > 0.);
      vertexField[i] /= volumeVertex[i];
    }

    return vertexField;
}

Coupling::~Coupling()
{
#if defined(DEBUG) && 0
  std::cout << "destroying '" << _name << "' coupling" << std::endl;
#endif  

  if (_isCoupledRank ) {

    delete _tmpVertexField;

    delete _tmpDistantField;

    delete _supportMesh;


  //
  // TODO: Recoder les coord bary 2D en c++
    if (_entitiesDim != 2) {
      delete[] _barycentricCoordinatesIndex;
      delete[] _barycentricCoordinates;
    }

    else {
      BFT_FREE(_barycentricCoordinatesIndex);
      BFT_FREE(_barycentricCoordinates);
    }
    
    fvm_parall_set_mpi_comm(_fvmComm);

    if (_fvmLocator != NULL)
      fvm_locator_destroy(_fvmLocator);

    if (_fvmWriter != NULL)
      fvm_writer_finalize(_fvmWriter);

    fvm_parall_set_mpi_comm(MPI_COMM_NULL);

  }
   
  if (_mergeComm != MPI_COMM_NULL)
    MPI_Comm_free(&_mergeComm); 

  if (_couplingComm != MPI_COMM_NULL)
    MPI_Comm_free(&_couplingComm); 
  
  if (_fvmComm != MPI_COMM_NULL)
    MPI_Comm_free(&_fvmComm);

  if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING && 
      !_isCoupledRank) {

    if (_locatedPoint != NULL) 
      delete []  _locatedPoint;
    
    if (_notLocatedPoint != NULL) 
      delete []  _notLocatedPoint ;
  }

}

void Coupling::_interpolate(double *referenceField,
                            std::vector<double>& interpolatedField,
                            const int stride)
{
  if (_solverType == CWIPI_SOLVER_CELL_CENTER)
    referenceField = &_extrapolate(referenceField)[0];

  switch(_entitiesDim) {

  case 1 :
    _interpolate1D(referenceField, interpolatedField, stride);
    break;

  case 2 :
    _interpolate2D(referenceField, interpolatedField, stride);
    break;

  case 3 :
    _interpolate3D(referenceField, interpolatedField, stride);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "'%i' bad entities dimension\n",_entitiesDim);
  }
}

void Coupling::_interpolate1D(double *referenceVertexField,
                              std::vector<double>& interpolatedField,
                              const int stride)
{
  const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
  const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
  const double *distantCoords   = fvm_locator_get_dist_coords(_fvmLocator);
  const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
  const int *eltsConnec        = _supportMesh->getEltConnectivity();
  const double *coords         = _supportMesh->getVertexCoords();

  for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
    int iel = distantLocation[ipoint] - 1;
    double coef1 = _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]];
    double coef2 = _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]+1];
    int pt1 = eltsConnecPointer[iel] - 1;
    int pt2 = eltsConnecPointer[iel+1] - 1;

    for (int k = 0; k < stride; k++)
      interpolatedField[stride*ipoint+k] = coef1 * referenceVertexField[stride*pt1+k] +
      coef2 * referenceVertexField[stride*pt2+k];
  }
}

void Coupling::_interpolate2D(double *vertexField,
                              std::vector<double>& interpolatedField,
                              const int stride)
{

  const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
  const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
  const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
  const int *eltsConnec        = _supportMesh->getEltConnectivity();

  for (int ipoint = 0; ipoint <nDistantPoint; ipoint++) {
    int iel = distantLocation[ipoint] - 1;
    int index = _barycentricCoordinatesIndex[ipoint];
    int nSom = _barycentricCoordinatesIndex[ipoint+1] - index;

    for (int k = 0; k < stride; k++)
      interpolatedField[stride*ipoint + k] = 0;
    for (int isom = 0; isom <  nSom; isom++) {
      for (int k = 0; k < stride; k++)
        interpolatedField[stride*ipoint+k] += vertexField[stride*(eltsConnec[eltsConnecPointer[iel]+isom]-1)+k]
                                                          *_barycentricCoordinates[index+isom];
    }
  }
}



void Coupling::_interpolate3D(double *vertexField,
                              std::vector<double>& interpolatedField,
                              const int stride)
{

  // TODO: Faire le calcul des coordonnees barycentriques pour les polyedres
  // TODO: Dans un premier temps faire le calcul pour les tétraèdres

  const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
  const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
  const double *distantCoords   = fvm_locator_get_dist_coords(_fvmLocator);
  const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
  const int *eltsConnec        = _supportMesh->getEltConnectivity();
  const int nStandardElt       = _supportMesh->getNElts() - _supportMesh->getNPolyhedra();
  const double *coords         = _supportMesh->getVertexCoords();
  double coeff[4];

  for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
    int iel = distantLocation[ipoint] - 1;
    if (iel < nStandardElt) {
      int index = eltsConnecPointer[iel];
      int nVertex = eltsConnecPointer[iel+1]-index;

      for (int k = 0; k < stride; k++) {
        double a[4][4] = {{0., 0., 0., 0.},
            {0., 0., 0., 0.},
            {0., 0., 0., 0.},
            {0., 0., 0., 0.}};
        double b[4] = {0., 0., 0., 0.};
        for (int i = 0; i < nVertex; i++) {
          int iVertex = eltsConnec[index+i]-1;
          double v_x = coords[3*iVertex];
          double v_y = coords[3*iVertex+1];
          double v_z = coords[3*iVertex+2];
          double v_f = vertexField[stride*iVertex+k];

          a[0][0] += v_x * v_x;
          a[0][1] += v_x * v_y;
          a[0][2] += v_x * v_z;
          a[0][3] += v_x;

          a[1][1] += v_y * v_y;
          a[1][2] += v_y * v_z;
          a[1][3] += v_y;

          a[2][2] += v_z * v_z;
          a[2][3] += v_z;

          a[3][3] += 1.;

          b[0] += v_x * v_f;
          b[1] += v_y * v_f;
          b[2] += v_z * v_f;
          b[3] += v_f;
        }
        a[1][0] = a[0][1];
        a[2][0] = a[0][2];
        a[3][0] = a[0][3];

        a[2][1] = a[1][2];
        a[3][1] = a[1][3];

        a[3][2] = a[2][3];

        if (solve_ax_b_4(a, b, coeff) == 0) {
          interpolatedField[stride*ipoint+k] = (coeff[0] * distantCoords[3*ipoint]
                                                                         + coeff[1] * distantCoords[3*ipoint+1]
                                                                                                    + coeff[2] * distantCoords[3*ipoint+2]
                                                                                                                               + coeff[3]);
        }
        else {
          interpolatedField[stride*ipoint+k] = vertexField[stride*nVertex+k]; /* last encountered value */
        }
      }
    }
    else {

      const int *polyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex();
      const int *polyhedraCellToFaceConnectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();
      const int *polyhedraFaceConnectivityIndex = _supportMesh->getPolyhedraFaceConnectivityIndex() ;
      const int *polyhedraFaceConnectivity = _supportMesh->getPolyhedraFaceConnectivity();

      int ipoly = iel - nStandardElt;
      int nFacePolyhedra = polyhedraFaceIndex[iel+1] - polyhedraFaceIndex[iel];
      int faceIndex = polyhedraCellToFaceConnectivity[iel];
      int nVertexFace = 0;

      for (int j = 0; j < nFacePolyhedra; j++) {
        int iface = polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
        nVertexFace += polyhedraFaceConnectivityIndex[iface+1] - polyhedraFaceConnectivityIndex[iface];
      }

      std::vector<int> vertexPoly(nVertexFace);
      for (int j = 0; j < nFacePolyhedra; j++) {
        int iface = polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
        int nVertexLocFace = polyhedraFaceConnectivityIndex[iface+1] - polyhedraFaceConnectivityIndex[iface];
        int vertexIndex = polyhedraFaceConnectivityIndex[iface];
        for (int k = 0; k < nVertexLocFace; k++)
          vertexPoly.push_back(polyhedraFaceConnectivity[vertexIndex+k]);
      }
      quickSort(&vertexPoly[0], 0, vertexPoly.size()-1, NULL);

      int iVertex = -1;
      double v_x;
      double v_y;
      double v_z;
      double v_f;

      for(int k = 0; k < stride; k++) {
        double a[4][4] = {{0., 0., 0., 0.},
            {0., 0., 0., 0.},
            {0., 0., 0., 0.},
            {0., 0., 0., 0.}};
        double b[4] = {0., 0., 0., 0.};

        for (int j = 0; j < vertexPoly.size(); j++) {
          if (iVertex < vertexPoly[j]) {
            iVertex = vertexPoly[j];
            v_x = coords[3*iVertex];
            v_y = coords[3*iVertex+1];
            v_z = coords[3*iVertex+2];
            v_f = vertexField[stride*iVertex+k];

            a[0][0] += v_x * v_x;
            a[0][1] += v_x * v_y;
            a[0][2] += v_x * v_z;
            a[0][3] += v_x;

            a[1][1] += v_y * v_y;
            a[1][2] += v_y * v_z;
            a[1][3] += v_y;

            a[2][2] += v_z * v_z;
            a[2][3] += v_z;

            a[3][3] += 1.;

            b[0] += v_x * v_f;
            b[1] += v_y * v_f;
            b[2] += v_z * v_f;
            b[3] += v_f;
          }
          a[1][0] = a[0][1];
          a[2][0] = a[0][2];
          a[3][0] = a[0][3];

          a[2][1] = a[1][2];
          a[3][1] = a[1][3];

          a[3][2] = a[2][3];

          if (solve_ax_b_4(a, b, coeff) == 0) {
            interpolatedField[stride*ipoint+k] = (coeff[0] * distantCoords[3*ipoint]
                                                                           + coeff[1] * distantCoords[3*ipoint+1]
                                                                                                      + coeff[2] * distantCoords[3*ipoint+2]
                                                                                                                                 + coeff[3]);
          }
          else {
            interpolatedField[stride*ipoint+k] = v_f; /* last encountered value */
          }
        }
      }
    }
  }
}


void Coupling::defineMesh(const int nVertex,
                          const int nElement,
                          const double coordinates[],
                          int connectivity_index[],
                          int connectivity[])
{
  if (_supportMesh  != NULL)
    bft_error(__FILE__, __LINE__, 0, "coupling mesh is already created\n");

  if (_isCoupledRank)
    _supportMesh = new Mesh(_fvmComm,
                            _entitiesDim,
                            nVertex,
                            nElement,
                            coordinates,
                            connectivity_index,
                            connectivity);
  else
    bft_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning," 
              " the coupling mesh must be defined only by the root rank\n");

}

void Coupling::setPointsToLocate(const int    n_points,
                                 double coordinate[])
{
  if (_isCoupledRank) {
    _nPointsToLocate = n_points;
    _coordsPointsToLocate = coordinate;
    _toLocate = true;
  }
  else
    bft_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning," 
              " the points to locate must be defined only by the root rank\n");
}


void Coupling::defineMeshAddPolyhedra(const int n_element,
                                      int face_index[],
                                      int cell_to_face_connectivity[],
                                      int face_connectivity_index[],
                                      int face_connectivity[])

{
  if (_isCoupledRank) {
    if (_supportMesh  == NULL)
      bft_error(__FILE__, __LINE__, 0, "No mesh to add elements\n");

    _supportMesh->addPolyhedra(n_element,
                               face_index,
                               cell_to_face_connectivity,
                               face_connectivity_index,
                               face_connectivity);
  }
  bft_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning," 
            " the coupling mesh must be defined only by the root rank\n");

}


void Coupling::updateLocation()
{
  if (_isCoupledRank)
    _toLocate = true;
  else
    bft_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning," 
              " updateLocation must be called only by the root rank\n");
}

cwipi_exchange_status_t Coupling::exchange(const char    *exchangeName,
                                               const int     stride,
                                               const int     timeStep,
                                               const double  timeValue,
                                               const char    *sendingFieldName,
                                               const double  *sendingField,
                                               char          *receivingFieldName,
                                               double        *receivingField,
                                               void          *ptFortranInterpolationFct)


{
  cwipi_exchange_status_t status = CWIPI_EXCHANGE_OK;
  
  const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();
  
  int rootRank;

  if (!_isCoupledRank && sendingField != NULL )
    bft_printf("Warning : sendingField != NULL, " 
              " only field defined by the root rank is sent\n");

  if (_isCoupledRank) {
    
    int currentRank;
    MPI_Comm_rank(_couplingComm, &currentRank);
    
    int lLocalName = _name.size() + 1;
    int lDistantName = 0;
    MPI_Status MPIStatus;

    if (currentRank == 0 || 
        currentRank == _coupledApplicationBeginningRankCouplingComm +
                       _coupledApplicationNRankCouplingComm) {

      //
      // Check coupling name

      MPI_Sendrecv(&lLocalName,   1, MPI_INT, 
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   &lDistantName, 1, MPI_INT, 
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   _couplingComm, &MPIStatus);

      char *distantCouplingName = new char[lDistantName];

      MPI_Sendrecv(const_cast <char*>(_name.c_str()), 
                   lLocalName, MPI_CHAR, 
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   distantCouplingName, lDistantName, MPI_CHAR,
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   _couplingComm, &MPIStatus);

      if (strcmp(_name.c_str(), distantCouplingName))
	bft_error(__FILE__, __LINE__, 0, "'%s' '%s' bad synchronization point\n",
		  _name.c_str(),
		  distantCouplingName);
      delete[] distantCouplingName;

      //
      // Check exchange name

      lLocalName = strlen(exchangeName)+1;
      MPI_Sendrecv(&lLocalName,   1, MPI_INT, 
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   &lDistantName, 1, MPI_INT, 
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   _couplingComm, &MPIStatus);

      char *distantExchangeName = new char[lDistantName];

      MPI_Sendrecv(const_cast <char*>(exchangeName),        
                   lLocalName, MPI_CHAR, 
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   distantExchangeName, lDistantName, MPI_CHAR, 
                   _coupledApplicationBeginningRankCouplingComm, 0,
		   _couplingComm, &MPIStatus);
      
      if (strcmp(exchangeName, distantExchangeName))
	bft_error(__FILE__, __LINE__, 0, "'%s' '%s' bad synchronization point\n",
		  exchangeName,
		  distantExchangeName);
      
      delete[] distantExchangeName;
    }
  }
    
  //
  // Locate
  
  locate();
  
  //
  // Prepare data (interpolate, extrapolate...)
  
  if (_isCoupledRank) {
    
    assert (_fvmLocator != NULL);
    
    const int nVertex                 = _supportMesh->getNVertex();
    const int nElts                   = _supportMesh->getNElts();
    const int nPoly                   = _supportMesh->getNPolyhedra();
    const int *localConnectivityIndex = _supportMesh->getEltConnectivityIndex();
    const int *localConnectivity      = _supportMesh->getEltConnectivity();
    const int *localPolyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex() ;
    const int *localPolyhedraCellToFaceConnectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();
    const int *localPolyhedraFaceConnectivity_index = _supportMesh->getPolyhedraFaceConnectivityIndex();
    const int *localPolyhedraFaceConnectivity       = _supportMesh->getPolyhedraFaceConnectivity();
    
    const int nDistantPoint     = fvm_locator_get_n_dist_points(_fvmLocator);
    const int *distantLocation  = fvm_locator_get_dist_locations(_fvmLocator);
    const double *distantCoords = fvm_locator_get_dist_coords(_fvmLocator);
    const int* interiorList     = fvm_locator_get_interior_list(_fvmLocator);
    const int nInteriorList     = fvm_locator_get_n_interior(_fvmLocator);
    
    int lDistantField = stride * nDistantPoint;
    int lReceivingField = stride * _nPointsToLocate;
    
    if (_tmpDistantField == NULL)
      _tmpDistantField = new std::vector<double> (lDistantField);
    
    std::vector<double>& tmpDistantField = *_tmpDistantField;
    if (tmpDistantField.size() < lDistantField)
      tmpDistantField.resize(lDistantField);
    
    //
    // Interpolation
    
    if (sendingField != NULL) {
      
      assert(!(_interpolationFct != NULL && ptFortranInterpolationFct != NULL));
      
      //
      // Callback Fortran
      
      if (ptFortranInterpolationFct != NULL)
        PROCF(callfortinterpfct, CALLFORTINTERPFCT) (
           const_cast <int *> (&_entitiesDim),
           const_cast <int *> (&nVertex),
           const_cast <int *> (&nElts),
           const_cast <int *> (&nPoly),
           const_cast <int *> (&nDistantPoint),
           const_cast <double *> (_supportMesh->getVertexCoords()),
           const_cast <int *> (localConnectivityIndex),
           const_cast <int *> (localConnectivity),
           const_cast <int *> (localPolyhedraFaceIndex),
           const_cast <int *> (localPolyhedraCellToFaceConnectivity),
           const_cast <int *> (localPolyhedraFaceConnectivity_index),
           const_cast <int *> (localPolyhedraFaceConnectivity),
           const_cast <double *> (distantCoords),
           const_cast <int *> (distantLocation),
           const_cast <int *> (_barycentricCoordinatesIndex),
           const_cast <double *> (_barycentricCoordinates),
           const_cast <int *> (&stride),
           const_cast <int *> ((const int *) &_solverType),
           const_cast <double *> (sendingField),
           const_cast <double *> (&tmpDistantField[0]),
           ptFortranInterpolationFct
          );

      //
      // Callback C

      else if (_interpolationFct != NULL)
        _interpolationFct(_entitiesDim,
                          nVertex,
                          nElts,
                          nPoly,
                          nDistantPoint,
                          _supportMesh->getVertexCoords(),
                          localConnectivityIndex,
                          localConnectivity,
                          localPolyhedraFaceIndex,
                          localPolyhedraCellToFaceConnectivity,
                          localPolyhedraFaceConnectivity_index,
                          localPolyhedraFaceConnectivity,
                          distantCoords,
                          distantLocation,
                          _barycentricCoordinatesIndex,
                          _barycentricCoordinates,
                          stride,
                          _solverType,
                          sendingField,
                          &tmpDistantField[0]);
      else
        _interpolate((double* )sendingField,
                     tmpDistantField,
                     stride);
    }

    //
    // Exchange

    double *ptSending = NULL;

    if (sendingField != NULL)
      ptSending = &tmpDistantField[0];

#ifdef NAN
    if (receivingField != NULL && nInteriorList > 0){
      const int idx = 0;
      receivingField[idx] = NAN;
    }
#endif

    fvm_parall_set_mpi_comm(_fvmComm);

    fvm_locator_exchange_point_var(_fvmLocator,
                                   (void *) ptSending,
                                   (void *) receivingField,
                                   NULL,
                                   sizeof(double),
                                   stride,
                                   0);

    fvm_parall_set_mpi_comm(MPI_COMM_NULL);

    //
    // Check receiving

#ifdef NAN
    if (receivingField != NULL && nInteriorList > 0) {
      const int idx = 0;
      if (isnan(receivingField[idx]))
        status = CWIPI_EXCHANGE_BAD_RECEIVING;
    }
#endif

    //
    // Visualization

    if (status == CWIPI_EXCHANGE_OK)

      _fieldsVisualization(exchangeName,
                           stride,
                           timeStep,
                           timeValue,
                           sendingFieldName,
                           sendingField,
                           receivingFieldName,
                           receivingField);
 
    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      
      MPI_Comm_rank(localComm, &rootRank);
      
      assert(rootRank == 0);
      
      MPI_Bcast( &status, 
                 1,
                 MPI_INT, 
                 0, 
                 localComm );

      if( receivingField != NULL ) 
        
        MPI_Bcast( receivingField, 
                   stride* (_nPointsToLocate-_nNotLocatedPoint),
                   MPI_DOUBLE, 
                   0, 
                   localComm );
      
    }
  }

  else if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {
    
    MPI_Bcast( &status, 
               1,
               MPI_INT, 
               0, 
               localComm );
    
    if ( receivingField != NULL ) 
      
      MPI_Bcast( receivingField, 
                 stride* (_nPointsToLocate-_nNotLocatedPoint),
                 MPI_DOUBLE, 
                 0, 
                 localComm );
  }

  return status;
}

void Coupling::_initVisualization()
{
  if (_fvmWriter == NULL && _outputFrequency > 0) {

    bft_file_mkdir_default("cwipi");

    std::string pathString = "cwipi/"+_name + "_" +
    _localApplicationProperties.getName()+"_"+
    _coupledApplicationProperties.getName();

    fvm_parall_set_mpi_comm(_fvmComm);

    _fvmWriter = fvm_writer_init("Chr",
                                 pathString.c_str(),
                                 _outputFormat.c_str(),
                                 _outputFormatOption.c_str(),
                                 FVM_WRITER_FIXED_MESH);

    int localCommSize;
    MPI_Comm_size(_fvmComm, &localCommSize);
    int localRank = 0;
    MPI_Comm_rank(_fvmComm, &localRank);

    fvm_writer_export_nodal(_fvmWriter, &(_supportMesh->getFvmNodal()));

    // Export sub domain

    if (localCommSize > 1) {

      const int nElts  = _supportMesh->getNElts();

      int *domLoc = new int [nElts];

      for (int i = 0; i < nElts; i++)
        domLoc[i] = localRank+1;

      fvm_writer_export_field(_fvmWriter,
                              const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                              "partitioning",
                              FVM_WRITER_PER_ELEMENT,
                              1,
                              FVM_NO_INTERLACE,
                              0,
                              NULL,
                              FVM_INT32,
                              1,
                              1.,
                              (const void *const *)  &domLoc);
      delete[] domLoc;
    }

    // TODO: A deplacer et a recreer en cas de maillage mobile

    int nNotLocatedPointSum;
    MPI_Allreduce (&_nNotLocatedPoint, &nNotLocatedPointSum, 
                   1, MPI_INT, MPI_SUM,
                   _fvmComm);

    if (nNotLocatedPointSum != 0 && _coordsPointsToLocate == NULL) {
      if (_solverType == CWIPI_SOLVER_CELL_CENTER) {
        const int nElts  = _supportMesh->getNElts();

        int *domLoc = new int [nElts];

        for (int i = 0; i < nElts; i++)
          domLoc[i] = 1;

        for (int i = 0; i < _nNotLocatedPoint; i++)
          domLoc[_notLocatedPoint[i]-1] = 0;

        fvm_writer_export_field(_fvmWriter,
                                const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                                "location",
                                FVM_WRITER_PER_ELEMENT,
                                1,
                                FVM_NO_INTERLACE,
                                0,
                                NULL,
                                FVM_INT32,
                                1,
                                1.,
                                (const void *const *)  &domLoc);
      }
      else {
        const int nVertex = _supportMesh->getNVertex();

        int *domLoc = new int [nVertex];

        for (int i = 0; i < nVertex; i++)
          domLoc[i] = 0;

        for (int i = 0; i < _nNotLocatedPoint; i++)
          domLoc[_notLocatedPoint[i]-1] = 1;

        fvm_writer_export_field(_fvmWriter,
                                const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                                "location",
                                FVM_WRITER_PER_NODE,
                                1,
                                FVM_NO_INTERLACE,
                                0,
                                NULL,
                                FVM_INT32,
                                1,
                                1.,
                                (const void *const *)  &domLoc);
        delete [] domLoc;
      }
    }

    fvm_parall_set_mpi_comm(MPI_COMM_NULL);

  }
}


void Coupling::_fieldsVisualization(const char *exchangeName,
                              const int stride,
                              const int timeStep,
                              const double timeValue,
                              const char  *sendingFieldName,
                              const void *sendingField,
                              const char  *receivingFieldName,
                              const void *receivingField)
{

  const double couplingDoubleMax = 1e33;

  std::string localName;

  if ((_outputFrequency > 0) && (timeStep % _outputFrequency == 0)) {

    assert(_fvmWriter != NULL);

    fvm_parall_set_mpi_comm(_fvmComm);

    fvm_writer_var_loc_t fvm_writer_var_loc;
    fvm_interlace_t fvm_interlace;
    int dim;

    if (_solverType == CWIPI_SOLVER_CELL_CENTER)
      fvm_writer_var_loc = FVM_WRITER_PER_ELEMENT;
    else
      fvm_writer_var_loc = FVM_WRITER_PER_NODE;

    //
    // Visualization if stride is scalar or vector

    if (stride == 1 || stride == 3) {

      dim = stride;
      if (stride == 1)
        fvm_interlace = FVM_NO_INTERLACE;
      else
        fvm_interlace = FVM_INTERLACE;

      if (sendingFieldName != NULL) {

        localName = "S_" + std::string(exchangeName) +
        "_" + std::string(sendingFieldName);

        if (sendingField != NULL) {

          std::vector<double> *cpSendingField = NULL;

          fvm_writer_export_field(_fvmWriter,
                                  const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                                  localName.c_str(),
                                  fvm_writer_var_loc,
                                  dim,
                                  fvm_interlace,
                                  0,
                                  NULL,
                                  FVM_DOUBLE,
                                  timeStep,
                                  timeValue,
                                  (const void *const *) &sendingField);
        }
      }

      if (receivingFieldName != NULL) {

        //
        // Not located point treatment

        if (receivingField != NULL && _coordsPointsToLocate == NULL) {

          localName = "R_" + std::string(exchangeName) +
          "_" + std::string(receivingFieldName);

          if (_nNotLocatedPoint != 0) {
            int lReceivingField = stride * _nPointsToLocate;
            std::vector<double> *cpReceivingField = new std::vector<double> (lReceivingField, couplingDoubleMax);

            const int nLocatedPoint = _nPointsToLocate - _nNotLocatedPoint;
            for (int i = 0; i < nLocatedPoint; i++) {
              for (int j = 0; j < stride; j++)
                (*cpReceivingField)[stride*(_locatedPoint[i]-1)+j] = ((double*) receivingField)[stride*i+j];
            }

            fvm_writer_export_field(_fvmWriter,
                                    const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                                    localName.c_str(),
                                    fvm_writer_var_loc,
                                    dim,
                                    fvm_interlace,
                                    0,
                                    NULL,
                                    FVM_DOUBLE,
                                    timeStep,
                                    timeValue,
                                    (const void *const *) cpReceivingField);
            delete cpReceivingField;
          }

          else
            fvm_writer_export_field(_fvmWriter,
                                    const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                                    localName.c_str(),
                                    fvm_writer_var_loc,
                                    dim,
                                    fvm_interlace,
                                    0,
                                    NULL,
                                    FVM_DOUBLE,
                                    timeStep,
                                    timeValue,
                                    (const void *const *) &receivingField);
        }
      }
    }
    fvm_parall_set_mpi_comm(MPI_COMM_NULL);
  }
  
}


void Coupling::locate()
{

  const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();
  
  
  if ( _toLocate ) {
    
    int rootRank = -1;

    if ( _isCoupledRank ) {
      
      fvm_parall_set_mpi_comm(_fvmComm);

      //
      // Create locator
      
      if (_fvmLocator == NULL)
        _fvmLocator = fvm_locator_create(_tolerance, 
                                         _couplingComm, 
                                         _coupledApplicationNRankCouplingComm, 
                                         _coupledApplicationBeginningRankCouplingComm);

      // TODO: Revoir les coordonnees des points a localiser (cas centres sommets + centres faces + autres points)
      // TODO: Ajouter un locator pour les sommets pour les centres faces,...
      
      double* coords = NULL;
      if (_coordsPointsToLocate != NULL)
        coords = _coordsPointsToLocate;
      
      else if(_solverType == CWIPI_SOLVER_CELL_CENTER) {
        _nPointsToLocate = _supportMesh->getNElts();
        coords = const_cast <double*> (&(_supportMesh->getCellCenterCoords()[0]));
      }
      
      else if(_solverType == CWIPI_SOLVER_CELL_VERTEX) {
        _nPointsToLocate = _supportMesh->getNVertex();
        coords = const_cast <double*> (_supportMesh->getVertexCoords());
      }
      
      fvm_locator_set_nodal(_fvmLocator,
                            &_supportMesh->getFvmNodal(),
                            0,
                            3,
                            _nPointsToLocate,
                            NULL,
                            coords);

      _toLocate = false;
      const int nLocatedPoint = fvm_locator_get_n_interior(_fvmLocator);
      const int nNotLocatedPoint = _nPointsToLocate - nLocatedPoint;
      const int* exteriorList = fvm_locator_get_exterior_list(_fvmLocator);
      const int* interiorList = fvm_locator_get_interior_list(_fvmLocator);
      const int* locationList = fvm_locator_get_dist_locations(_fvmLocator);
      const int nExterior = fvm_locator_get_n_exterior(_fvmLocator);
      assert(nNotLocatedPoint == nExterior);
      
      _notLocatedPoint = const_cast<int *> (exteriorList);
      _locatedPoint = const_cast<int *> (interiorList);
      _nNotLocatedPoint = nNotLocatedPoint;
      _location = const_cast<int *> (locationList);
      _nDistantpoint = fvm_locator_get_n_dist_points(_fvmLocator);
      
      if (_barycentricCoordinatesIndex != NULL) {
        if (_entitiesDim != 2) {
          delete[] _barycentricCoordinatesIndex;
          delete[] _barycentricCoordinates;
        }
        
        else {
          BFT_FREE(_barycentricCoordinatesIndex);
          BFT_FREE(_barycentricCoordinates);
        }
      }
      
      //
      // TODO: Prevoir une fabrique pour supprimer les tests if sur _entitiesDim
      //       Le calcul des coordonnees barycentriques se fera dans cette fabrique
      
      if (_barycentricCoordinatesIndex == NULL) {
        if (_entitiesDim == 1) {
          const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
          const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
          const double *distantCoords   = fvm_locator_get_dist_coords(_fvmLocator);
          
          const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
          const double *localCoords    = _supportMesh->getVertexCoords();
          
          _barycentricCoordinatesIndex = new int[nDistantPoint+1];
          _barycentricCoordinates = new double[2*nDistantPoint];
          _barycentricCoordinatesIndex[0] = 0;
          
          for (int ipoint = 0; ipoint < nDistantPoint ; ipoint++) {
            int iel = distantLocation[ipoint] - 1;
            _barycentricCoordinatesIndex[ipoint+1] = _barycentricCoordinatesIndex[ipoint] + 2;
            int index = eltsConnecPointer[iel];
            int nVertex = eltsConnecPointer[iel+1] - eltsConnecPointer[iel];
            assert(nVertex == 2);
            int pt1 = eltsConnecPointer[iel] - 1;
            int pt2 = eltsConnecPointer[iel+1] - 1;
            double coef1 = sqrt((localCoords[3*pt1]-distantCoords[3*ipoint])*(localCoords[3*pt1]-distantCoords[3*ipoint])+
                                (localCoords[3*pt1+1]-distantCoords[3*ipoint+1])*(localCoords[3*pt1+1]-distantCoords[3*ipoint+1])+
                                (localCoords[3*pt1+2]-distantCoords[3*ipoint+2])*(localCoords[3*pt1+2]-distantCoords[3*ipoint+2]));
            double coef2 = sqrt((localCoords[3*pt2]-distantCoords[3*ipoint])*(localCoords[3*pt2]-distantCoords[3*ipoint])+
                                (localCoords[3*pt2+1]-distantCoords[3*ipoint+1])*(localCoords[3*pt2+1]-distantCoords[3*ipoint+1])+
                                (localCoords[3*pt2+2]-distantCoords[3*ipoint+2])*(localCoords[3*pt2+2]-distantCoords[3*ipoint+2]));
            _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]] = coef1/(coef1+coef2);
            _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]+1] = coef2/(coef1+coef2);
          }
        }
        
        else if (_entitiesDim == 2) {
          int nPoints;
          coo_baryc(_fvmLocator,
                    _supportMesh->getNVertex(),
                    _supportMesh->getVertexCoords(),
                    _supportMesh->getNElts(),
                    _supportMesh->getEltConnectivityIndex(),
                    _supportMesh->getEltConnectivity(),
                    &nPoints,
                    &_barycentricCoordinatesIndex,
                    &_barycentricCoordinates);
        }
        
        else if (_entitiesDim == 3) {
          //TODO: calcul des coord barycentriques 3D
          int nPoints;
          //         coo_baryc(_fvmLocator,
          //                   _supportMesh->getNVertex(),
          //                   _supportMesh->getVertexCoords(),
          //                   _supportMesh->getNElts(),
          //                   _supportMesh->getEltConnectivityIndex(),
          //                   _supportMesh->getEltConnectivity(),
          //                   &nPoints,
          //                   &_barycentricCoordinatesIndex,
          //                   &_barycentricCoordinates);
        }
        
      }
      
      _initVisualization();
      
      //
      // Send to local proc  
      
      if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {
        
        MPI_Comm_rank(localComm, &rootRank);
        
        assert(rootRank == 0);
        
        MPI_Bcast( &_nPointsToLocate, 1, MPI_INT, rootRank, localComm );
        MPI_Bcast( &_nNotLocatedPoint, 1, MPI_INT, rootRank, localComm );
        
        MPI_Bcast( _locatedPoint,
                   _nPointsToLocate-_nNotLocatedPoint,
                   MPI_INT, 
                   rootRank, 
                   localComm );

        MPI_Bcast( _notLocatedPoint, _nNotLocatedPoint, MPI_INT, rootRank, localComm );
        
      }

      fvm_parall_set_mpi_comm(MPI_COMM_NULL);

    }

    //
    // Receive location  
    
    else if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {
      
      rootRank = 0;
      _toLocate = false;
     
      MPI_Bcast( &_nPointsToLocate, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast( &_nNotLocatedPoint, 1, MPI_INT, rootRank, localComm );

      if (_locatedPoint != NULL) 
        delete []  _locatedPoint;
      
      if (_notLocatedPoint != NULL) 
        delete []  _notLocatedPoint ;
      
      _locatedPoint = new int[_nPointsToLocate-_nNotLocatedPoint];
      _notLocatedPoint = new int[_nNotLocatedPoint];
      
      MPI_Bcast( _locatedPoint,
                 _nPointsToLocate-_nNotLocatedPoint,
                 MPI_INT, 
                 rootRank, 
                 localComm );
      MPI_Bcast( _notLocatedPoint, _nNotLocatedPoint, MPI_INT, rootRank, localComm );
    
    }

  }
}

void Coupling::_createCouplingComm()
{

  const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

  if( _couplingComm == MPI_COMM_NULL ) {

    const int localBeginningRank = _localApplicationProperties.getBeginningRank();
    const int localEndRank       = _localApplicationProperties.getEndRank();

    const MPI_Comm& globalComm = _localApplicationProperties.getGlobalComm();
    int currentRank;
    MPI_Comm_rank(globalComm, &currentRank);

    const int distantBeginningRank = _coupledApplicationProperties.getBeginningRank();
    const int distantEndRank = _coupledApplicationProperties.getEndRank();

    //
    // Define coupledApplicationComm

    int nLocalRank = localEndRank - localBeginningRank + 1;
    int nDistantRank = distantEndRank - distantBeginningRank + 1;
    
    assert(localBeginningRank != distantBeginningRank);

    //
    // TODO: Trouver un tag unique pour une paire de codes (Risque de pb dans MPI_Intercomm_create)
    //       Lors de la creation de 2 objets couplages de maniere simultanee
 
    const int tag = 1;

    MPI_Comm tmpInterComm;

    MPI_Intercomm_create(localComm, 0, globalComm, distantBeginningRank, tag, &tmpInterComm);

    int order;
    if (localBeginningRank < distantBeginningRank) {
      order = 0;
     _coupledApplicationBeginningRankCouplingComm = nLocalRank ;
    }
    else {
      order = 1;
     _coupledApplicationBeginningRankCouplingComm = 0 ;
    }

    MPI_Intercomm_merge(tmpInterComm, order, &_mergeComm);

    MPI_Comm_free(&tmpInterComm);

    int coupledApplicationCommSize;

    MPI_Comm_size(_mergeComm, &coupledApplicationCommSize);

    assert(coupledApplicationCommSize == nLocalRank + nDistantRank);

    //
    // Exchange coupling type

    cwipi_coupling_type_t* couplingTypes = 
      new cwipi_coupling_type_t[coupledApplicationCommSize];

    MPI_Allgather((void*)& _couplingType,
                  1, 
                  MPI_INT,
                  couplingTypes,
                  1,
                  MPI_INT,
                  _mergeComm);

    //
    // Check coupling type

    int begRank = 0;
    int endRank = 0;
    
    if (localBeginningRank < distantBeginningRank) {
      begRank = 0;
      endRank = nLocalRank;
    }
    else {
      begRank = nDistantRank;
      endRank = nLocalRank + nDistantRank;
    }

    for (int i = begRank; i < endRank; i ++) 
      if (couplingTypes[i] != _couplingType)
        bft_error(__FILE__, __LINE__, 0, 
                  "Two different coupling types for the '%s' application\n",
                  _localApplicationProperties.getName().c_str());


    _coupledApplicationNRankCouplingComm = nDistantRank;

    _isCoupledRank = _localApplicationProperties.getBeginningRank() == currentRank ||
       _couplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING;
      
    cwipi_coupling_type_t distantCouplingType;

    if (localBeginningRank < distantBeginningRank) 
      distantCouplingType = couplingTypes[nLocalRank];
    else
      distantCouplingType = couplingTypes[0];

    delete [] couplingTypes;

    //
    // Build coupling communicator

    if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING || 
        distantCouplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {

      int *rankList = new int[coupledApplicationCommSize];
      int nRankList;

      if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING && 
          distantCouplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {

        nRankList = 2;
        rankList[0] = 0;

        _coupledApplicationNRankCouplingComm = 1; 
        if (localBeginningRank < distantBeginningRank) { 
          rankList[1] = nLocalRank;
          _coupledApplicationBeginningRankCouplingComm = 1;
        }
        else { 
          rankList[1] = nDistantRank;
          _coupledApplicationBeginningRankCouplingComm = 0;
        }
      }
      
      else if (distantCouplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
        nRankList = 1 + nDistantRank;
        _coupledApplicationNRankCouplingComm = nDistantRank; 
        if (localBeginningRank < distantBeginningRank) {
          rankList[0] = 0;
          _coupledApplicationBeginningRankCouplingComm = 1;
          for (int i = 0; i < nDistantRank; i++) 
            rankList[i+1] = nLocalRank + i;
        }
        else {
          _coupledApplicationBeginningRankCouplingComm = 0;
          for (int i = 0; i < nDistantRank; i++) 
            rankList[i] = i;
          rankList[nDistantRank] = nDistantRank;
        }
      }

      else if (_couplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
        nRankList = 1 + nLocalRank;
        _coupledApplicationNRankCouplingComm = 1; 
        if (localBeginningRank < distantBeginningRank) {
          _coupledApplicationBeginningRankCouplingComm = nLocalRank;
          for (int i = 0; i < nLocalRank; i++) 
            rankList[i] = i;
          rankList[nLocalRank] = nLocalRank;
        }
       
        else {
          rankList[0] = 0;
          _coupledApplicationBeginningRankCouplingComm = 0;
          for (int i = 0; i < nLocalRank; i++) 
            rankList[i+1] = nDistantRank + i;
        }
      }

      else {
        
        bft_error(__FILE__, __LINE__, 0, 
                  "Error in 'build coupling communicator'\n");
      }

      MPI_Group mergeGroup = MPI_GROUP_NULL;
      MPI_Group couplingGroup = MPI_GROUP_NULL;

      MPI_Comm_group(_mergeComm, &mergeGroup);

      MPI_Group_incl(mergeGroup, nRankList, rankList, &couplingGroup);

      MPI_Comm_create(_mergeComm, couplingGroup, &_couplingComm);

      delete [] rankList;

    }
    else
      MPI_Comm_dup(_mergeComm, &_couplingComm);
  }


  if (_fvmComm == MPI_COMM_NULL) {
    
    if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
      
      int list1 = 0;
      MPI_Group localGroup = MPI_GROUP_NULL;
      MPI_Group fvmGroup = MPI_GROUP_NULL;
      
      MPI_Comm dupLocalComm = MPI_COMM_NULL;
      MPI_Comm_dup(localComm, &dupLocalComm);

      MPI_Comm_group(dupLocalComm, &localGroup);
      MPI_Group_incl(localGroup, 1, &list1, &fvmGroup);
      MPI_Comm_create(localComm,
                      fvmGroup,
                      &_fvmComm);
      MPI_Comm_free(&dupLocalComm);
    }
    
    else 
      
      MPI_Comm_dup(localComm, &_fvmComm);
 
  }
}


} // namespace cwipi

