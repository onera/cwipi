
#include <cassert>
#include <cmath>
#include <sstream>

#include <bft_error.h>
#include <bft_file.h>
#include <bft_printf.h>

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "mesh.hxx"
#include "applicationProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "coo_baryc.h"
#include "couplings.h"

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
      int *stored,
      int *global_num,
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
namespace couplings {

Coupling::Coupling(const std::string& name,
                   const ApplicationProperties& localApplicationProperties,
                   const ApplicationProperties& coupledApplicationProperties,
                   const int entitiesDim,
                   const double tolerance,
                   const couplings_solver_type_t solverType,
                   const int    outputFrequency,
                   const char  *outputFormat,
                   const char  *outputFormatOption)
:_localApplicationProperties(localApplicationProperties),
_coupledApplicationProperties(coupledApplicationProperties),
_entitiesDim(entitiesDim),_tolerance(tolerance), _solverType(solverType),
_outputFormat(outputFormat), _outputFormatOption(outputFormatOption),
_fvmWriter(NULL), _outputFrequency(outputFrequency), _name(name),
_nDistantpoint(0)

{
  _tmpVertexField = NULL;
  _tmpDistantField = NULL;
  _supportMesh = NULL;
  _coordsPointsToLocate = NULL;
  _fvmLocator = NULL;
  _interpolationFct = NULL;
  _toLocate = true;
  _barycentricCoordinatesIndex = NULL;
  _barycentricCoordinates = NULL;
  _nNotLocatedPoint = 0;
  _nPointsToLocate = 0;
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
        if (_supportMesh->getParentNum() != NULL)
          vertexField[vertex]  += cellCenterField[(*_supportMesh->getParentNum())[i]-1] * cellVolume[i];
        else
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
              if (_supportMesh->getParentNum() != NULL)
                vertexField[ivertex - 1]  += cellCenterField[(*_supportMesh->getParentNum())[i]-1] * cellVolume[i];
              else
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
  // BUG pas bug : mis ca dans bft au lieu de std::cout
  bft_printf( "destroying '" );
  bft_printf( _name.c_str() );
  bft_printf( "' coupling\n" );
  delete _tmpVertexField;
  delete _tmpDistantField;
  delete _supportMesh;

  if (_entitiesDim != 2) {
    delete[] _barycentricCoordinatesIndex;
    delete[] _barycentricCoordinates;
  }

  else {
    BFT_FREE(_barycentricCoordinatesIndex);
    BFT_FREE(_barycentricCoordinates);
  }

  if (_fvmLocator != NULL)
    fvm_locator_destroy(_fvmLocator);
  if (_fvmWriter != NULL)
    fvm_writer_finalize(_fvmWriter);
}

void Coupling::_interpolate(double *referenceField,
                            std::vector<double>& interpolatedField,
                            const int stride)
{
  if (_solverType == COUPLINGS_SOLVER_CELL_CENTER)
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
  // TODO: Dans un premier temps faire le calcul pour les tétraèdres

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
    bft_error(__FILE__, __LINE__, 0, "mesh is already created\n");

  _supportMesh = new Mesh(_entitiesDim,
                          nVertex,
                          nElement,
                          coordinates,
                          connectivity_index,
                          connectivity);

}

void Coupling::setPointsToLocate(const int    n_points,
                                 double coordinate[])
{
  _nPointsToLocate = n_points;
  _coordsPointsToLocate = coordinate;
  _toLocate = true;
}


void Coupling::defineMeshAddPolyhedra(const int n_element,
                                      int face_index[],
                                      int cell_to_face_connectivity[],
                                      int face_connectivity_index[],
                                      int face_connectivity[])

{
  if (_supportMesh  == NULL)
    bft_error(__FILE__, __LINE__, 0, "No mesh to add elements\n");

  _supportMesh->addPolyhedra(n_element,
                             face_index,
                             cell_to_face_connectivity,
                             face_connectivity_index,
                             face_connectivity);
}


void Coupling::updateLocation()
{
  _toLocate = true;
}

couplings_exchange_status_t Coupling::exchange(const char                          *exchangeName,
                                               const int                            stride,
                                               //const couplings_field_dimension_t    fieldDimension,
                                               const int                            timeStep,
                                               const double                         timeValue,
                                               const char                          *sendingFieldName,
                                               const double                        *sendingField,
                                               char                                *receivingFieldName,
                                               double                              *receivingField,
                                               void                                *ptFortranInterpolationFct)


{

//  std::cout << "couplings_exchange" << std::endl;
  couplings_exchange_status_t status = COUPLINGS_EXCHANGE_OK;

  //
  // Check exchange_name

  const int localBeginningRank = _localApplicationProperties.getBeginningRank();
  const int distantBeginningRank = _coupledApplicationProperties.getBeginningRank();
  const MPI_Comm& globalComm = _localApplicationProperties.getGlobalComm();
  const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

  int currentRank;
  MPI_Comm_rank(globalComm, &currentRank);

  int lLocalName = _name.size() + 1;
  int lDistantName = 0;
  MPI_Status MPIStatus;

  if (currentRank == localBeginningRank) {

    //
    // Check coupling name

    MPI_Sendrecv(&lLocalName,   1, MPI_INT, distantBeginningRank, 0,
                 &lDistantName, 1, MPI_INT, distantBeginningRank, 0,
                 globalComm, &MPIStatus);

    char *distantCouplingName = new char[lDistantName];

    MPI_Sendrecv(const_cast <char*>(_name.c_str()),        lLocalName, MPI_CHAR, distantBeginningRank, 0,
                 distantCouplingName, lDistantName, MPI_CHAR, distantBeginningRank, 0,
                 globalComm, &MPIStatus);

    if (strcmp(_name.c_str(), distantCouplingName))
      bft_error(__FILE__, __LINE__, 0, "'%s' '%s' bad synchronisation point\n",
                _name.c_str(),
                distantCouplingName);
    delete[] distantCouplingName;

    //
    // Check exchange name

    lLocalName = strlen(exchangeName)+1;
    MPI_Sendrecv(&lLocalName,   1, MPI_INT, distantBeginningRank, 0,
                 &lDistantName, 1, MPI_INT, distantBeginningRank, 0,
                 globalComm, &MPIStatus);

    char *distantExchangeName = new char[lDistantName];

    MPI_Sendrecv(const_cast <char*>(exchangeName),        lLocalName, MPI_CHAR, distantBeginningRank, 0,
                 distantExchangeName, lDistantName, MPI_CHAR, distantBeginningRank, 0,
                 globalComm, &MPIStatus);

    if (strcmp(exchangeName, distantExchangeName))
      bft_error(__FILE__, __LINE__, 0, "'%s' '%s' bad synchronisation point\n",
                exchangeName,
                distantExchangeName);
    delete[] distantExchangeName;
  }

  //
  // Locate
//  std::cout << "avant locate" << std::endl;
  locate();
//  std::cout << "apres locate" << std::endl;
  //
  // Prepare data (interpolate, extrapolate...)

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

  int lDistantField = 0;
  int lReceivingField = 0;

  lDistantField = stride * nDistantPoint;
  lReceivingField = stride * _nPointsToLocate;

  if (_tmpDistantField == NULL)
    _tmpDistantField = new std::vector<double> (lDistantField);

    std::vector<double>& tmpDistantField = *_tmpDistantField;
    if (tmpDistantField.size() < lDistantField)
      tmpDistantField.resize(lDistantField);

    //
    // Interpolation

//    std::cout << "avant interpolation" << std::endl;
    if (sendingField != NULL) {

      assert(!(_interpolationFct != NULL && ptFortranInterpolationFct != NULL));

      // Callback Fortran

      if (ptFortranInterpolationFct != NULL) {

        if (_supportMesh->getParentNum() == NULL) {

          int stored = 1;
          PROCF(callfortinterpfct, CALLFORTINTERPFCT) (const_cast <int *> (&_entitiesDim),
              const_cast <int *> (&nVertex),
              const_cast <int *> (&nElts),
              const_cast <int *> (&nPoly),
              const_cast <int *> (&nDistantPoint),
              const_cast <double *> (_supportMesh->getVertexCoords()),
              const_cast <int *> (&stored),
              const_cast <int *> (&stored),
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
        }

        else {

          int stored = 0;
          PROCF(callfortinterpfct, CALLFORTINTERPFCT) (const_cast <int *> (&_entitiesDim),
              const_cast <int *> (&nVertex),
              const_cast <int *> (&nElts),
              const_cast <int *> (&nPoly),
              const_cast <int *> (&nDistantPoint),
              const_cast <double *> (_supportMesh->getVertexCoords()),
              const_cast <int *> (&stored),
              const_cast <int *> (&(*_supportMesh->getParentNum())[0]),
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
        }
      }

      // Callback C

      else if (_interpolationFct != NULL) {

        if (_supportMesh->getParentNum() == NULL)
          _interpolationFct(_entitiesDim,
                            nVertex,
                            nElts,
                            nPoly,
                            nDistantPoint,
                            _supportMesh->getVertexCoords(),
                            NULL,
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
          _interpolationFct(_entitiesDim,
                            nVertex,
                            nElts,
                            nPoly,
                            nDistantPoint,
                            _supportMesh->getVertexCoords(),
                            &(*_supportMesh->getParentNum())[0],
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
      }
      else
        _interpolate((double* )sendingField,
                     tmpDistantField,
                     stride);
    }

//    std::cout << "apres interpolation" << std::endl;
    //
    // Exchange

    double *ptSending = NULL;

    if (sendingField != NULL)
      ptSending = &tmpDistantField[0];

//    std::cout << "avant createNaN" << std::endl;
    if (receivingField != NULL && nInteriorList > 0){
      const int idx = 0;
      receivingField[idx] = _createNan();
    }

//    fvm_locator_dump( _fvmLocator );
//    std::cout << "avant fvm_locator_exchange_point_var" << std::endl;
    fvm_locator_exchange_point_var(_fvmLocator,
                                   (void *) ptSending,
                                   (void *) receivingField,
                                   NULL,
                                   sizeof(double),
                                   stride,
                                   0);

    //
    // Check receiving

//    std::cout << "apres echange" << std::endl;
    if (receivingField != NULL && nInteriorList > 0) {
      std::ostringstream os;
      const int idx = 0;
      os << receivingField[idx];
      if ((os.str()[0] == 'n') || (os.str()[0] == 'N'))
        status = COUPLINGS_EXCHANGE_BAD_RECEIVING;
    }

//    std::cout << "avant check" << std::endl;
    //
    // Not located point treatment
    // TODO: A supprimer si on ne définit les champs qu'aux points localises

    if (receivingField != NULL) {

      if (_nNotLocatedPoint != 0 && status == COUPLINGS_EXCHANGE_OK) {
        std::vector<double> cpReceivingField(lReceivingField);
        for (int i = 0; i < lReceivingField; i++)
          cpReceivingField[i] = receivingField[i];

        const int nLocatedPoint = _nPointsToLocate - _nNotLocatedPoint;
        for (int i = 0; i < nLocatedPoint; i++) {
          for (int j = 0; j < stride; j++)
            receivingField[stride*(_locatedPoint[i]-1)+j] = cpReceivingField[stride*i+j];
        }
        for (int i = 0; i < _nNotLocatedPoint; i++) {
          for (int j = 0; j < stride; j++)
            receivingField[stride*(_notLocatedPoint[i]-1)+j] = _createNan();
        }
      }
    }

    //
    // Visualization

//    std::cout << "avant visualization" << std::endl;
    _visualization(exchangeName,
                   stride,
                   timeStep,
                   timeValue,
                   sendingFieldName,
                   sendingField,
                   receivingFieldName,
                   receivingField);

//    std::cout << "fin fonction exchange" << std::endl;
    return status;
}


void Coupling::_visualization(const char *exchangeName,
                              const int stride,
                              const int timeStep,
                              const double timeValue,
                              const char  *sendingFieldName,
                              const void *sendingField,
                              const char  *receivingFieldName,
                              const void *receivingField)
{

  std::string localName;

  if ((_outputFrequency > 0) && (timeStep % _outputFrequency == 0)) {

    bft_file_mkdir_default("couplings");

    std::string pathString = "couplings/"+_name + "_" +
    _localApplicationProperties.getName()+"_"+
    _coupledApplicationProperties.getName();

    if (_fvmWriter == NULL) {
      _fvmWriter = fvm_writer_init("Chr",
                                   pathString.c_str(),
                                   _outputFormat.c_str(),
                                   _outputFormatOption.c_str(),
                                   FVM_WRITER_FIXED_MESH);

      int localCommSize = 0;
      const MPI_Comm &localComm = _localApplicationProperties.getLocalComm();
      MPI_Comm_size(localComm, &localCommSize);
      int localRank = 0;
      MPI_Comm_rank(localComm, &localRank);

      if (localCommSize > 1) {

        int *allNElts     = new int[localCommSize];

        //
        // global element num

        int nElts = _supportMesh->getNElts();
        unsigned int *globalEltNum = new unsigned int[nElts];

        MPI_Allgather(&nElts,
                      1,
                      MPI_INT,
                      allNElts,
                      1,
                      MPI_INT,
                      localComm);

        int nGlobal = 0;
        for(int i = 0; i < localRank; i++)
          nGlobal += allNElts[i];

        for(int i = 0; i < nElts; i++)
          globalEltNum[i] = nGlobal + i + 1;

        fvm_nodal_order_faces(&(_supportMesh->getFvmNodal()), globalEltNum);
        fvm_nodal_init_io_num(&(_supportMesh->getFvmNodal()), globalEltNum, 2);

        delete[] globalEltNum;

        //
        // global vertex num

        int nVertex = _supportMesh->getNVertex();
        unsigned int *globalVertexNum = new unsigned int[nVertex];

        MPI_Allgather(&nVertex,
                      1,
                      MPI_INT,
                      allNElts,
                      1,
                      MPI_INT,
                      localComm);

        nGlobal = 0;
        for(int i = 0; i < localRank; i++)
          nGlobal += allNElts[i];

        for(int i = 0; i < nVertex; i++)
          globalVertexNum[i] = nGlobal + i + 1;

        fvm_nodal_order_vertices(&(_supportMesh->getFvmNodal()), globalVertexNum);
        fvm_nodal_init_io_num(&(_supportMesh->getFvmNodal()), globalVertexNum, 0);

        delete[] globalVertexNum;
        delete[] allNElts;

      }

#if defined(DEBUG) && 0

      fvm_nodal_dump(&(_supportMesh->getFvmNodal()));

#endif
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

      if (_notLocatedPoint != NULL && _coordsPointsToLocate == NULL) {
        if (_solverType == COUPLINGS_SOLVER_CELL_CENTER) {
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
    }

    fvm_writer_var_loc_t fvm_writer_var_loc;
    fvm_interlace_t fvm_interlace;
    int dim;

    if (_solverType == COUPLINGS_SOLVER_CELL_CENTER)
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

          if (_supportMesh->getParentNum() != NULL && _solverType == COUPLINGS_SOLVER_CELL_CENTER) {
            int lSendingField = 0;
            lSendingField = stride * _nPointsToLocate;

            cpSendingField = new std::vector<double>(lSendingField);
            std::vector<double> &refCpSendingField = *cpSendingField;

            for (int i = 0; i < _nPointsToLocate; i++) {
              for (int j = 0; j < stride; j++)
                refCpSendingField[stride*i+j] = ((double*)sendingField)[stride*((*_supportMesh->getParentNum())[i]-1)+j];
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
                                    (const void *const *) &refCpSendingField);
            delete cpSendingField;
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
                                    (const void *const *) &sendingField);
        }
      }

      if (receivingFieldName != NULL) {
        if (receivingField != NULL && _coordsPointsToLocate == NULL) {
          std::vector<double> *cpReceivingField = NULL;

          localName = "R_" + std::string(exchangeName) +
          "_" + std::string(receivingFieldName);

          if (_supportMesh->getParentNum() != NULL && _solverType == COUPLINGS_SOLVER_CELL_CENTER) {
            int lReceivingField = 0;
            lReceivingField = stride * _nPointsToLocate;

            cpReceivingField = new std::vector<double>(lReceivingField);
            std::vector<double> &refCpReceivingField = *cpReceivingField;

            for (int i = 0; i < _nPointsToLocate; i++) {
              for (int j = 0; j < stride; j++)
                refCpReceivingField[stride*i+j] =
                  ((double*)receivingField)[stride*((*_supportMesh->getParentNum())[i]-1)+j];
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
                                    (const void *const *) &cpReceivingField[0]);
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
  }
}


void Coupling::locate()
{
  if (_fvmLocator == NULL || _toLocate) {

    if (_fvmLocator == NULL) {
      const MPI_Comm& globalComm = _localApplicationProperties.getGlobalComm();
      const int begRank = _coupledApplicationProperties.getBeginningRank();
      const int nRank = _coupledApplicationProperties.getEndRank() - begRank + 1;
      _fvmLocator = fvm_locator_create(_tolerance, globalComm, nRank, begRank);
    }

    double* coords = NULL;
    if (_coordsPointsToLocate != NULL)
      coords = _coordsPointsToLocate;

    else if(_solverType == COUPLINGS_SOLVER_CELL_CENTER) {
      _nPointsToLocate = _supportMesh->getNElts();
      if (_supportMesh->getParentNum() != NULL) {
        const std::vector<int> & parentNum = *_supportMesh->getParentNum();
        const std::vector<double> & couplingCellCoords = _supportMesh->getCellCenterCoords();
        coords = new double[3*_nPointsToLocate];
        for (int i = 0; i < _nPointsToLocate; i++) {
          coords[3*(parentNum[i]-1)]   = couplingCellCoords[3*i];
          coords[3*(parentNum[i]-1)+1] = couplingCellCoords[3*i+1];
          coords[3*(parentNum[i]-1)+2] = couplingCellCoords[3*i+2];
        }
      }
      else
        coords = const_cast <double*> (&(_supportMesh->getCellCenterCoords()[0]));
    }

    else if(_solverType == COUPLINGS_SOLVER_CELL_VERTEX) {
      _nPointsToLocate = _supportMesh->getNVertex();
      coords = const_cast <double*> (_supportMesh->getVertexCoords());
    }

    //fvm_nodal_dump(&(_supportMesh->getFvmNodal()));
    //std::cout << "toto" <<std::endl;
    //std::cout << "nPointsToLocate=" << _nPointsToLocate << std::endl;
    //std::cout << "coords(" << coords << ")="; for( int l=0; l<_nPointsToLocate; l++) std::cout << coords[3*l+0] << " " << coords[3*l+1] << " "<< coords[3*l+2] << std::endl;
    //for(int l=0; l<_nPointsToLocate; l++) std::cout << _supportMesh->getCellCenterCoords()[l] << " "; std::cout << std::endl;

    fvm_locator_set_nodal(_fvmLocator,
                          &_supportMesh->getFvmNodal(),
                          0,
                          3,
                          _nPointsToLocate,
                          NULL,
                          coords);
    //std::cout << "milieu locate" << std::endl;
    if(_solverType == COUPLINGS_SOLVER_CELL_CENTER && _supportMesh->getParentNum() != NULL)
      delete[] coords;

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

    if (_barycentricCoordinatesIndex == NULL) {
      if (_entitiesDim == 1) {
        const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
        const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
        const double *distantCoords   = fvm_locator_get_dist_coords(_fvmLocator);

        const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
        const double *localCoords         = _supportMesh->getVertexCoords();

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
    }
  }
}

double  Coupling::_createNan()
{
// Creation artificielle d'un Nan
//  double big = 3e99999999;
//  double userNan = big/big;
//  std::ostringstream os;
//  os << userNan;
//  if ((os.str()[0] != 'n') && (os.str()[0] != 'N'))
//    bft_error(__FILE__, __LINE__, 0, "'%f' %s bad nan detection\n", userNan, os.str().c_str());
//  return userNan;
  return 0;
}


} // namespace couplings

