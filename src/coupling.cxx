
#include <cassert>
#include <cmath>
#include <bft_error.h>

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "mesh.hxx"
#include "applicationProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "coo_baryc.h"
#include "couplings.h"

namespace couplings {

  Coupling::Coupling(const std::string& name,
                     const ApplicationProperties& localApplicationProperties,
                     const ApplicationProperties& coupledApplicationProperties,
                     const int entitiesDim,
                     const int tolerance,
                     const couplings_solver_type_t solverType,
                     const int    outputFrequency,
                     const char  *outputFormat,
                     const char  *outputFormatOption)
    :_localApplicationProperties(localApplicationProperties),
     _coupledApplicationProperties(coupledApplicationProperties),
     _entitiesDim(entitiesDim),_tolerance(tolerance), _solverType(solverType),
     _outputFormat(outputFormat), _outputFormatOption(outputFormatOption),
     _fvmWriter(NULL), _outputFrequency(outputFrequency)
    
  {
    _tmpVertexField = NULL;
    _supportMesh = NULL;
    _coordsPointsToLocate = NULL;
    _fvmLocator = NULL;
    _interpolationFct = NULL;
    _toLocate = true;
    _barycentricCoordinatesIndex = NULL;
    _barycentricCoordinates = NULL;
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
    delete _tmpVertexField;
    delete _supportMesh;
    delete _coordsPointsToLocate;
    fvm_locator_destroy(_fvmLocator);
  }

  void Coupling::_interpolate(double *referenceField, std::vector<double>& interpolatedField)
  {
    if (_solverType == COUPLINGS_SOLVER_CELL_CENTER)
      referenceField = &_extrapolate(referenceField)[0];

    switch(_entitiesDim) {

    case 1 :
      _interpolate1D(referenceField, interpolatedField);
      break;
      
    case 2 :
      _interpolate2D(referenceField, interpolatedField);
      break;
      
    case 3 :
      _interpolate3D(referenceField, interpolatedField);
      break;
      
    default:
      bft_error(__FILE__, __LINE__, 0, "'%i' bad entities dimension\n",_entitiesDim);
    }
  }

  void Coupling::_interpolate1D(double *referenceVertexField, std::vector<double>& interpolatedField)
  {
    const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
    const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
    const double *distantCoords   = fvm_locator_get_dist_coords(_fvmLocator);
    const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
    const int *eltsConnec        = _supportMesh->getEltConnectivity();
    const double *coords         = _supportMesh->getVertexCoords();

    if (_barycentricCoordinatesIndex == NULL) {
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
        double coef1 = sqrt((coords[3*pt1]-distantCoords[3*ipoint])*(coords[3*pt1]-distantCoords[3*ipoint])+
                            (coords[3*pt1+1]-distantCoords[3*ipoint+1])*(coords[3*pt1+1]-distantCoords[3*ipoint+1])+
                            (coords[3*pt1+2]-distantCoords[3*ipoint+2])*(coords[3*pt1+2]-distantCoords[3*ipoint+2]));
        double coef2 = sqrt((coords[3*pt2]-distantCoords[3*ipoint])*(coords[3*pt2]-distantCoords[3*ipoint])+
                            (coords[3*pt2+1]-distantCoords[3*ipoint+1])*(coords[3*pt2+1]-distantCoords[3*ipoint+1])+
                            (coords[3*pt2+2]-distantCoords[3*ipoint+2])*(coords[3*pt2+2]-distantCoords[3*ipoint+2]));
        _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]] = coef1/(coef1+coef2); 
        _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]+1] = coef2/(coef1+coef2);
      }
    }
    for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
      int iel = distantLocation[ipoint] - 1;
      double coef1 = _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]];
      double coef2 = _barycentricCoordinates[_barycentricCoordinatesIndex[ipoint]+1];
      int pt1 = eltsConnecPointer[iel] - 1;
      int pt2 = eltsConnecPointer[iel+1] - 1;

      interpolatedField[ipoint] = coef1 * referenceVertexField[pt1] + coef2 * referenceVertexField[pt2]; 
    }
  }
  
  void Coupling::_interpolate2D(double *vertexField, std::vector<double>& interpolatedField)
  {
    if (_barycentricCoordinatesIndex == NULL) {
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
    
    const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
    const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
    const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
    const int *eltsConnec        = _supportMesh->getEltConnectivity();

    for (int ipoint = 0; ipoint <nDistantPoint; ipoint++) {
      int iel = distantLocation[ipoint] - 1;
      interpolatedField[ipoint] = 0;
      int index = _barycentricCoordinatesIndex[ipoint];
      int nSom = _barycentricCoordinatesIndex[ipoint+1] - index;
      for (int isom = 0; isom <  nSom; isom++) {
        interpolatedField[ipoint] += vertexField[eltsConnec[eltsConnecPointer[iel]+isom]-1]
                                      *_barycentricCoordinates[index+isom];
      }
    }
  }

  void Coupling::_interpolate3D(double *vertexField, std::vector<double>& interpolatedField)
  {

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
          double v_f = vertexField[iVertex];

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
          interpolatedField[ipoint] = (coeff[0] *  distantCoords[3*ipoint]
                                       + coeff[1] * distantCoords[3*ipoint+1]
                                       + coeff[2] * distantCoords[3*ipoint+2]
                                       + coeff[3]);
        }
        else {
          interpolatedField[ipoint] = vertexField[nVertex]; /* last encountered value */
        }
      }
      else {
        const int *polyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex();
        const int *polyhedraCellToFaceConnectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();
        const int *polyhedraFaceConnectivityIndex = _supportMesh->getPolyhedraFaceConnectivityIndex() ;
        const int *polyhedraFaceConnectivity = _supportMesh->getPolyhedraFaceConnectivity();
        double a[4][4] = {{0., 0., 0., 0.},
                          {0., 0., 0., 0.},
                          {0., 0., 0., 0.},
                          {0., 0., 0., 0.}};
        double b[4] = {0., 0., 0., 0.};
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
        for (int j = 0; j < vertexPoly.size(); j++) {
          if (iVertex < vertexPoly[j]) {
            iVertex = vertexPoly[j];
            v_x = coords[3*iVertex];
            v_y = coords[3*iVertex+1];
            v_z = coords[3*iVertex+2];
            v_f = vertexField[iVertex];

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
            interpolatedField[ipoint] = (coeff[0] *  distantCoords[3*ipoint]
                                         + coeff[1] * distantCoords[3*ipoint+1]
                                         + coeff[2] * distantCoords[3*ipoint+2]
                                         + coeff[3]);
          }
          else {
            interpolatedField[ipoint] = v_f; /* last encountered value */
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

  int Coupling::exchange(const char                          *exchangeName,
                         const couplings_field_dimension_t    fieldDimension, 
                         const int                            timeStep, 
                         const double                         timeValue,
                         const char                          *sendingFieldName,
                         const double                        *sendingField, 
                         char                                *receivingFieldName,
                         double                              *receivingField)
  {
    
    int isReceived = 1;

    //
    // Check exchange_name

    const int localBeginningRank = _localApplicationProperties.getBeginningRank();
    const int distantBeginningRank = _coupledApplicationProperties.getBeginningRank();
    const MPI_Comm& globalComm = _localApplicationProperties.getGlobalComm();
    const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

    int currentRank;
    MPI_Comm_rank(globalComm, &currentRank);

    int lLocalName = strlen(exchangeName)+1;
    int lDistantName = 0;
    MPI_Status status;

    if (currentRank == localBeginningRank) {

      MPI_Sendrecv(&lLocalName,   1, MPI_INT, distantBeginningRank, 0,
                   &lDistantName, 1, MPI_INT, distantBeginningRank, 0,
                   globalComm, &status);
      
      char *distantExchangeName = new char[lDistantName];
    
      MPI_Sendrecv(const_cast <char*>(exchangeName),        lLocalName, MPI_CHAR, distantBeginningRank, 0,
                   distantExchangeName, lDistantName, MPI_CHAR, distantBeginningRank, 0,
                   globalComm, &status);

      if (strcmp(exchangeName, distantExchangeName))
        bft_error(__FILE__, __LINE__, 0, "'%s' '%s' bad synchronisation point\n", 
                  exchangeName, 
                  distantExchangeName);
    }

    //
    // Locate
    
    _locate();

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

    const int nDistantPoint      = fvm_locator_get_n_dist_points(_fvmLocator);
    const int *distantLocation   = fvm_locator_get_dist_locations(_fvmLocator);
    const double *distantCoords   = fvm_locator_get_dist_coords(_fvmLocator);

    int lDistantField = 0;
    int lReceivingField = 0;
    if (fieldDimension == COUPLINGS_FIELD_DIMENSION_INTERLACED_VECTOR) {
      lDistantField = 3*nDistantPoint;
      lReceivingField = 3*_nPointsToLocate;
    }
    else if (fieldDimension == COUPLINGS_FIELD_DIMENSION_SCALAR){
      lDistantField = nDistantPoint;
      lReceivingField = _nPointsToLocate;
    }

    if (_tmpDistantField == NULL) 
      _tmpDistantField = new std::vector<double> (lDistantField);
    
    std::vector<double>& tmpDistantField = *_tmpDistantField;
    if (tmpDistantField.size() < lDistantField)
      tmpDistantField.resize(lDistantField);

    //
    // Interpolation

    if (sendingField != NULL) {

      if (_interpolationFct != NULL) {

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
                            fieldDimension,
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
                            fieldDimension,
                            _solverType,
                            sendingField,
                            &tmpDistantField[0]);
      }
      else
        _interpolate((double* )sendingField, tmpDistantField);
    }

    //
    // Exchange

    double *ptSending = NULL;

    if (sendingField == NULL)
      ptSending = NULL;
    else
      ptSending = &tmpDistantField[0];

    if (receivingField != NULL) {
      double big = 3e99999999;
      receivingField[0] = big/big;

      if (receivingField[0] == receivingField[0])
        bft_error(__FILE__, __LINE__, 0, "'%f' bad entities dimension\n",receivingField[0]);
    }

    switch(fieldDimension) {

    case (COUPLINGS_FIELD_DIMENSION_SCALAR):
      fvm_locator_exchange_point_var(_fvmLocator,
                                     (void *) ptSending,
                                     (void *) receivingField,
                                     NULL,
                                     sizeof(double),
                                     1,
                                     0);
      break;

    case (COUPLINGS_FIELD_DIMENSION_INTERLACED_VECTOR):
      fvm_locator_exchange_point_var(_fvmLocator,
                                     (void *) ptSending,
                                     (void *) receivingField,
                                     NULL,
                                     sizeof(double),
                                     3,
                                     0);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "'%i' bad field dimension\n", 
                fieldDimension);
    }
    
    if (receivingField == NULL)
      isReceived = -1;
    else if ((_nPointsToLocate > 0) && 
             (_nPointsToLocate > _nNotLocatedPoint) &&
             (receivingField[0] != receivingField[0]))
      isReceived =  0;
    else if (_nPointsToLocate == _nNotLocatedPoint)
      isReceived = -2;
      
    int localCommSize = 0;
    MPI_Comm_size(localComm, &localCommSize);

    int *allIsReceived = new int[localCommSize];

    MPI_Allgather(&isReceived, 
                  1, 
                  MPI_INT, 
                  allIsReceived, 
                  1, 
                  MPI_INT, 
                  localComm);
    
    isReceived = 1;
    for (int i = 0 ; i < localCommSize ; i++) {
      if (allIsReceived[i] == -1) {
        isReceived = 0;
        break;
      }
    }
    
    //
    // Not located point treatment
 
    if (_nNotLocatedPoint != 0 && isReceived == 1) {
      std::vector<double> cpReceivingField(lReceivingField);
      for (int i = 0; i < lReceivingField; i++) 
        cpReceivingField[i] = receivingField[i];

      for (int i = 0; i < _nPointsToLocate; i++) {

        switch(fieldDimension) {
          
        case (COUPLINGS_FIELD_DIMENSION_SCALAR):
          receivingField[(*_supportMesh->getParentNum())[i]-1] =  cpReceivingField[i];          
          break;
          
        case (COUPLINGS_FIELD_DIMENSION_INTERLACED_VECTOR):
          receivingField[3*((*_supportMesh->getParentNum())[i]-1)]   = cpReceivingField[3*i];          
          receivingField[3*((*_supportMesh->getParentNum())[i]-1)+1] = cpReceivingField[3*i+1];          
          receivingField[3*((*_supportMesh->getParentNum())[i]-1)+2] = cpReceivingField[3*i+2];          
          break;

        }
      }
    }

    //
    // Visualization 
    
    _visualization(exchangeName,
                   fieldDimension, 
                   timeStep, 
                   timeValue,
                   sendingFieldName,
                   sendingField,
                   receivingFieldName,
                   receivingField);

    return isReceived;
  }

  void Coupling::_visualization(const char *exchangeName,
                                const couplings_field_dimension_t fieldDimension, 
                                const int timeStep, 
                                const double timeValue,
                                const char  *sendingFieldName,
                                const void *sendingField,
                                const char  *receivingFieldName,
                                const void *receivingField)
  {
    if ((_outputFrequency > 0) && (timeStep % _outputFrequency == 0)) {

      std::string pathString = "./Couplings/"+_name+"/"+
        _localApplicationProperties.getName()+"->"+
        _coupledApplicationProperties.getName();
      
      if (_fvmWriter == NULL) {
        _fvmWriter = fvm_writer_init("Chr",
                                     pathString.c_str(),
                                     _outputFormat.c_str(),
                                     _outputFormatOption.c_str(),
                                     FVM_WRITER_FIXED_MESH); 
        fvm_writer_export_nodal(_fvmWriter, &(_supportMesh->getFvmNodal()));
      }
      
      std::string localName = "S_"+std::string(sendingFieldName);
      fvm_writer_var_loc_t fvm_writer_var_loc;
      fvm_interlace_t fvm_interlace;
      int dim;
      
      if (_solverType == COUPLINGS_SOLVER_CELL_CENTER)
        fvm_writer_var_loc = FVM_WRITER_PER_ELEMENT;
      else
        fvm_writer_var_loc = FVM_WRITER_PER_NODE;
      
      if (fieldDimension == COUPLINGS_FIELD_DIMENSION_SCALAR) {
        dim = 1;
        fvm_interlace = FVM_NO_INTERLACE;
      }
      else if (fieldDimension == COUPLINGS_FIELD_DIMENSION_INTERLACED_VECTOR){
        dim = 3;
        fvm_interlace = FVM_INTERLACE;
      }
      else
        bft_error(__FILE__, __LINE__, 0, "'%i' bad field dimension\n", fieldDimension);
    
      if (sendingField != NULL)

        fvm_writer_export_field(_fvmWriter,
                                const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                                localName.c_str(),
                                fvm_writer_var_loc,
                                1,
                                fvm_interlace,
                                0,
                                NULL,
                                FVM_DOUBLE,
                                timeStep,
                                timeValue,
                                (const void *const *) &sendingField);

      if (receivingField != NULL) {
    
        localName = "R_"+std::string(receivingFieldName);
        fvm_writer_export_field(_fvmWriter,
                                const_cast<fvm_nodal_t *> (&_supportMesh->getFvmNodal()),
                                localName.c_str(),
                                fvm_writer_var_loc,
                                1,
                                FVM_NO_INTERLACE,
                                0,
                                NULL,
                                FVM_DOUBLE,
                                timeStep,
                                timeValue,
                                (const void *const *) &receivingField);
      }
    }
  }
  
  void Coupling::_locate()
  {
    if (_fvmLocator == NULL || _toLocate) {
      
      if (_fvmLocator == NULL) {
        const MPI_Comm& globalComm = _localApplicationProperties.getGlobalComm();
        const int begRank = _coupledApplicationProperties.getBeginningRank();
        const int endRank = _coupledApplicationProperties.getBeginningRank();
        _fvmLocator = fvm_locator_create(_tolerance, globalComm, begRank, endRank);
      }
      
      double* coords;
      if (_coordsPointsToLocate != NULL) 
        coords = _coordsPointsToLocate;
      
      else if(_solverType == COUPLINGS_SOLVER_CELL_CENTER) {
        _nPointsToLocate = _supportMesh->getNElts();
        if (_supportMesh->getParentNum() != NULL) {
          const std::vector<int> & parentNum = *_supportMesh->getParentNum();
          const std::vector<double> & couplingCellCoords = _supportMesh->getCellCenterCoords();
          coords = new double[_nPointsToLocate];
          for (int i = 0; i < _nPointsToLocate; i++) {
            coords[3*(parentNum[3*i]-1)]   = couplingCellCoords[3*i];
            coords[3*(parentNum[3*i]-1)+1] = couplingCellCoords[3*i+1];
            coords[3*(parentNum[3*i]-1)+2] = couplingCellCoords[3*i+2];
          }
        }
        else
          coords = const_cast <double*> (&(_supportMesh->getCellCenterCoords()[0]));
      }
      
      else if(_solverType == COUPLINGS_SOLVER_CELL_VERTEX) {
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
      const int nExterior = fvm_locator_get_n_exterior(_fvmLocator);
      assert(nNotLocatedPoint == nExterior);

      //
      // Renumbering not located points

      if (_supportMesh->getParentNum() != NULL && 
          _solverType == COUPLINGS_SOLVER_CELL_CENTER &&
          _coordsPointsToLocate != NULL){

        const std::vector<int>& parentNum = *_supportMesh->getParentNum();

        if (_notLocatedPoint != NULL && nNotLocatedPoint > _nNotLocatedPoint)
          delete _notLocatedPoint;
        
        if (_notLocatedPoint == NULL)
          _notLocatedPoint = new int[nNotLocatedPoint];
        
        for (int i = 0; i < nNotLocatedPoint; i++) 
          _notLocatedPoint[i] = parentNum[exteriorList[i]-1];
      }
      else
        _notLocatedPoint = const_cast<int *> (exteriorList);

      _nNotLocatedPoint = nNotLocatedPoint;
    }
  }
}

