#ifndef __COUPLING_PROPERTIES_H__
#define __COUPLING_PROPERTIES_H__

#include <string>
#include <vector>

#include <fvm_locator.h>
#include <fvm_writer.h>
#include <fvm_nodal.h>
#include "couplings.h"

namespace couplings {

  class ApplicationProperties;

  class Mesh;

  class Coupling {

  public:

    Coupling(const std::string& name,
             const couplings_coupling_type_t couplingType,
             const ApplicationProperties& localApplicationProperties,
             const ApplicationProperties& coupledApplicationProperties,
             const int entitiesDim,
             const double tolerance,
             const couplings_solver_type_t solverType,
             const int    outputFrequency,
             const char  *outputFormat,
             const char  *outputFormatOption);

    virtual ~Coupling();

    void defineMesh(const int nVertex,
                    const int nElement,
                    const double coordinates[],
                    int connectivity_index[],
                    int connectivity[]);

    void defineMeshAddPolyhedra(const int n_element,
                                int face_index[],
                                int cell_to_face_connectivity[],
                                int face_connectivity_index[],
                                int face_connectivity[]);

    couplings_exchange_status_t exchange(const char                          *exchange_name,
                                         const int                            stride,
                                         const int                            time_step,
                                         const double                         time_value,
                                         const char                          *sending_field_name,
                                         const double                        *sending_field,
                                         char                                *receiving_field_name,
                                         double                              *receiving_field,
                                         void                                *ptFortranInterpolationFct);

    void updateLocation();

    void setPointsToLocate(const int n_points,
                           double coordinate[]);

    inline void set_interpolation_function(couplings_interpolation_fct_t *fct);

    inline const int & getNNotlocatedPoint() const;

    inline const int *getNotlocatedPoint() const;

    inline int getNLocatedPoint() const;

    inline const int *getLocatedPoint() const;

    inline const int *getDistantLocation() const;

    inline int getNDistantPoint() const;

    inline const int *getDistantBarycentricCoordinatesIndex() const;

    inline const double *getDistantBarycentricCoordinates() const;

    void locate();

  private:

    Coupling();

    Coupling &operator=(const Coupling &other);

     std::vector<double> & _extrapolate(double *cellCenterField);

    //
    // TODO: Dans l'avenir créer une fabrique abstraite qui permet de definir differentes methodes d'interpolation

    void _interpolate(double *vertexField,
                      std::vector<double>& interpolatedField,
                      const int stride);

    void _interpolate1D(double *vertexField,
                        std::vector<double>& interpolatedField,
                        const int stride);

    void _interpolate2D(double *vertexField,
                        std::vector<double>& interpolatedField,
                        const int stride);

    void _interpolate3D(double *vertexField,
                        std::vector<double>& interpolatedField,
                        const int stride);

    void _fieldsVisualization(const char *exchangeName,
                              const int stride,
                              const int timeStep,
                              const double timeValue,
                              const char  *sendingFieldName,
                              const void *sendingField,
                              const char  *receivingFieldName,
                              const void *receivingField);

    void _initVisualization();

    void _createCouplingComm();

  private:
    const std::string   _name;
    const couplings_coupling_type_t _couplingType;
    const ApplicationProperties& _localApplicationProperties;
    const ApplicationProperties& _coupledApplicationProperties;
    const int            _entitiesDim;
    const double            _tolerance;
    const couplings_solver_type_t  _solverType;
    const std::string   _outputFormat;
    const std::string   _outputFormatOption;
    const int           _outputFrequency;
  private:
    int                  _nPointsToLocate;
    Mesh                *_supportMesh;
    double              *_coordsPointsToLocate;
    fvm_locator_t       *_fvmLocator;
    fvm_writer_t        *_fvmWriter;
    MPI_Comm             _couplingComm;
    int                  _coupledApplicationNRankCouplingComm;
    int                  _coupledApplicationBeginningRankCouplingComm;
    bool                 _isCoupledRank;
    MPI_Comm             _mergeComm;
    MPI_Comm             _fvmComm;
    couplings_interpolation_fct_t *_interpolationFct;
    bool                 _toLocate;
    int                 *_barycentricCoordinatesIndex;
    double              *_barycentricCoordinates;
    int                  _nNotLocatedPoint;
    int                 *_notLocatedPoint;
    int                  _nDistantpoint;
    int                 *_locatedPoint;
    int                 *_location;
  private:
    std::vector<double> *_tmpVertexField; //Evite une allocation a chaque extrapolation
    std::vector<double> *_tmpDistantField;
  };

}

#endif //__COUPLING_PROPERTIES_H__
