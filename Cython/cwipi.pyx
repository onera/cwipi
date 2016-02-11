#-----------------------------------------------------------------------------
# This file is part of the CWIPI library. 
#
# Copyright (C) 2011  ONERA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

cimport cython

import numpy as np
cimport numpy as np

cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *
from libc.stdlib cimport malloc, free

interp_f={}
current_cpl = "" 

cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)

cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass

cdef extern from "cwipi.h":

    ctypedef enum cwipi_coupling_type_t:
        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING
        CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING
        CWIPI_COUPLING_SEQUENTIAL

    ctypedef enum cwipi_mesh_type_t:
        CWIPI_STATIC_MESH
        CWIPI_MOBILE_MESH
        CWIPI_CYCLIC_MESH

    ctypedef enum cwipi_solver_type_t:
        CWIPI_SOLVER_CELL_CENTER
        CWIPI_SOLVER_CELL_VERTEX

    ctypedef enum cwipi_located_point_info_t:
        CWIPI_BASIC_INFO
        CWIPI_DISTANT_MESH_INFO

    ctypedef enum cwipi_exchange_status_t:
        CWIPI_EXCHANGE_OK
        CWIPI_EXCHANGE_BAD_RECEIVING

    ctypedef void (*cwipi_interpolation_fct_t) (int entities_dim,
                                               int n_local_vertex,
                                               int n_local_element,
                                               int n_local_polhyedra,
                                               int n_distant_point,
                                               double local_coordinates[],
                                               int local_connectivity_index[],
                                               int local_connectivity[],
                                               int local_polyhedra_face_index[],
                                               int local_polyhedra_cell_to_face_connectivity[],
                                               int local_polyhedra_face_connectivity_index[],
                                               int local_polyhedra_face_connectivity[],
                                               double distant_points_coordinates[],
                                               int distant_points_location[],
                                               float distant_points_distance[],
                                               int distant_points_barycentric_coordinates_index[],
                                               double distant_points_barycentric_coordinates[],
                                               int stride,
                                               cwipi_solver_type_t  solver_type,
                                               void *local_field,
                                               void *distant_field)


#
# CWIPI bases functions

    void cwipi_init(MPI_Comm common_comm, char* application_name, MPI_Comm* application_comm)
    void cwipi_set_output_listing(FILE* output_listing)
    void cwipi_dump_application_properties()
    void cwipi_finalize()

#
# Functions about control parameters

    void cwipi_add_local_int_control_parameter(char* name, int initial_value)
    void cwipi_add_local_double_control_parameter(char* name, double initial_value)
    void cwipi_add_local_string_control_parameter(char* name, char* initial_value)
    void cwipi_set_local_int_control_parameter(char* name, int value)
    void cwipi_set_local_double_control_parameter(char* name, double value)
    void cwipi_set_local_string_control_parameter(char* name, char* value)
    int cwipi_get_local_int_control_parameter(char* name)
    double cwipi_get_local_double_control_parameter(char* name)
    char* cwipi_get_local_string_control_parameter(char* name)
    void cwipi_delete_local_int_control_parameter(char* name)
    void cwipi_delete_local_double_control_parameter(char* name)
    void cwipi_delete_local_string_control_parameter(char* name)
    int cwipi_get_distant_int_control_parameter(char* application_name, char* name)
    double cwipi_get_distant_double_control_parameter(char* application_name, char* name)
    char* cwipi_get_distant_string_control_parameter(char* application_name, char* name)
    void cwipi_synchronize_control_parameter(char* application_name)

#
# Coupling basic function

    void cwipi_create_coupling(char* coupling_name, cwipi_coupling_type_t coupling_type, char* coupled_application,
                               int entitiesDim, double tolerance, cwipi_mesh_type_t mesh_type, cwipi_solver_type_t solver_type, 
                               int output_frequency, char* output_format, char* output_format_option, ...)
    void cwipi_delete_coupling(char* coupling_id)
    void cwipi_set_points_to_locate(char* coupling_id, int n_points, double coordinate[])
    void cwipi_define_mesh(char* coupling_id, int n_vertex, int n_element, double* coordinates, int* connectivity_index, int* connectivity)
#    void cwipi_shared_fvm_nodal(char* coupling_name,
#                                fvm_nodal_t* fvm_nodal)
    void cwipi_add_polyhedra(char* coupling_id, int n_element, int face_index[], int cell_to_face_connectivity[], int n_face, int face_connectivity_index[],
                             int face_connectivity[])
    void cwipi_locate (char* coupling_id)
    void cwipi_update_location (char* coupling_id)
    cwipi_exchange_status_t cwipi_exchange(char* coupling_id, char* exchange_name, int stride, int time_step, double time_value,
                                           char* sending_field_name, double* sending_field, char* receiving_field_name, double* receiving_field, int* n_not_located_points)
    void cwipi_issend(char *coupling_name, char *exchange_name, int tag, int stride,
                      int time_step, double time_value, char *sending_field_name,
                      double *sending_field,int *request)
    void cwipi_irecv(char *coupling_name, char *exchange_name, int tag, int stride,
                     int time_step, double time_value, char *receiving_field_name,
                     double *receiving_field, int *request)
    void cwipi_wait_issend(char *coupling_name, int request)
    void cwipi_wait_irecv(char *coupling_name, int request)
    void cwipi_set_interpolation_function(char* coupling_id, cwipi_interpolation_fct_t fct)
    int* cwipi_get_not_located_points(char* coupling_id)
    int* cwipi_get_located_points(char* coupling_id)
    int cwipi_get_n_located_points(char* coupling_id)
    int cwipi_get_n_not_located_points(char* coupling_id)

    #
    # info about mpi rank, mesh and element where local points are located 

    int* cwipi_get_distant_location (char* coupling_id)
    float* cwipi_get_distant_distance (char* coupling_id)
    double* cwipi_get_distant_coordinates(char* coupling_id)
    int* cwipi_get_distant_barycentric_coordinates_index (char* coupling_id)
    double* cwipi_get_distant_barycentric_coordinates (char* coupling_id)
    int cwipi_get_n_distant_points(char* coupling_id)
    int cwipi_get_n_distant_ranks(char *coupling_id)
    int *cwipi_get_distant_distribution(char *coupling_id)
    int *cwipi_get_located_points_distribution(char *coupling_id)
 
    int cwipi_has_int_parameter(char *application_name, char *name)
    int cwipi_has_double_parameter(char *application_name, char *name)
    int cwipi_has_string_parameter(char *application_name, char *name)

    int cwipi_get_n_int_parameters(char *application_name)
    int cwipi_get_n_double_parameters(char *application_name)
    int cwipi_get_n_string_parameters(char *application_name)

    char ** cwipi_get_list_int_parameters(char *application_name)
    char ** cwipi_get_list_double_parameters(char *application_name)
    char ** cwipi_get_list_string_parameters(char *application_name)

    void cwipi_set_location_index(char *coupling_name, int index)
    void cwipi_load_location(char *coupling_name)
    void cwipi_save_location(char *coupling_name)

    void cwipi_open_location_file(char *coupling_name, char *filename, char *mode)
    void cwipi_close_location_file(char *coupling_name)

COUPLING_PARALLEL_WITH_PARTITIONING = CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING
COUPLING_PARALLEL_WITHOUT_PARTITIONING = CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING
COUPLING_SEQUENTIAL = CWIPI_COUPLING_SEQUENTIAL

STATIC_MESH = CWIPI_STATIC_MESH
MOBILE_MESH = CWIPI_MOBILE_MESH
CYCLIC_MESH = CWIPI_CYCLIC_MESH

SOLVER_CELL_CENTER = CWIPI_SOLVER_CELL_CENTER
SOLVER_CELL_VERTEX = CWIPI_SOLVER_CELL_VERTEX

BASIC_INFO = CWIPI_BASIC_INFO
DISTANT_MESH_INFO = CWIPI_DISTANT_MESH_INFO

EXCHANGE_OK = CWIPI_EXCHANGE_OK
EXCHANGE_BAD_RECEIVING = CWIPI_EXCHANGE_BAD_RECEIVING

#
# Basic fonctions
# ---------------

def  init(MPI.Comm common_comm, char* application_name):
    """
     Initialize the cwipi library and create 
     the current communicator application from 'common_comm'.

     parameters:
       common_comm       <-- Common MPI communicator
       application_name  <-- Current application name
       application_comm  --> Internal MPI communicator for the current
                             application
     It is a synchronization point between all applications
    """
    cdef MPI_Comm c_common_comm = common_comm.ob_mpi

    cdef MPI.Comm application_comm = MPI.Comm()

    cwipi_init(c_common_comm, application_name, &(application_comm.ob_mpi))

    return  application_comm


def set_output_listing(file output_listing):
    """
    Set up the file used for the output listing

    parameters:
      output_listing      <-- Output listing file (C function)
    """
    cdef FILE* c_file = PyFile_AsFile(output_listing)
    cwipi_set_output_listing(c_file)     


def dump_application_properties():
    """
    Dump application properties
    """
    cwipi_dump_application_properties()


def finalize():
    """
    Finalize
    """
    cwipi_finalize()

#
# Control parameters
# ------------------

def add_local_int_control_parameter(char* name, int initial_value):
    """
    Add a integer control parameter

    parameters
      name           <-- parameter name
      initial_value  <-- initial value
    """
    cwipi_add_local_int_control_parameter(name, initial_value)


def add_local_double_control_parameter(char* name, double initial_value):
    """
    Add a double control parameter

    parameters
      name           <-- parameter name
      initial_value  <-- initial value
    """
    cwipi_add_local_double_control_parameter(name, initial_value)
 

def add_local_string_control_parameter(char* name, char* initial_value):
    """
    Add a string control parameter

    parameters
      name           <-- parameter name
      initial_value  <-- initial value
    """
    cwipi_add_local_string_control_parameter(name, initial_value)


def set_local_int_control_parameter(char* name, int value):
    """
    Set a integer control parameter

    parameters
      name   <-- parameter name
      value  <-- value
    """
    cwipi_set_local_int_control_parameter(name, value)


def set_local_double_control_parameter(char* name, double value):
    """
    Set a double control parameter

    parameters
      name   <-- parameter name
      value  <-- value
    """
    cwipi_set_local_double_control_parameter(name, value)


def set_local_string_control_parameter(char* name, char* value):
    """
    Set a string control parameter

    parameters
      name     <-- parameter name
      value    <-- value
    """
    cwipi_set_local_string_control_parameter(name, value)


def get_local_int_control_parameter(char* name):
    """
    Get a local integer control parameter

    parameters
      name           <-- parameter name
    """
    return cwipi_get_local_int_control_parameter(name)


def get_local_double_control_parameter(char* name):
    """
    Get a local double control parameter

    parameters
      name           <-- parameter name
    """
    return cwipi_get_local_double_control_parameter(name)


def get_local_string_control_parameter(char* name):
    """
    Get a local string control parameter

    parameters
      name           <-- parameter name
    """
    return cwipi_get_local_string_control_parameter(name)


def delete_local_int_control_parameter(char* name):
    """
    Delete a local integer control parameter

    parameters
      name           <-- parameter name
    """
    cwipi_delete_local_int_control_parameter(name)


def delete_local_double_control_parameter(char* name):
    """
    Delete a local double control parameter

    parameters
      name           <-- parameter name
    """
    cwipi_delete_local_double_control_parameter(name)


def delete_local_string_control_parameter(char* name):
    """
    Delete a local string control parameter

    parameters
      name           <-- parameter name
    """
    cwipi_delete_local_string_control_parameter(name)


def get_distant_int_control_parameter(char* application_name, char* name):
    """
    Get a distant integer control parameter

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    """
    return cwipi_get_distant_int_control_parameter(application_name, name)


def get_distant_double_control_parameter(char* application_name, char* name):
    """
    Get a distant double control parameter

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    """
    return cwipi_get_distant_double_control_parameter(application_name, name)


def get_distant_string_control_parameter(char* application_name, char* name):
    """
    Get a distant string control parameter

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    """
    return cwipi_get_distant_string_control_parameter(application_name, name)


def has_int_parameter(char* application_name, char* name):
    """
    Has this int parameter ?

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    return
      boolean
    """
    return (cwipi_has_int_parameter(application_name, name) == 1)


def has_double_parameter(char* application_name, char* name):
    """
    Has this double parameter ?

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    return
      boolean
    """
    return (cwipi_has_double_parameter(application_name, name) == 1)


def has_string_parameter(char* application_name, char* name):
    """
    Has this double parameter ?

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    return
      boolean
    """
    return (cwipi_has_string_parameter(application_name, name) == 1)


def get_list_int_parameter(char* application_name):
    """
    return int parameters names

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    return
      list
    """
    i_parameters = []
    cdef int n_parameters = cwipi_get_n_int_parameters(application_name)
    cdef char** c_parameters = cwipi_get_list_int_parameters(application_name)

    for i in range(n_parameters) :
        i_parameters.append(str(c_parameters[i]))
        free(c_parameters[i])

    free(c_parameters)

    return i_parameters


def get_list_double_parameter(char* application_name):
    """
    return double parameters names

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    return
      list
    """
    d_parameters = []
    cdef int n_parameters = cwipi_get_n_double_parameters(application_name)
    cdef char** c_parameters = cwipi_get_list_double_parameters(application_name)

    for i in range(n_parameters) :
        d_parameters.append(str(c_parameters[i]))
        free(c_parameters[i])

    free(c_parameters)
    return d_parameters


def get_list_string_parameter(char* application_name):
    """
    return string parameters names

    parameters
      application_name <-- distant application name
      name             <-- parameter name
    return
      list
    """
    s_parameters = []
    cdef int n_parameters = cwipi_get_n_string_parameters(application_name)
    cdef char** c_parameters = cwipi_get_list_string_parameters(application_name)

    for i in range(n_parameters) :
        s_parameters.append(str(c_parameters[i]))
        free(c_parameters[i])

    free(c_parameters)

    return s_parameters


def synchronize_control_parameter(char* application_name):
    """
    Synchronize local control parameters with an other application.
    It is a synchronization point with this second application

    parameters
      application_name    <-- application name
    """
    cwipi_synchronize_control_parameter(application_name)

#
# Class coupling
# --------------

cdef class Coupling (object):

    """
    Coupling
    """
    cdef char* name

    def __init__(self, 
                 char* coupling_name, 
                 cwipi_coupling_type_t coupling_type, 
                 char *coupled_application, 
                 int entitiesDim,
                 double tolerance, 
                 cwipi_mesh_type_t mesh_type, 
                 cwipi_solver_type_t solver_type,
                 int output_frequency, 
                 char *output_format,
                 char *output_format_option,
                 nb_locations = None):
        """
        Init
        """
        global current_cpl
        self.name = coupling_name
        current_cpl = self.name
        cdef int _nb_locations
        if (nb_locations is None):
           _nb_locations = 1
        else :
           _nb_locations = <int> nb_locations
        cwipi_create_coupling(coupling_name, 
                              coupling_type, 
                              coupled_application,
                              entitiesDim, 
                              tolerance,  
                              mesh_type, 
                              solver_type, 
                              output_frequency, 
                              output_format, 
                              output_format_option,
                              _nb_locations)
        current_cpl = ""

    def __dealloc__(self):
        cwipi_delete_coupling(self.name) 

#
# Class coupling : basic functions
# --------------------------------

    def define_mesh(self, 
                    int n_vertex, 
                    int n_element, 
                    np.ndarray[np.double_t] coordinates not None, 
                    np.ndarray[np.int32_t] connectivity_index not None, 
                    np.ndarray[np.int32_t] connectivity not None):
        """
        Define mesh
        """
        global current_cpl
        current_cpl = self.name
        assert (3 * n_vertex) <= coordinates.size
        assert (n_element + 1) <= connectivity_index.size
        cwipi_define_mesh(self.name,
                          n_vertex, 
                          n_element, 
                          <double *> coordinates.data, 
                          <int *> connectivity_index.data, 
                          <int *> connectivity.data)
        current_cpl = ""

 

    def set_points_to_locate(self,
                             int n_points, 
                             np.ndarray[np.double_t] coordinates not None):
        """
        Set points to locate
        """ 
        global current_cpl
        current_cpl = self.name
        cwipi_set_points_to_locate(self.name,
                                   n_points,
                                   <double *> coordinates.data)
        current_cpl = ""


    def add_polyhedra(self, 
                      int n_element, 
                      np.ndarray[np.int32_t] face_index not None, 
                      np.ndarray[np.int32_t] cell_to_face_connectivity not None, 
                      int n_face, 
                      np.ndarray[np.int32_t] face_connectivity_index not None,
                      np.ndarray[np.int32_t] face_connectivity not None):
        """
        Add polyhedra
        """
        global current_cpl
        current_cpl = self.name
        assert face_index.size == (n_element + 1)
        cwipi_add_polyhedra(self.name, 
                            n_element,
                            <int *> face_index.data, 
                            <int *> cell_to_face_connectivity.data,
                            n_face, 
                            <int *> face_connectivity_index.data,
                            <int *> face_connectivity.data)

        current_cpl = ""


    def locate(self):
        """
        Locate
        """
        global current_cpl
        current_cpl = self.name
        cwipi_locate (self.name)
        current_cpl = ""


    def set_location_index(self, int index):
        """
        Set location index
        """
        global current_cpl
        current_cpl = self.name
        cwipi_set_location_index(self.name, index)
        current_cpl = ""


    def load_location(self):
        """
        Set location index
        """
        global current_cpl
        current_cpl = self.name
        cwipi_load_location(self.name)
        current_cpl = ""


    def save_location(self):
        """
        Set location index
        """
        global current_cpl
        current_cpl = self.name
        cwipi_save_location(self.name)
        current_cpl = ""


    def open_location_file(self, char *filename, char *mode):
        """
        Set location index
        """
        global current_cpl
        current_cpl = self.name
        cwipi_open_location_file(self.name, filename, mode)
        current_cpl = ""


    def close_location_file(self):
        """
        Set location index
        """
        global current_cpl
        current_cpl = self.name
        cwipi_close_location_file(self.name)
        current_cpl = ""


    def update_location(self):
        """
        Locate
        """
        global current_cpl
        current_cpl = self.name
        cwipi_update_location (self.name)
        current_cpl = ""


    def exchange(self, 
                 char* exchange_name, 
                 int stride, 
                 int time_step, 
                 double time_value,
                 char* sending_field_name, 
                 np.ndarray[np.double_t] sending_field, 
                 char* receiving_field_name, 
                 np.ndarray[np.double_t] receiving_field): 
        """
        Exchange
        """
        global current_cpl
        current_cpl = self.name

        cdef int c_n_not_located_points
        cdef cwipi_exchange_status_t status

        if (sending_field is None) and  (receiving_field is not None):
          status = cwipi_exchange(self.name, 
                                   exchange_name, 
                                   stride, 
                                   time_step, 
                                   time_value,
                                   sending_field_name, 
                                   NULL,
                                   receiving_field_name, 
                                   <double*> receiving_field.data, 
                                   &c_n_not_located_points)
        elif (sending_field is not None) and  (receiving_field is None):
          status = cwipi_exchange(self.name, 
                                   exchange_name, 
                                   stride, 
                                   time_step, 
                                   time_value,
                                   sending_field_name, 
                                   <double*> sending_field.data, 
                                   receiving_field_name, 
                                   NULL, 
                                   &c_n_not_located_points)
        elif (sending_field is not None) and  (receiving_field is not None):
          status = cwipi_exchange(self.name, 
                                   exchange_name, 
                                   stride, 
                                   time_step, 
                                   time_value,
                                   sending_field_name, 
                                   <double*> sending_field.data, 
                                   receiving_field_name, 
                                   <double*> receiving_field.data, 
                                   &c_n_not_located_points)
        else :
          status = cwipi_exchange(self.name, 
                                   exchange_name, 
                                   stride, 
                                   time_step, 
                                   time_value,
                                   sending_field_name, 
                                   NULL,
                                   receiving_field_name, 
                                   NULL,
                                   &c_n_not_located_points)

        current_cpl = ""
        return {'status':status, 'n_not_located_points':c_n_not_located_points}


    def issend(self, 
               char* exchange_name,
               int tag, 
               int stride, 
               int time_step, 
               double time_value,
               char* sending_field_name, 
               np.ndarray[np.double_t] sending_field): 
        """
        Issend
        """
        cdef int request

        global current_cpl
        current_cpl = self.name

        if sending_field is None:
            cwipi_issend(self.name, 
                         exchange_name, 
                         tag,
                         stride, 
                         time_step, 
                         time_value,
                         sending_field_name, 
                         NULL,
                         &request)
        else:
            cwipi_issend(self.name, 
                         exchange_name, 
                         tag,
                         stride, 
                         time_step, 
                         time_value,
                         sending_field_name, 
                         <double*> sending_field.data, 
                         &request)


        current_cpl = ""
        return {'request':request}


    def irecv(self, 
              char* exchange_name,
              int tag, 
              int stride, 
              int time_step, 
              double time_value,
              char* receiving_field_name, 
              np.ndarray[np.double_t] receiving_field): 
        """
        Irecv
        """
        cdef int request

        global current_cpl
        current_cpl = self.name

        if receiving_field is None:
            cwipi_irecv(self.name, 
                        exchange_name, 
                        tag,
                        stride, 
                        time_step, 
                        time_value,
                        receiving_field_name, 
                        NULL,
                        &request)
        else:
            cwipi_irecv(self.name, 
                        exchange_name, 
                        tag,
                        stride, 
                        time_step, 
                        time_value,
                        receiving_field_name, 
                        <double*> receiving_field.data, 
                        &request)


        current_cpl = ""
        return {'request':request}


    def wait_issend(self,
                    int request):
        """
        Wait issend
        """
        global current_cpl
        current_cpl = self.name
        cwipi_wait_issend(self.name, request)

        current_cpl = ""


    def wait_irecv(self,
                    int request):
        """
        Wait irecv
        """
        global current_cpl
        current_cpl = self.name
        cwipi_wait_irecv(self.name, request)

        current_cpl = ""


    def get_n_located_points(self):
        """
        Get number of located points
        """
        return cwipi_get_n_located_points(self.name)


    def get_n_not_located_points(self):
        """
        Get number of not located points
        """
        return cwipi_get_n_not_located_points(self.name)


    def get_not_located_points(self):
        """
        Get not located points
        """
        np.import_array()
        cdef np.npy_intp dims = <np.npy_intp> cwipi_get_n_not_located_points(self.name)
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_INT32,
                                                 <void *> cwipi_get_not_located_points(self.name)) 


    def get_located_points(self):
        """
        Get not located points
        """
        np.import_array()
        cdef np.npy_intp dims = <np.npy_intp> cwipi_get_n_located_points(self.name)
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_INT32,
                                                 <void *> cwipi_get_located_points(self.name)) 

    def get_distant_location(self):
        """
        Get distant point location 
        """
        np.import_array()
        cdef np.npy_intp dims = <np.npy_intp> cwipi_get_n_distant_points(self.name)
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_INT32,
                                                 <void *> cwipi_get_distant_location(self.name))

 
    def get_distant_distance(self):
        """
        Get distant points distance to location element
        """
        np.import_array()
        cdef np.npy_intp dims = <np.npy_intp> cwipi_get_n_distant_points(self.name)
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_FLOAT,
                                                 <void *> cwipi_get_distant_distance(self.name))


    def get_distant_coordinates(self):
        """
        Get distant points coordinates
        """
        np.import_array()
        cdef np.npy_intp dims = 3 * <np.npy_intp> cwipi_get_n_distant_points(self.name)
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_DOUBLE,
                                                 <void *> cwipi_get_distant_coordinates(self.name))


    def get_distant_barycentric_coordinates_index(self):
        """
        Get distant points barycentric coordinates index
        """
        np.import_array()
        cdef np.npy_intp dims = <np.npy_intp> cwipi_get_n_distant_points(self.name) + 1
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_INT32,
                                                 <void *> cwipi_get_distant_barycentric_coordinates_index(self.name))


    def get_distant_barycentric_coordinates(self):
        """
        Get distant points barycentric coordinates
        """
        np.import_array()
        cdef np.npy_intp dims1 = <np.npy_intp> cwipi_get_n_distant_points(self.name) + 1
        cdef np.npy_intp dims = <np.npy_intp> (cwipi_get_distant_barycentric_coordinates_index(self.name)[dims1])
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_DOUBLE,
                                                 <void *> cwipi_get_distant_barycentric_coordinates(self.name))


    def get_n_distant_points(self):
        """
        Get number of distant points
        """
        return  cwipi_get_n_distant_points(self.name)


    def get_n_distant_ranks(self):
        """
        Get number of distant ranks
        """
        return  cwipi_get_n_distant_ranks(self.name)


    def get_distant_distribution(self):
        """
        Get distant points distribution on distant ranks
        """
        np.import_array()
        cdef np.npy_intp dims = <np.npy_intp> cwipi_get_n_distant_ranks(self.name) + 1
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_INT32,
                                                 <void *> cwipi_get_distant_distribution(self.name))


    def get_located_points_distribution(self):
        """
        Get located points distribution on distant ranks
        """
        np.import_array()
        cdef np.npy_intp dims = <np.npy_intp> cwipi_get_n_distant_ranks(self.name) + 1
        if (dims == 0):
            return None
        else :
            return np.PyArray_SimpleNewFromData(1, 
                                                 &dims, 
                                                 np.NPY_INT32,
                                                 <void *> cwipi_get_located_points_distribution(self.name))


    def set_interpolation_function(self, f):
        """
        """
        global current_cpl
        global interp_f
        current_cpl = self.name
        interp_f[self.name]=f
        cwipi_set_interpolation_function(self.name, callback)
        current_cpl = ""


    def getName(self):
        """
        """
        return self.name

#
# Class coupling : additional functions (activated by set_info function)
# ----------------------------------------------------------------------

# TODO:
  
#
# info about mpi rank, mesh and element where local points are located 

#    int* cwipi_get_distant_location (char* coupling_id)
#    double* cwipi_get_distant_coordinates(char* coupling_id)
#    int* cwipi_get_distant_barycentric_coordinates_index (char* coupling_id)
#    double* cwipi_get_distant_barycentric_coordinates (char* coupling_id)
#    int cwipi_get_n_distant_points(char* coupling_id)


cdef void callback(int entities_dim,
                   int n_local_vertex,
                   int n_local_element,
                   int n_local_polhyedra,
                   int n_distant_point,
                   double local_coordinates[],
                   int local_connectivity_index[],
                   int local_connectivity[],
                   int local_polyhedra_face_index[],
                   int local_polyhedra_cell_to_face_connectivity[],
                   int local_polyhedra_face_connectivity_index[],
                   int local_polyhedra_face_connectivity[],
                   double distant_points_coordinates[],
                   int distant_points_location[],
                   float distant_points_distance[],
                   int distant_points_barycentric_coordinates_index[],
                   double distant_points_barycentric_coordinates[],
                   int stride,
                   cwipi_solver_type_t  solver_type,
                   void *local_field,
                   void *distant_field):
    """
    """
    global current_cpl
    global inter_f
    np.import_array()

    cdef np.npy_intp dims = 0

    if (local_coordinates == NULL or n_local_vertex == 0):
        local_coordinates_a = None
    else:   
        dims = <np.npy_intp>(3 * n_local_vertex)
        local_coordinates_a = np.PyArray_SimpleNewFromData(1, 
                                                       &dims, 
                                                       np.NPY_DOUBLE,
                                                       <void *> local_coordinates) 

    if (local_connectivity_index == NULL or n_local_element == 0):
        local_connectivity_index_a = None 
        local_connectivity_a = None
    else:
        dims = <np.npy_intp>(n_local_element + 1)    
        local_connectivity_index_a = np.PyArray_SimpleNewFromData(1, 
                                                                  &dims, 
                                                                  np.NPY_INT32,
                                                                  <void *> local_connectivity_index) 
        dims = <np.npy_intp>(local_connectivity_index[n_local_element])    
        local_connectivity_a = np.PyArray_SimpleNewFromData(1, 
                                                            &dims, 
                                                            np.NPY_INT32,
                                                            <void *> local_connectivity) 

    if (local_polyhedra_face_index == NULL or n_local_polhyedra == 0):
        local_polyhedra_face_index_a = None
        local_polyhedra_cell_to_face_a = None
        local_polyhedra_face_connectivity_index_a = None
        local_polyhedra_face_connectivity_a = None
    else:
        dims = <np.npy_intp>(n_local_polhyedra+1)    
        local_polyhedra_face_index_a = np.PyArray_SimpleNewFromData(1, 
                                                                    &dims, 
                                                                    np.NPY_INT32,
                                                                    <void *> local_polyhedra_face_index) 

        dims = <np.npy_intp>(local_polyhedra_face_index[n_local_polhyedra])    
        local_polyhedra_cell_to_face_a = np.PyArray_SimpleNewFromData(1, 
                                                          &dims, 
                                                          np.NPY_INT32,
                                                          <void *> local_polyhedra_cell_to_face_connectivity) 

        dims = <np.npy_intp> (0)
        for i in local_polyhedra_cell_to_face_a:
            dims = <np.npy_intp> max(dims, i)
        dims = dims + 1
        local_polyhedra_face_connectivity_index_a = np.PyArray_SimpleNewFromData(1, 
                                                                                 &dims, 
                                                                                 np.NPY_INT32,
                                                                                 <void *> local_polyhedra_face_connectivity_index) 
    
        dims = <np.npy_intp> (local_polyhedra_cell_to_face_connectivity[<int>dims -1])
        local_polyhedra_face_connectivity_a = np.PyArray_SimpleNewFromData(1, 
                                                                          &dims, 
                                                                          np.NPY_INT32,
                                                                         <void *> local_polyhedra_face_connectivity) 

    if (distant_points_location == NULL or n_distant_point == 0):
        distant_points_location_a = None
        distant_points_coordinates_a = None
        distant_points_distance_a = None
        distant_points_barycentric_coordinates_index_a = None
        distant_points_barycentric_coordinates_a = None
    else:

        dims = <np.npy_intp>(n_distant_point)    
        distant_points_location_a = np.PyArray_SimpleNewFromData(1, 
                                                                 &dims, 
                                                                 np.NPY_INT,
                                                                 <void *> distant_points_location) 

        distant_points_distance_a = np.PyArray_SimpleNewFromData(1, 
                                                                 &dims, 
                                                                 np.NPY_FLOAT,
                                                                 <void *> distant_points_distance) 
        dims = <np.npy_intp>(3 * n_distant_point)    
        distant_points_coordinates_a = np.PyArray_SimpleNewFromData(1, 
                                                                    &dims, 
                                                                    np.NPY_DOUBLE,
                                                                    <void *> distant_points_coordinates) 
        if (distant_points_barycentric_coordinates_index == NULL or distant_points_barycentric_coordinates == NULL):
            distant_points_barycentric_coordinates_index_a = None
            distant_points_barycentric_coordinates_a = None
        else:
            dims = <np.npy_intp>(n_distant_point + 1)    
            distant_points_barycentric_coordinates_index_a = np.PyArray_SimpleNewFromData(1, 
                                                                                          &dims, 
                                                                                          np.NPY_INT32,
                                                                                          <void *>distant_points_barycentric_coordinates_index) 

            dims = <np.npy_intp>(distant_points_barycentric_coordinates_index[n_distant_point])    
            distant_points_barycentric_coordinates_a = np.PyArray_SimpleNewFromData(1, 
                                                                                    &dims, 
                                                                                    np.NPY_DOUBLE,
                                                                                    <void *>distant_points_barycentric_coordinates) 


    dims = <np.npy_intp>(0)
    if (solver_type == SOLVER_CELL_CENTER):
        dims = <np.npy_intp> (n_local_element * stride)
    elif (solver_type == SOLVER_CELL_VERTEX):
        dims = <np.npy_intp> (n_local_vertex * stride)

    if (dims == 0 or local_field == NULL):
        local_field_a = None
    else:
        local_field_a = np.PyArray_SimpleNewFromData(1, 
                                                     &dims, 
                                                     np.NPY_DOUBLE,
                                                     <void *> local_field) 

    dims = <np.npy_intp>(0)
    if (solver_type == SOLVER_CELL_CENTER):
        dims = <np.npy_intp> (n_distant_point * stride)
    elif (solver_type == SOLVER_CELL_VERTEX):
        dims = <np.npy_intp> (n_distant_point * stride)

    if (dims == 0 or distant_field == NULL):
        distant_field_a = None
    else:
        distant_field_a = np.PyArray_SimpleNewFromData(1, 
                                                       &dims, 
                                                       np.NPY_DOUBLE,
                                                       <void *> distant_field) 

    (<object>interp_f[current_cpl])(entities_dim,
                                    n_local_vertex,
                                    n_local_element,
                                    n_local_polhyedra,
                                    n_distant_point,
                                    local_coordinates_a,
                                    local_connectivity_index_a,
                                    local_connectivity_a,
                                    local_polyhedra_face_index_a,
                                    local_polyhedra_cell_to_face_a,
                                    local_polyhedra_face_connectivity_index_a,
                                    local_polyhedra_face_connectivity_a,
                                    distant_points_coordinates_a,
                                    distant_points_location_a,
                                    distant_points_distance_a,
                                    distant_points_barycentric_coordinates_index_a,
                                    distant_points_barycentric_coordinates_a,
                                    stride,
                                    solver_type,
                                    local_field_a,
                                    distant_field_a)

