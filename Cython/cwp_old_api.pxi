# cython: c_string_type=str, c_string_encoding=ascii
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

"""
cwipi - Coupling With Interpolation Parallel Interface library.
"""


# 2018-06-15: python3 modifications:
#   - FILE* no longer available
#   - str vs. bytes (1st line of this file)
#   - tiny API modifications (irecv, isend return request directly)
#   - const in callback (otherwise needs -fpermissive)
#   -  __dealloc__ -> __del__
#   - some safety checks (contiguity, orphan nodes)
#   - add __version__
#   - docstring

import numpy as np
cimport numpy as np

cimport mpi4py.MPI as MPI
# from mpi4py.@mpi4pylibmpi@ cimport *

from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE, fdopen
from cpython.object cimport PyObject_AsFileDescriptor

# fixme: internal global variables: shoud _prefix them, or even remove them?
interp_f={}
interp_ho_loc_f={}
interp_ho_bas_f={}
current_cpl = ""

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

    ctypedef enum cwipi_element_t:
        CWIPI_NODE
        CWIPI_EDGE2
        CWIPI_EDGEHO
        CWIPI_FACE_TRIA3
        CWIPI_FACE_TRIAHO
        CWIPI_FACE_QUAD4
        CWIPI_FACE_QUADHO
        CWIPI_FACE_POLY
        CWIPI_CELL_TETRA4
        CWIPI_CELL_TETRAHO
        CWIPI_CELL_HEXA8
        CWIPI_CELL_HEXAHO
        CWIPI_CELL_PRISM6
        CWIPI_CELL_PRISMHO
        CWIPI_CELL_PYRAM5
        CWIPI_CELL_PYRAMHO
        CWIPI_CELL_POLY


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


    ctypedef void (*cwipi_user_interp_ho_fct_t) (int entities_dim,
                                                 int order,
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
                                                 int distant_points_weights_index[],
                                                 double distant_points_weights[],
                                                 double distant_points_uvw[],
                                                 int stride,
                                                 cwipi_solver_type_t solver_type,
                                                 void *local_field,
                                                 void *distant_field)


    ctypedef double (*cwipi_ho_location_fct_t) (int entities_dim,
                                                int order,
                                                int n_nodes,
                                                double *nodes_coords,
                                                double *point_coords,
                                                double *projected_coords,
                                                double *projected_uvw)


    ctypedef void (*cwipi_ho_basis_fct_t) (int entities_dim,
                                           int order,
                                           int n_nodes,
                                           int n_vtx,
                                           double *uvw,
                                           double *weights)


#
# CWIPI bases functions

    void cwipi_init(MPI.MPI_Comm common_comm, char* application_name, MPI.MPI_Comm* application_comm)
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
    void cwipi_ho_define_mesh(char *coupling_id, int n_vertex, int n_element, int order, double *coordinates, int *connectivity_index, int *connectivity)
    void cwipi_ho_options_set (char *coupling_id, char *option, char *value)
    void cwipi_ho_ordering_from_IJK_set (char *coupling_id, cwipi_element_t t_elt, int n_nodes, int *uvw_grid)
    void cwipi_ho_user_elt_set (cwipi_element_t elt_type, cwipi_ho_basis_fct_t element_basis, cwipi_ho_location_fct_t location_in_element)
    void cwipi_ho_set_interpolation_function (char *coupling_id, cwipi_user_interp_ho_fct_t fct)

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

def init(MPI.Comm common_comm, char* application_name):
    """init(common_comm, application_name)
    Initialize the cwipi library and create the current communicator application from ``common_comm``.

    Parameters:
        common_comm(MPI.Comm):   Common MPI communicator between all applications
        application_name(str):   Current local application name

    Returns:
        application_comm (MPI.Comm): Internal MPI communicator for the current application group.

    It is a synchronization point between all applications.
    """
    cdef MPI.MPI_Comm c_common_comm = common_comm.ob_mpi

    cdef MPI.Comm application_comm = MPI.Comm()

    cwipi_init(c_common_comm, application_name, &(application_comm.ob_mpi))

    return  application_comm


def set_output_listing(output_listing):
    """set_output_listing(output_listing)
    Set up the file used for the output listing.

    Parameters:
        output_listing(file or str):  Output listing file

    Note:
        If this file is written by both cwipi and python, I/O buffering may mix up the output.
        In this case, calling output_listing.flush() before any cwipi call may help.
    """
    cdef int fd
    cdef FILE* c_file
    cdef str mode
    if (isinstance(output_listing, str)):
        output_listing = open(output_listing, 'w')
    # next line will most likely fail if output_listing is not
    # an io.TextIOWrapper(in py3) or a file(in py2)
    mode = output_listing.mode
    # should test if writable?

    fd = PyObject_AsFileDescriptor(output_listing)
    c_file = fdopen(fd, mode)
    cwipi_set_output_listing(c_file)


def dump_application_properties():
    """
    Dump application properties.
    """
    cwipi_dump_application_properties()


def finalize():
    """
    Finalize this module.

    After this call, no other function from this module can be called.
    It is a synchronization point between all applications.
    """
    cwipi_finalize()

#
# Control parameters
# ------------------

## Pythonic function:
def add_local_parameter(name, value):
    """add_local_parameter(name, value)
    Add an integer, double or string parameter to the coupling.

    Parameters:
      name (str) : parameter's name
      value (int, float, str): its initial value

    Note:
      The parameter's type is guessed from ``type(value)``. So the following are different:

        >>> add_local_parameter('time', 10)
        >>> add_local_parameter('time', 10.0)
    """
    if isinstance(value, int):
        add_local_int_control_parameter(name, value)
    elif isinstance(value, float):
        add_local_double_control_parameter(name, value)
    elif isinstance(value, str) or isinstance(value, bytes):  # fixme: bytes ou str ? ou les deux ?
        add_local_string_control_parameter(name, value)
    else:
        raise TypeError("Argument 'value' should be a string, an integer or a float, not '%s'"%type(value))

def set_local_parameter(name, value):
    """set_local_parameter(name, value)
    Set a parameter's value.

    Parameters:
      name (str) : an (existing) parameter's name
      value (int, float, str): its new value

    See also: :func:`add_local_parameter`.
    """
    if isinstance(value, int):
        set_local_int_control_parameter(name, value)
    elif isinstance(value, float):
        set_local_double_control_parameter(name, value)
    elif isinstance(value, str) or isinstance(value, bytes):  # fixme: bytes ou str ? ou les deux ?
        set_local_string_control_parameter(name, value)
    else:
        raise TypeError("Argument 'value' should be a string, an integer or a float, not '%s'"%type(value))

# Even more pythonic access would be with a dictionary(?). So we could retrieve the current value of a parameter?
# For example:
# print (cwipi.local_control_parameters)  # --> {}
# cpl.local_control_parameters['time'] = 2.0
# cpl.local_control_parameters['iteration'] = 200
# cpl.local_control_parameters['quantity'] = 'pressure'
# cpl.synchronize_control_parameters()    # now cpl.distant_control_parameters is up-to-date
# print (cpl.global_control_parameters['time'])  # --> whatever distant put in it
# print (cpl.global_control_parameters['foobar'])  # raise KeyError('No such dist param')
## even better, derive a dict, that has in addition 1 function and 1 attr:
# cpl.local_control_parameters['iteration'] = 300
# print (cwipi.local_control_parameters)  # {'time':20, 'iteration': 300, 'quantity': 'pressure'}
# cpl.control_parameters.synchronized     # is False, because I just changed it
## and same for distant_control_parameters

## Low-level accesses

def add_local_int_control_parameter(char* name, int initial_value):
    """add_local_int_control_parameter(name, initial_value)
    """
    cwipi_add_local_int_control_parameter(name, initial_value)


def add_local_double_control_parameter(char* name, double initial_value):
    """add_local_double_control_parameter(name, initial_value)
    """
    cwipi_add_local_double_control_parameter(name, initial_value)


def add_local_string_control_parameter(char* name, char* initial_value):
    """add_local_string_control_parameter(name, initial_value)
    """
    cwipi_add_local_string_control_parameter(name, initial_value)


def set_local_int_control_parameter(char* name, int value):
    """set_local_int_control_parameter(name, value)
    """
    cwipi_set_local_int_control_parameter(name, value)


def set_local_double_control_parameter(char* name, double value):
    """set_local_double_control_parameter(name, value)
    """
    cwipi_set_local_double_control_parameter(name, value)


def set_local_string_control_parameter(char* name, char* value):
    """set_local_string_control_parameter(name, value)
    """
    cwipi_set_local_string_control_parameter(name, value)


def get_local_int_control_parameter(char* name):
    """get_local_int_control_parameter(name)
    Return the local integer control parameter called `name`.
    """
    return cwipi_get_local_int_control_parameter(name)


def get_local_double_control_parameter(char* name):
    """get_local_double_control_parameter(name)
    Return the local double control parameter called `name`.
    """
    return cwipi_get_local_double_control_parameter(name)


def get_local_string_control_parameter(char* name):
    """get_local_string_control_parameter(name)
    Return the local string control parameter called `name`.
    """
    return cwipi_get_local_string_control_parameter(name)


def delete_local_int_control_parameter(char* name):
    """delete_local_int_control_parameter(name)
    Delete a local integer control parameter.
    """
    cwipi_delete_local_int_control_parameter(name)


def delete_local_double_control_parameter(char* name):
    """delete_local_double_control_parameter(name)
    Delete a local double control parameter.
    """
    cwipi_delete_local_double_control_parameter(name)


def delete_local_string_control_parameter(char* name):
    """delete_local_string_control_parameter(name)
    Delete a local string control parameter.
    """
    cwipi_delete_local_string_control_parameter(name)


def get_distant_int_control_parameter(char* application_name, char* name):
    """get_distant_int_control_parameter(application_name, name)
    Return a distant integer control parameter.
    """
    return cwipi_get_distant_int_control_parameter(application_name, name)


def get_distant_double_control_parameter(char* application_name, char* name):
    """get_distant_double_control_parameter(application_name, name)
    Return a distant double control parameter.
    """
    return cwipi_get_distant_double_control_parameter(application_name, name)


def get_distant_string_control_parameter(char* application_name, char* name):
    """get_distant_string_control_parameter(application_name, name)
    Return a distant string control parameter.
    """
    return cwipi_get_distant_string_control_parameter(application_name, name)


def has_int_parameter(char* application_name, char* name):
    """has_int_parameter(application_name, name)
    Return whether the distant application has this integer parameter.
    """
    return (cwipi_has_int_parameter(application_name, name) == 1)


def has_double_parameter(char* application_name, char* name):
    """has_double_parameter(application_name, name)
    Return whether the distant application has this double parameter.
    """
    return (cwipi_has_double_parameter(application_name, name) == 1)


def has_string_parameter(char* application_name, char* name):
    """has_string_parameter(application_name, name)
    Return whether the distant application has this string parameter.
    """
    return (cwipi_has_string_parameter(application_name, name) == 1)


def get_list_int_parameter(char* application_name):
    """get_list_int_parameter(application_name)
    Return the list of integer parameters in the given distant application.
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
    """get_list_double_parameter(application_name)
    Return the list of double parameters in the given distant application.
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
    """get_list_string_parameter(application_name)
    Return the list of string parameters in the given distant application.
    """
    s_parameters = []
    cdef int n_parameters = cwipi_get_n_string_parameters(application_name)
    cdef char** c_parameters = cwipi_get_list_string_parameters(application_name)

    for i in range(n_parameters) :
        s_parameters.append(str(c_parameters[i]))
        free(c_parameters[i])

    free(c_parameters)

    return s_parameters

# fixme: sould be plural : synchronize_control_parameters
def synchronize_control_parameter(char* application_name):
    """synchronize_control_parameter(application_name)
    Synchronize local control parameters with an other application.

    It is a synchronization point with this second application.
    """
    cwipi_synchronize_control_parameter(application_name)

#
# Class coupling
# --------------

cdef class Coupling (object):

    """
    Create a CWIPI coupling object.
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
        """__init__(coupling_name, coupling_type, coupled_application, entitiesDim, tolerance, mesh_type, solver_type, output_frequency, output_format, output_format_option, nb_locations=None)

        Args:
            coupling_name(str): a name for this object
            coupling_type(cwipi_coupling_type_t): either one of: COUPLING_PARALLEL_WITH_PARTITIONING, COUPLING_PARALLEL_WITHOUT_PARTITIONING, COUPLING_SEQUENTIAL
            coupled_application(str): the identifier of the distant code
            entitiesDim(int): dimension of the declared geometry
            tolerance(double): geometric matching tolerance
            mesh_type(cwipi_mesh_type_t): either one of: STATIC_MESH, MOBILE_MESH, CYCLIC_MESH
            solver_type(cwipi_solver_type_t): either one of: SOLVER_CELL_CENTER, SOLVER_CELL_VERTEX
            output_frequency(int): output frequency (0 disables output)
            output_format(str): either one of "EnSight Gold", "MED_fichier", "CGNS"
            output_format_option(str): extra output options
            nb_locations(int, optional):
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

    # Need a __del__, not __dealloc__ , because it requires self.name, that is defined in __init__
    #      (and in fact __dealloc__ is the counterpart of __cinit__)
    # def __dealloc__(self):
    #     cwipi_delete_coupling(self.name)
    def __del__(self):
        cwipi_delete_coupling(self.name)
    # Anyways, wouldn't it be better for this .pyx file to handle the C++ pointer, without going through the C-API?


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
        define_mesh(n_vertex, n_element, coordinates, connectivity_index, connectivity)
        Define the mesh: nodes, elements and connectivity.

        Args:
            n_vertex(int): number of vertices/nodes
            n_element(int): number of elements
            coordinates(double ndarray): the nodes coordinates, an array whose shape is either (n_vertex*3) or (n_vertex, 3)
            connectivity_index (int32 ndarray): the index of the first node of each element in ``connectivity``; shape is ``(n_element+1)``, indices start at zero.
            connectivity (int32 ndarray): the array defining elements' connectivity; indices start at 1 (i.e. 1 is the first node).
               For example the nodes of element ``iel`` are given by the slice ``connectivity[connectivity_index[iel]:connectivity_index[iel+1]]``.
               The array's size is the sum over elements of el.nb_nodes.
        """
        global current_cpl
        current_cpl = self.name
        assert (3 * n_vertex) <= coordinates.size
        assert (n_element + 1) <= connectivity_index.size, "(%d) <= (%d)"%((n_element + 1), connectivity_index.size)

        ## safety verifications:
        # n_vertex = 0 may indeed happen, is even recommended to provide idle processes for load balancing
        if n_vertex >1:
            # Contiguity
            assert coordinates.flags['C_CONTIGUOUS'], "coordinates must be C-contiguous"
            assert connectivity_index.flags['C_CONTIGUOUS'], "connectivity_index must be C-contiguous"
            assert connectivity.flags['C_CONTIGUOUS'], "connectivity must be C-contiguous"
            # Check for orphan nodes, which currently (0.9.5) segfaults in cpl.locate
            assert (np.min(connectivity) == 1), "connectivity values must start at 1"
            assert (np.max(connectivity) == n_vertex), "connectivity values must be <= n_vertex"
            assert (len(np.unique(connectivity)) == n_vertex), "Cannot have orphan vertex (current limitation)"

        cwipi_define_mesh(self.name,
                          n_vertex,
                          n_element,
                          <double *> coordinates.data,
                          <int *> connectivity_index.data,
                          <int *> connectivity.data)
        current_cpl = ""


    def ho_define_mesh(self,
                       int n_vertex,
                       int n_element,
                       int order,
                       np.ndarray[np.double_t] coordinates not None,
                       np.ndarray[np.int32_t] connectivity_index not None,
                       np.ndarray[np.int32_t] connectivity not None):
        """ho_define_mesh(n_vertex, n_element, order, coordinates, connectivity_index, connectivity)
        Define high order mesh.

        Args:
            order(int): the mesh's interpolation order.
        All other arguments are the same as in :func:`define_mesh`.
        """
        global current_cpl
        current_cpl = self.name
        assert (3 * n_vertex) <= coordinates.size
        assert (n_element + 1) <= connectivity_index.size
        cwipi_ho_define_mesh(self.name,
                             n_vertex,
                             n_element,
                             order,
                             <double *> coordinates.data,
                             <int *> connectivity_index.data,
                             <int *> connectivity.data)
        current_cpl = ""


    def ho_options_set(self,
                       char *option,
                       char *value):
        """
        Set option for high order meshes
        """
        global current_cpl
        current_cpl = self.name
        cwipi_ho_options_set(self.name,
                             option,
                             value)
        current_cpl = ""


    def ho_ordering_from_IJK_set(self,
                                 cwipi_element_t t_elt,
                                 int n_nodes,
                                 np.ndarray[np.int32_t] ijk_grid not None):
        """
        Define ho element ordering from the location in the (I, J, K) grid.
        """
        global current_cpl
        current_cpl = self.name
        cwipi_ho_ordering_from_IJK_set(self.name,
                                       <cwipi_element_t> t_elt,
                                       n_nodes,
                                        <int *> ijk_grid.data)
        current_cpl = ""


    def set_points_to_locate(self,
                             int n_points,
                             np.ndarray[np.double_t] coordinates not None):
        """set_points_to_locate(n_points, coordinates)

        Set points to locate, e.g. the receiving points.

        Args:
           n_points (int): number of points
           coords (double ndarray): the coordinates, an array whose shape is (n_points, dim), with C-ordering.
        """
        global current_cpl
        current_cpl = self.name
        assert (3 * n_points) <= coordinates.size
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
        """locate()
        Start the localization process (this is a synchronization point).
        """
        global current_cpl
        current_cpl = self.name
        cwipi_locate (self.name)
        current_cpl = ""


    def set_location_index(self, int index):
        """
        Set location index.
        """
        global current_cpl
        current_cpl = self.name
        cwipi_set_location_index(self.name, index)
        current_cpl = ""


    def load_location(self):
        """
        """
        global current_cpl
        current_cpl = self.name
        cwipi_load_location(self.name)
        current_cpl = ""


    def save_location(self):
        """
        """
        global current_cpl
        current_cpl = self.name
        cwipi_save_location(self.name)
        current_cpl = ""


    def open_location_file(self, char *filename, char *mode):
        """open_location_file(char *filename, char *mode)
        """
        global current_cpl
        current_cpl = self.name
        cwipi_open_location_file(self.name, filename, mode)
        current_cpl = ""


    def close_location_file(self):
        """
        """
        global current_cpl
        current_cpl = self.name
        cwipi_close_location_file(self.name)
        current_cpl = ""


    def update_location(self):
        """
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
        """exchange(exchange_name, stride, time_step, time_value, sending_field_name, np.ndarray[np.double_t] sending_field, receiving_field_name, np.ndarray[np.double_t] receiving_field)
        Exchange quantities.
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
        Non blocking send.

        Returns:
            the request (that you should later pass to :func:`wait_issend`).

        Args:
            exchange_name(str):
            tag(int): an arbitraty MPI tag, that should match the corresponding recv.
            stride(int): the stride, i.e. number of quantities received.
            time_step(int): the time step (only used in output files).
            time_value(double): the time (only used in output files).
            sending_field_name(str): a name for this quantity (only used in output files).
            sending_field(double ndarray): the data to send.
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
        return request


    def irecv(self,
              char* exchange_name,
              int tag,
              int stride,
              int time_step,
              double time_value,
              char* receiving_field_name,
              np.ndarray[np.double_t] receiving_field):
        """irecv(exchange_name, tag, stride, time_step, time_value, receiving_field_name, receiving_field)
        Non blocking receive.

        Returns:
            the request (that you should later pass to :func:`wait_irecv`).

        Args:
            exchange_name(str):
            tag(int): an arbitraty MPI tag, that should match the corresponding send.
            stride(int): the stride, i.e. number of quantities received.
            time_step(int): the time step (only used in output files).
            time_value(double): the time (only used in output files).
            receiving_field_name(str): a name for this quantity (only used in output files).
            receiving_field(double ndarray): an array large enough (size should be at least n*stride, where n is the number of receiving points).
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
        return request


    def wait_issend(self,
                    int request):
        """
        wait_issend(request)
        """
        global current_cpl
        current_cpl = self.name
        cwipi_wait_issend(self.name, request)

        current_cpl = ""


    def wait_irecv(self,
                    int request):
        """
        wait_irecv(request)
        """
        global current_cpl
        current_cpl = self.name
        cwipi_wait_irecv(self.name, request)

        current_cpl = ""


    def get_n_located_points(self):
        """
        Get number of located points.
        """
        return cwipi_get_n_located_points(self.name)


    def get_n_not_located_points(self):
        """
        Get number of not located points.
        """
        return cwipi_get_n_not_located_points(self.name)


    def get_not_located_points(self):
        """
        Get not located points ranks.
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
        Get located points ranks.
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
        Get distant point location.
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
        Get distant points distance to location element.
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
        Get distant points coordinates.
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
        Get distant points barycentric coordinates index.
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
        Get distant points barycentric coordinates.
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
        Get number of distant points.
        """
        return  cwipi_get_n_distant_points(self.name)


    def get_n_distant_ranks(self):
        """
        Get number of distant ranks.
        """
        return  cwipi_get_n_distant_ranks(self.name)


    def get_distant_distribution(self):
        """
        Get distant points distribution on distant ranks.
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
        Get located points distribution on distant ranks.
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


    def set_ho_interpolation_function(self, f):
        """
        """
        global current_cpl
        global interp_ho_f
        current_cpl = self.name
        interp_ho_f[self.name]=f
        cwipi_ho_set_interpolation_function(self.name, ho_callback)
        current_cpl = ""


    def getName(self):
        """
        """
        return self.name


def ho_user_elt_set(elt_type, f_basis, f_loc):
    """
    """
    global interp_ho_loc_f
    global interp_ho_bas_f
    interp_ho_bas_f[0]=f_basis
    interp_ho_loc_f[0]=f_loc
    cwipi_ho_user_elt_set (<cwipi_element_t> elt_type,
                           ho_bas_callback,
                           ho_loc_callback)


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
                   const double local_coordinates[],
                   const int local_connectivity_index[],
                   const int local_connectivity[],
                   const int local_polyhedra_face_index[],
                   const int local_polyhedra_cell_to_face_connectivity[],
                   const int local_polyhedra_face_connectivity_index[],
                   const int local_polyhedra_face_connectivity[],
                   const double distant_points_coordinates[],
                   const int distant_points_location[],
                   const float distant_points_distance[],
                   const int distant_points_barycentric_coordinates_index[],
                   const double distant_points_barycentric_coordinates[],
                   int stride,
                   cwipi_solver_type_t  solver_type,
                   const void *local_field,
                   void *distant_field):
    """
    """
    global current_cpl
    global interp_f
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



cdef void ho_callback(int entities_dim,
                      int order,
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
                      int distant_points_weights_index[],
                      double distant_points_weights[],
                      double distant_points_uvw[],
                      int stride,
                      cwipi_solver_type_t  solver_type,
                      void *local_field,
                      void *distant_field):
    """
    """
    global current_cpl
    global interp_ho_f
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
        distant_points_weights_index_a = None
        distant_points_weights_a = None
        distant_points_uvw_a = None

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
        if (distant_points_weights_index == NULL or distant_points_weights == NULL):
            distant_points_weights_index_a = None
            distant_points_weights_a = None
        else:
            dims = <np.npy_intp>(n_distant_point + 1)
            distant_points_weights_index_a = np.PyArray_SimpleNewFromData(1,
                                                                          &dims,
                                                                          np.NPY_INT32,
                                                                          <void *>distant_points_weights_index)

            dims = <np.npy_intp>(distant_points_weights_index[n_distant_point])
            distant_points_weights_a = np.PyArray_SimpleNewFromData(1,
                                                                    &dims,
                                                                    np.NPY_DOUBLE,
                                                                    <void *>distant_points_weights)

        if (distant_points_uvw == NULL):
            distant_points_uvw_a = None
        else:
            dims = <np.npy_intp>(n_distant_point * entities_dim)
            distant_points_uvw_a = np.PyArray_SimpleNewFromData(1,
                                                                &dims,
                                                                np.NPY_DOUBLE,
                                                                <void *>distant_points_uvw)


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

    (<object>interp_ho_f[current_cpl])(entities_dim,
                                       order,
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
                                       distant_points_weights_index_a,
                                       distant_points_weights_a,
                                       distant_points_uvw_a,
                                       stride,
                                       solver_type,
                                       local_field_a,
                                       distant_field_a)



cdef void ho_bas_callback(int entities_dim,
                          int order,
                          int n_nodes,
                          int n_pts,
                          double *uvw,
                          double *weights):
    """
    """
    global interp_ho_bas_f
    np.import_array()

    cdef np.npy_intp dims = 0

    if (uvw == NULL):
        uvw_a = None
    else:
        dims = <np.npy_intp>(n_pts * entities_dim)
        uvw_a = np.PyArray_SimpleNewFromData(1,
                                             &dims,
                                             np.NPY_DOUBLE,
                                              <void *> uvw)

    if (weights == NULL):
        weights_a = None
    else:
        dims = <np.npy_intp>(n_pts * n_nodes)
        weights_a = np.PyArray_SimpleNewFromData(1,
                                                 &dims,
                                                 np.NPY_DOUBLE,
                                                 <void *> weights)

    (<object>interp_ho_bas_f[0])(entities_dim,
                                 order,
                                 n_nodes,
                                 n_pts,
                                 uvw_a,
                                 weights_a)



cdef double ho_loc_callback(int entities_dim,
                            int order,
                            int n_nodes,
                            double *nodes_coords,
                            double *point_coords,
                            double *projected_coords,
                            double *projected_uvw):
    """
    """
    global interp_ho_loc_f
    np.import_array()

    cdef np.npy_intp dims = 0

    if (nodes_coords == NULL):
        nodes_coords_a = None
    else:
        dims = <np.npy_intp>(n_nodes * 3)
        nodes_coords_a = np.PyArray_SimpleNewFromData(1,
                                                      &dims,
                                                      np.NPY_DOUBLE,
                                                      <void *> nodes_coords)

    if (point_coords == NULL):
        point_coords_a = None
    else:
        dims = <np.npy_intp>(3)
        point_coords_a = np.PyArray_SimpleNewFromData(1,
                                                      &dims,
                                                      np.NPY_DOUBLE,
                                                      <void *> point_coords)


    if (projected_coords == NULL):
        projected_coords_a = None
    else:
        dims = <np.npy_intp>(3)
        projected_coords_a = np.PyArray_SimpleNewFromData(1,
                                                          &dims,
                                                          np.NPY_DOUBLE,
                                                          <void *> projected_coords)



    if (projected_uvw == NULL):
        projected_uvw_a = None
    else:
        dims = <np.npy_intp>(entities_dim)
        projected_uvw_a = np.PyArray_SimpleNewFromData(1,
                                                          &dims,
                                                          np.NPY_DOUBLE,
                                                          <void *> projected_uvw)

    return (<object>interp_ho_loc_f[0])(entities_dim,
                                        order,
                                        n_nodes,
                                        nodes_coords_a,
                                        point_coords_a,
                                        projected_coords_a,
                                        projected_uvw_a)
