#!/usr/bin/env python
    
import mpi4py.MPI as MPI
import numpy as np
import sys

f=None

def userInterp(entities_dim,
                 n_local_vertex,
                 n_local_element,
                 n_local_polhyedra,
                 n_distant_point,
                 local_coordinates,
                 local_connectivity_index,
                 local_connectivity,
                 local_polyhedra_face_index,
                 local_polyhedra_cell_to_face,
                 local_polyhedra_face_connectivity_index,
                 local_polyhedra_face_connectivity,
                 distant_points_coordinates,
                 distant_points_location,
                 distant_points_barycentric_coordinates_index,
                 distant_points_barycentric_coordinates,
                 stride,
                 solver_type,
                 local_field,
                 distant_field):
    """
    User interpolation method
    """
    global f

    f.write('user interpolation')

    if (distant_field != None):
        for i in distant_field:
            i = 1.234


def runTest():
    """
    Run Python API test
    """
    global f
    commloc=MPI.Comm()
    commworld = MPI.COMM_WORLD

    rank = commworld.rank
    size = commworld.size

    if (rank == 0):
        print "\nSTART: python_api.py"

    if (size != 2):
        if rank == 0:
            print "      Not executed : only available for 2 processus"
        return

    applis = ["proc0","proc1"]
    fname = ["proc0","proc1"]

    try:
        import cwipi.CWIPI as CWIPI
    except:
        if rank == 0:
            print "      Error : CWIPI module not found (update PYTHONPATH variable)"
        sys.exit(1)

    #
    # Define a python file to write CWIPI outputs 

    if (rank == 0):
        print "        Output redirection"

    f=open(applis[rank]+".txt",'w')
    CWIPI.set_output_listing(f)
    CWIPI.init(MPI.COMM_WORLD, applis[rank], commloc)

    #
    # Control parameters

    if (rank == 0):
        print "        Control parameters"

    CWIPI.add_local_int_control_parameter("param_1", 1)
    CWIPI.synchronize_control_parameter(applis[(rank + 1) % 2])
    CWIPI.dump_application_properties()

    param_1 = CWIPI.get_local_int_control_parameter("param_1")
    param_2 = CWIPI.get_distant_int_control_parameter(applis[(rank + 1) % 2],"param_1")
    f.write('parametres: {param_1}, {param_2}'.format(param_1=param_1, param_2=param_2))

    #
    # Class coupling

    if (rank == 0):
        print "        Create coupling"

    # Constructor

    cpl = CWIPI.Coupling("cpl", 
                         CWIPI.COUPLING_PARALLEL_WITH_PARTITIONING, 
                         applis[(rank + 1) % 2] , 
                         2, 
                         0.1,  
                         CWIPI.STATIC_MESH, 
                         CWIPI.SOLVER_CELL_VERTEX, 
                         1, 
                         "Ensight", 
                         "txt")

    # Mesh

    if (rank == 0):
        print "        Create mesh"

    coord = np.array([-1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0], dtype=np.double)
    connec_idx = np.array([0, 4], dtype=np.int32)
    connec = np.array([1, 2, 3, 4], dtype=np.int32)

    cpl.define_mesh(4, 1, coord, connec_idx, connec)

    # Only send 

    if (rank == 0):
        print "        Exchange Proc 0 -> Proc 1"

    f.write("send :")

    sendField=np.array([0.1, 0.2, 0.3, 0.4], dtype=np.double)
    recvField=np.arange(4, dtype=np.double)

    if rank == 0:
        result = cpl.exchange("ech1", 
                              1, 
                              1, 
                              0.1, 
                              "field_s", sendField, 
                              "field_r", None)
    else:
        result = cpl.exchange("ech1", 
                              1, 
                              1, 
                              0.1, 
                              "field_s", None, 
                              "field_r", recvField)

    f.write('  - status : {param_1}'.format(param_1=result["status"]))
    if rank == 1:
        f.write("  - number of not located points : {param}".format(param=result["n_not_located_points"]))


    # Send and receive with user interpolation

    f.write("send receive :")

    if (rank == 0):
        print "        Define user interpolation"

    cpl.set_interpolation_function(userInterp) # User interpolation activation

    if (rank == 0):
        print "        Exchange Proc 0 <-> Proc 1 with user interpolation"

    result = cpl.exchange("ech2", 
                          1, 
                          1, 
                          0.1, 
                          "field_s2", sendField, 
                          "field_r2", recvField)

    f.write("  - status : {param_1}".format(param_1=result["status"]))
    f.write("  - number of not located points : {param}".format(param=result["n_not_located_points"]))

    # Delete coupling object

    if (rank == 0):
        print "        Delete coupling"

    del cpl

if __name__ == '__main__':
    runTest()

