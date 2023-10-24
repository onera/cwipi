---
jupytext:
  text_representation:
    extension: '.md'
    format_name: myst
    format_version: '0.7'
    jupytext_version: 1.4.0+dev
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Exercise 1 : a basic coupling

After having seen the core concepts to set up a coupling with CWIPI, we will discover the associated function calls in this very first basic coupling.

+++

*(Load custom magics)*

```{code-cell}
import os, sys
module_path = os.path.abspath(os.path.join('../../utils'))
if module_path not in sys.path:
    sys.path.append(module_path)
```

```{code-cell}
%reload_ext visu_magics
%reload_ext code_magics
```

+++

CWIPI has been written to function in a massivelly parallel distributed environement. Thus, the first thing to do, is the initialize the MPI environment:

```{code-cell}
%%code_block -p exercise_1 -i 1

#!/usr/bin/env python

import mpi4py.MPI as MPI

# Initialize MPI
comm = MPI.COMM_WORLD
i_rank = comm.rank
n_rank = comm.size
```
In the Python interface of CWIPI, all physical and geometric data (fields, meshes) are stored in Numpy arrays. Thus we import Numpy:

```{code-cell}
%%code_block -p exercise_1 -i 2

import numpy as np
```
The Python module of CWIPI since version 1.0 is called `pycwp`.
Let us import it and assure it has been found, that will come handy later when we will be executing this code:

```{code-cell}
%%code_block -p exercise_1 -i 3

# pycwp
try:
    from pycwp import pycwp
    if i_rank == 0:
        print("Yes, we have found the module :D")
except:
    if i_rank == 0:
        print("Oh no, couldn't find the module :'(")
```

### Initialization

Now we will start using, CWIPI functions.
Please refer to the API referenced [here](https://numerics.gitlab-pages.onera.net/coupling/cwipi/dev/index.html).

The function to start a CWIPI coupling between two codes is **init**. It takes the MPI communicator that includes the MPI ranks of all the coupled codes.
In this basic exemple, `code 1` will be running on the MPI rank 0 and `code 2` on the MPI rank 1. Thus, CWIPI will get the MPI communicator composed of MPI rank 0 and 1. Why do we provide the name of the solver as an array? Well, because since version 1.0 CWIPI allows several solvers to run on the same MPI rank. In this basic case, we only have one code per MPI rank. In real life applications the solvers run on more than one MPI rank. Since all MPI ranks calling the **init** function are supposed to take part in the CWIPI computations, it could come handy to force CWIPI not to use certain MPI ranks. That is what the argument is_active_rank is for. At initialization, CWIPI provides each solver the MPI communicators giving the processors the communicator to communicate through the ranks of that solver.
In our basic case, `code 1` gets a communicator with only MPI rank 0 and `code 2` get the communicator with only MPI rank 1.

```{code-cell}
%%code_block -p exercise_1 -i 4

n_code = 1
if (i_rank == 0):
    code_name = ["code1"]
if (i_rank == 1):
    code_name = ["code2"]
is_active_rank = True
intra_comm = pycwp.init(comm,
                        code_name,
                        is_active_rank)
```
### Coupling definition

Since a solver can take part in several couplings, the Coupling object creation allows to define the interaction between two fixed solvers. Let use a metafor to be more clear.

<span style="color:blue">*Oscar and Marie are two engineers and their boss assigned then to the CWIPI project to work in pairs. They don't know each other. During the first work session, they are each assigned to a desk in the working room. It is time to introduce themselves. Oscar is on the yellow desk and says "I am Oscar working on the CWIPI project with Marie. I am 28 years old and I live in Châtillon". Marie is on the blue desk and says "I am Marie working on the CWIPI project with Oscar. I am 54 years old and I live in Palaiseau".*</span>

In a similar way, at this step, we will introduce`code 1` and `code 2` to each other. On the MPI rank on which the solver is running, it will create a **Coupling** object telling which solver is running there, through which coupling it wants to communicate with which other solver. Then it describes itself in more detail.
First it provides the dimension of the coupling interface, if it is partitionned, the spatial interpolation algorithm it wants to use, the number of paritions on that MPI rank, if the coupling interface moves and that it is not an interpolation in time (temporal interpolation is not yet implemented in CWIPI).

```{code-cell}
%%code_block -p exercise_1 -i 5

if (i_rank == 0):
    coupled_code_name = ["code2"]
if (i_rank == 1):
    coupled_code_name = ["code1"]
n_part = 1
cpl = pycwp.Coupling(code_name[0],
                     "code1_code2",
                     coupled_code_name[0],
                     pycwp.INTERFACE_SURFACE,
                     pycwp.COMM_PAR_WITH_PART,
                     pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                     n_part,
                     pycwp.DYNAMIC_MESH_STATIC,
                     pycwp.TIME_EXCH_USER_CONTROLLED)
```

### Vizualisation

Let us take a pause in our coupling definition, to talk about the **visu_set** function. It allows to activate the Ensight ASCII output of the coupling interface with the exchanged fields and the partitionning. Those outputs can easily be read with Paraview.
When setting up a coupling, you will certainly have some tunning work to do. To be able to visualize the what CWIPI does will come handy to debug.

```{code-cell}
%%code_block -p exercise_1 -i 6

cpl.visu_set(1,
             pycwp.VISU_FORMAT_ENSIGHT,
             "text")
```

### Coupling interface

Let us go on with describing the coupling between `code 1` and `code 2`. What caracterizes the mesh we work on?

![alt text](mesh.png)

It is composed of several types of elements. To start of easy, let's just say it is composed of polygons. To be more precise 5 elements. We can also see 11 vertices on this mesh.

To define the coupling interface mesh in CWIPI, we first tell that we have vertices soup. It is just a set of coordinates of which we can make no sense. Then we create sense why telling CWIPI how to connect these vertices to form our polygons. Finally, CWIPI has to digest the information we provided it. Well, how does this translate in termes of code?

#### Set the mesh vertices coordinates

We start defining our vertices soup using the **mesh_interf_vtx_set** from the Coupling class. The coordinate system in CWIPI is always 3D, so we allocate an array of 3 times the number of vertices (11 here) to set the coordinates in. The coordinates are interlaced (x0, y0, z0, x1, y1, z1, ..., xn, yn, zn). The None argument will be explained later.

```{code-cell}
%%code_block -p exercise_1 -i 7

coords = np.array([0,0,0,  1,0,0,  2,0,0,  3,0,0,  0,1,0,  2,1,0, \
          3,1,0,  1,2,0,  0,3,0,  2,3,0,  3,3,0], dtype=np.double)
cpl.mesh_interf_vtx_set(0,
                        coords,
                        None)
```

#### Set the mesh polygons connectivity

Let us create sense in that vertices soup. The function **mesh_interf_block_add** allows us to tell that in that soup vertices are connected as polygons (CWP_BLOCK_FACE_POLY). Then we use the function **mesh_interf_f_poly_block_set** which allows to describe the 5 polygons of our 2D mesh. An index array (connec_idx) of size n_elts+1 contains the information of the number of vertices per polygon. The first index is always 0, from there we add up the number of vertices per element. Here one triangle, 2 quadrangles and 2 pentagons.
The connectivity between elements and vertices is an array of size connec_idx(n_elts+1) (here 21).

```{code-cell}
%%code_block -p exercise_1 -i 8

block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)
connec_idx = np.array([0,3,7,11,16,21], dtype=np.int32)
connec = np.array([1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8], dtype=np.int32)
cpl.mesh_interf_f_poly_block_set(0,
                                 block_id,
                                 connec_idx,
                                 connec,
                                 None)
```

#### Finalize mesh

This is when CWIPI digests the information we just provided it using the function **mesh_interf_finalize**. Indeed, CWIPI hides the parallelism for users but inside the code it needs to know the global numbering of the mesh entities. The None arguments given earlier allow the user to provide this global numbering.
If not given this numbering is generated by CWIPI, as well as the underlying mesh data structure

```{code-cell}
%%code_block -p exercise_1 -i 9

cpl.mesh_interf_finalize()
```

### Field definition

Now we know the mesh we work with. Let us define the fields of the solvers that are exchanged. As said earlier, here to simplify we will only send a field from `code 1` to `code 2`.

#### Create the field

The first step is to create a Field object attached to the Coupling object associated to the coupling between `code 1` and `code 2`. The numerical method of both solvers use vertex centered fields (DOF_LOCATION_NODE). For `code 1` we tell that this `super fancy field` will be send (FIELD_EXCH_SEND) and that `code 2` will receive it (FIELD_EXCH_RECV). In this basic coupling the `super fancy field` that will be send has only one component which is the x coordinate of the mesh coordinates. For each field we tell that we want to visualize it in the Ensight ASCII output (STATUS_ON).

```{code-cell}
%%code_block -p exercise_1 -i 10

n_components = 1
# for code1
if (i_rank == 0):
  field = cpl.field_create("a super fancy field",
                            pycwp.DOUBLE,
                            pycwp.FIELD_STORAGE_INTERLACED,
                            n_components,
                            pycwp.DOF_LOCATION_NODE,
                            pycwp.FIELD_EXCH_SEND,
                            pycwp.STATUS_ON)
# for code2
if (i_rank == 1):
  field = cpl.field_create("a super fancy field",
                           pycwp.DOUBLE,
                           pycwp.FIELD_STORAGE_INTERLACED,
                           n_components,
                           pycwp.DOF_LOCATION_NODE,
                           pycwp.FIELD_EXCH_RECV,
                           pycwp.STATUS_ON)
```

#### Set the field values

The function **data_set** of the Field class is used here to set the arrays associated to the fields. `code 1` fills an array with the data that it wants to send to `code 2`.
`code 2` has to provide an array in which the field data from `code 1` will be stored.

```{code-cell}
%%code_block -p exercise_1 -i 11

n_vtx = len(coords)//3
send_field_data = np.arange(n_vtx*n_components, dtype=np.double)
recv_field_data = np.arange(n_vtx*n_components, dtype=np.double)

for i in range(n_vtx):
  send_field_data[i] = coords[3*i]

# for code1
if (i_rank == 0):
  field.data_set(0,
                 pycwp.FIELD_MAP_SOURCE,
                 send_field_data)

# for code2
if(i_rank == 1):
  field.data_set(0,
                 pycwp.FIELD_MAP_TARGET,
                 recv_field_data)
```

### Begin time step

In this basic exemple, only one solver iteration during which an exchange occurs will be done. The begin and the end of an iteration have to be marked for CWIPI using **time_step_beg** and **time_step_end** function for each solver. This information allows CWIPI for instance to sort the visualization ouput of the fields per iteration.
Note, that is mandatory to create the coupling and the associated fields before starting the first time step.

```{code-cell}
%%code_block -p exercise_1 -i 12

pycwp.time_step_beg(code_name[0],
                    0.0);
```

### Compute interpolation weights

Since we use the spatial interpolation algorithm locating a set of points (vertices of `code 2`) in a mesh (coupling interface of `code 1`), to ensure all points are located a tolerence can be set using the function **spatial_interp_property_set** (optional).
Before doing any exchange, it is mandatory to compute the spatial interpolation weights using **spatial_interp_weights_compute**.

```{code-cell}
%%code_block -p exercise_1 -i 13

cpl.spatial_interp_property_set("tolerance",
                                pycwp.DOUBLE,
                                "0.1")

cpl.spatial_interp_weights_compute()
```

### Exchange field values between codes

For `code 1` to send its Field data array to `code 2`, the non-blocking **issend** should be called. Similarly, `code 2` should call **irecv** to tell `code 1` that is wants to receive the Field data array. After that, the solvers can do other work while the exchange is being done. Once you want to be sure the send operation has completed in `code 1`, use **wait_issend**.
The interpolated Field data array has completely arrived for `code 2` once the call to **wait_irecv** is completed.

```{code-cell}
%%code_block -p exercise_1 -i 14

# for code1
if (i_rank == 0):
  field.issend()

# for code2
if (i_rank == 1):
  field.irecv()

# for code1
if (i_rank == 0):
  field.wait_issend()

# for code2
if (i_rank == 1):
  field.wait_irecv()
```

### Check interpolation

As said earlier, one can set a tolerence to ensure all points are located. To check if that tolerence was large enougth, the function **n_uncomputed_tgts_get** can be called to retreive the number of unlocated vertices of the coupling interface of `code 2`.
To know which vertices were unlocated the **uncomputed_tgts_get** is called.

```{code-cell}
%%code_block -p exercise_1 -i 15

if (i_rank == 1):
  n_uncomputed_tgts = field.n_uncomputed_tgts_get(0);
  uncomputed_tgts   = field.uncomputed_tgts_get(0);
```

### End time step and clean up

At the end of each solver iteration **time_step_end** is called to inform CWIPI that the time step has terminated. When there are no CWIPI exchanges left to be done, all Field and Coupling objects can be deleted (that is done automatically by the garabage collector once there are no references left on it).
Still the coupling interface should be manually deleted calling **mesh_interf_del** on the Coupling object.

```{code-cell}
%%code_block -p exercise_1 -i 16

pycwp.time_step_end(code_name[0])

# Delete field

cpl.mesh_interf_del()

# Delete the coupling
```
### End CWIPI

This call terminates the use of CWIPI by cleaning up the internal structures CWIPI created.

```{code-cell}
%%code_block -p exercise_1 -i 17

pycwp.finalize()
```

### End MPI environment

At the end of the code the MPI environment should be terminated.

```{code-cell}
%%code_block -p exercise_1 -i 18

MPI.Finalize()
```

## Execution and visualization

Run the following cells to execute to program you just wrote and visualize the basic coupling you implemented.

```{code-cell}
%merge_code_blocks -l python -p exercise_1 -n 2 -v -c
```

```{code-cell}
%%visualize
visu/CODE1_CODE2.case
```
