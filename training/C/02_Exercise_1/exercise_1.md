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
To help you with this, you are encouraged to look at the [documentation](https://numerics.gitlab-pages.onera.net/coupling/cwipi/dev/new_cwipi/new_cwipi.html#c-api-documentation).

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

CWIPI has been written to function in a massivelly parallel distributed environement.
Thus, the first thing to do, is the initialize the MPI environment:

```{code-cell}
%%code_block -p exercise_1_code_1 -i 1

#include "cwp.h"
#include "cwp_priv.h"

int
main(int argc, char *argv[]) {

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);
```

### Initialization

Now we will start using CWIPI functions !
Please refer to the API referenced [here](https://numerics.gitlab-pages.onera.net/coupling/cwipi/dev/index.html).

The function to start a CWIPI coupling between two codes is **CWP_Init**. It takes the MPI communicator that includes the MPI ranks of all the coupled codes.
In this basic exemple, `code 1` will be running on the MPI rank 0 and `code 2` on the MPI rank 1.
Thus, CWIPI will get the MPI communicator composed of MPI rank 0 and 1. Why do we provide the name of the solver as an array?
Well, because since version 1.0 CWIPI allows several solvers to run on the same MPI rank.
In this basic case, we only have one code per MPI rank. In real life applications the solvers run on more than one MPI rank.
Since all MPI ranks calling the **CWP_Init** function are supposed to take part in the CWIPI computations, it could come handy
to force CWIPI not to use certain MPI ranks. That is what the argument is_active_rank is for.
At initialization, CWIPI provides each solver the MPI communicators giving the processors the communicator to communicate through
the ranks of that solver.
In our basic case, `code 1` gets a communicator with only MPI rank 0 and `code 2` get the communicator with only MPI rank 1.

*Remark : In this exercise you will be doing the CWIPI calls only for `code 1`. We already implemented the calls for `code 2`
in `exercise_1_code2.c` in this folder. There is no point in cheating, you are here to learn.*

```{code-cell}
%%code_block -p exercise_1_code_1 -i 2

  int n_code = 1;
  const char  **code_name      = malloc(sizeof(char *) * n_code);
  CWP_Status_t  is_active_rank = CWP_STATUS_ON;
  MPI_Comm     *intra_comm     = malloc(sizeof(MPI_Comm) * n_code);

  code_name[0] = "code1";

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);
```
### Coupling definition

Since a solver can take part in several couplings, the Coupling object creation allows to define the interaction between two fixed solvers. Let use a metafor to be more clear.

<span style="color:blue">*Oscar and Marie are two engineers and their boss assigned then to the CWIPI project to work in pairs. They don't know each other. During the first work session, they are each assigned to a desk in the working room. It is time to introduce themselves. Oscar is on the yellow desk and says "I am Oscar working on the CWIPI project with Marie. I am 28 years old and I live in Châtillon". Marie is on the blue desk and says "I am Marie working on the CWIPI project with Oscar. I am 54 years old and I live in Palaiseau".*</span>

In a similar way, at this step, we will introduce`code 1` and `code 2` to each other. On the MPI rank on which the solver is running, it will create a coupling structure telling which solver is running there, through which coupling it wants to communicate with which other solver. Then it describes itself in more detail.
First it provides the dimension of the coupling interface, if it is partitionned, the spatial interpolation algorithm it wants to use, the number of paritions on that MPI rank, if the coupling interface moves and that it is not an interpolation in time (temporal interpolation is not yet implemented in CWIPI).

```{code-cell}
%%code_block -p exercise_1_code_1 -i 3

  int n_part = 1;
  const char  *coupling_name     = "code1_code2";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);

  coupled_code_name[0] = "code2";

  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);
```

### Vizualisation

Let us take a pause in our coupling definition, to talk about the **CWP_Visu_set** function.
It allows to activate the Ensight ASCII output of the coupling interface with the exchanged fields and the partitionning.
Those outputs can easily be read with Paraview.
When setting up a coupling, you will certainly have some tunning work to do. To be able to visualize the what CWIPI does will come handy to debug.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 4

  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");
```

### Coupling interface

Let us go on with a description of the coupling between `code 1` and `code 2`. What caracterizes the mesh we work on?

![alt text](mesh_code1.png)

It is a basic cartesian grid mesh composed of 9 squares and 16 vertices.

The coupling interface mesh of `code 2` looks like this.

![alt text](mesh_code1.png)

It is composed of several types of elements.
To start of easy, let's just say it is composed of polygons.
To be more precise 5 elements.
We can also see 11 vertices on this mesh.

We would like to emphasise that the meshes do not have to be coincident in order to couple using CWIPI.

To define the coupling interface mesh in CWIPI, we first tell that we have vertices soup.
It is just a set of coordinates of which we can make no sense. Then we create sense why telling CWIPI how to connect these vertices to form our polygons.
Finally, CWIPI has to digest the information we provided it. Well, how does this translate in terms of code?

#### Set the mesh vertices coordinates

We start defining our vertices soup using the **CWP_Mesh_interf_vtx_set** from the Coupling class.
The coordinate system in CWIPI is always 3D, so we allocate an array of 3 times the number of vertices (16 here) to set the coordinates in.
The coordinates are interlaced (x0, y0, z0, x1, y1, z1, ..., xn, yn, zn). The None argument will be explained later.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 5

  int    n_vtx = 16;
  double coords[48] = {0,0,0,  1,0,0,  2,0,0,  3,0,0,
                       0,1,0,  1,1,0,  2,1,0,  3,1,0,
                       0,2,0,  1,2,0,  2,2,0,  3,2,0,
                       0,3,0,  1,3,0,  2,3,0,  3,3,0};
  CWP_Mesh_interf_vtx_set(code_name[0],
                          coupling_name,
                          0,
                          n_vtx,
                          coords,
                          NULL);
```

#### Set the mesh polygons connectivity

Let us create sense in that vertices soup. The function **CWP_Mesh_interf_block_add** allows us to tell that in that soup vertices are connected as polygons (CWP_BLOCK_FACE_POLY).
Then we use the function **CWP_Mesh_interf_f_poly_block_set** which allows to describe the 9 polygons of our 2D mesh. An index array (connec_idx) of size n_elts+1 contains the information of the number of vertices per polygon.
The first index is always 0, from there we add up the number of vertices per element. Here the mesh is composed only of squares (4 vertices).
The connectivity between elements and vertices is an array of size connec_idx(n_elts+1) (here 36).

```{code-cell}
%%code_block -p exercise_1_code_1 -i 6

  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           coupling_name,
                                           CWP_BLOCK_FACE_POLY);

  int n_elts = 9;
  int connec_idx[10] = {0,4,8,12,16,20,24,28,32,36};
  int connec[36]     = {1,2,6,5,     2,3,7,6,      3,4,8,7,
                        5,6,10,9,    6,7,11,10,    7,8,12,11,
                        9,10,14,13,  10,11,15,14,  11,12,16,15};
  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   coupling_name,
                                   0,
                                   block_id,
                                   n_elts,
                                   connec_idx,
                                   connec,
                                   NULL);
```

#### Finalize mesh

This is when CWIPI digests the information we just provided it using the function **CWP_Mesh_interf_finalize**. Indeed, CWIPI hides the parallelism for users but inside the code it needs to know the global numbering of the mesh entities. The None arguments given earlier allow the user to provide this global numbering.
If not given this numbering is generated by CWIPI, as well as the underlying mesh data structure

```{code-cell}
%%code_block -p exercise_1_code_1 -i 7

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);
```

### Field definition

Now we know the mesh we work with. Let us define the fields of the solvers that are exchanged. As said earlier, here to simplify we will only send a field from `code 1` to `code 2`.

#### Create the field

The first step is to create a Field object attached to the Coupling object associated to the coupling between `code 1` and `code 2`. The numerical method of both solvers use vertex centered fields (DOF_LOCATION_NODE). For `code 1` we tell that this `super fancy field` will be send (FIELD_EXCH_SEND) and that `code 2` will receive it (FIELD_EXCH_RECV). In this basic coupling the `super fancy field` that will be send has only one component which is the x coordinate of the mesh coordinates. For each field we tell that we want to visualize it in the Ensight ASCII output (STATUS_ON).

```{code-cell}
%%code_block -p exercise_1_code_1 -i 8

  const char *field_name      = "a super fancy field";
  int         n_components    = 1;

  CWP_Field_create(code_name[0],
                   coupling_name,
                   field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   CWP_DOF_LOCATION_NODE,
                   CWP_FIELD_EXCH_SEND,
                   CWP_STATUS_ON);
```

#### Set the field values

The function **CWP_Field_data_set** of the Field class is used here to set the arrays associated to the fields. `code 1` fills an array with the data that it wants to send to `code 2`.
`code 2` has to provide an array in which the field data from `code 1` will be stored.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 9

  double *send_field_data = malloc(sizeof(double) * n_vtx * n_components);

  for (int i = 0; i < n_vtx; i++) {
    send_field_data[i] = coords[3*i];
  }
  CWP_Field_data_set(code_name[0],
                     coupling_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_SOURCE,
                     send_field_data);
```

### Begin time step

In this basic exemple, only one solver iteration during which an exchange occurs will be done.
The begin and the end of an iteration have to be marked for CWIPI using **CWP_Time_step_beg** and **CWP_Time_step_end** function for each solver. This information allows CWIPI for instance to sort the visualization ouput of the fields per iteration.
Note, that is mandatory to create the coupling and the associated fields before starting the first time step.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 10

  CWP_Time_step_beg(code_name[0],
                    0.0);
```

### Compute interpolation weights

Since we use the spatial interpolation algorithm locating a set of points (vertices of `code 2`) in a mesh (coupling interface of `code 1`), to ensure all points are located a tolerence can be set using the function **CWP_Spatial_interp_property_set** (optional).
Before doing any exchange, it is mandatory to compute the spatial interpolation weights using **CWP_Spatial_interp_weights_compute**.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 11

  CWP_Spatial_interp_property_set(code_name[0],
                                  coupling_name,
                                  "tolerance",
                                  CWP_DOUBLE,
                                  "0.001");
  CWP_Spatial_interp_weights_compute(code_name[0],
                                     coupling_name);
```

### Exchange field values between codes

For `code 1` to send its Field data array to `code 2`, the non-blocking **CWP_Field_issend** should be called. Similarly, `code 2` should call **CWP_Field_irecv** to tell `code 1` that is wants to receive the Field data array. After that, the solvers can do other work while the exchange is being done. Once you want to be sure the send operation has completed in `code 1`, use **CWP_Field_wait_issend**.
The interpolated Field data array has completely arrived for `code 2` once the call to **CWP_Field_wait_irecv** is completed.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 12

  CWP_Field_issend(code_name[0],
                   coupling_name,
                   field_name);

  CWP_Field_wait_issend(code_name[0],
                        coupling_name,
                        field_name);
```

### End time step and clean up

At the end of each solver iteration **CWP_Time_step_end** is called to inform CWIPI that the time step has terminated.
When there are no CWIPI exchanges left to be done, all field and coupling structures can be deleted (**CWP_Field_del** and **CWP_Cpl_del**).
Still the coupling interface should be manually deleted calling **CWP_Mesh_interf_del** on the Coupling object.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 14

  CWP_Time_step_end(code_name[0]);

  CWP_Field_del(code_name[0],
                coupling_name,
                field_name);

  CWP_Mesh_interf_del(code_name[0],
                      coupling_name);

  CWP_Cpl_del(code_name[0],
              coupling_name);
```

### End CWIPI

This call terminates the use of CWIPI by cleaning up the internal structures CWIPI created.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 15

  free(send_field_data);
  free(code_name);
  free(intra_comm);
  free(coupled_code_name);

  CWP_Finalize();

```

### End MPI environment

At the end of the code the MPI environment should be terminated.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 16

  MPI_Finalize();

  return EXIT_SUCCESS;
}
```

## Execution and visualization

Run the following cells to execute to program you just wrote and visualize the basic coupling you implemented.

```{code-cell}
%merge_code_blocks -l c -p exercise_1_code_1 -n 1 -v -c
```

```{code-cell}
%%visualize
cwipi_writer/code1_code2_code1_code2/CHR.case : s_a~super~fancy~field1
cwipi_writer/code1_code2_code2_code1/CHR.case : r_a~super~fancy~field1
```

# Bonus : a coupling with conservative interpolation

If you read this, it means that you quickly finished the first exercise. Congratulations !
As you have seen in the introduction, from version 1.x on, CWIPI has several spatial interpolation algorithms.
To go further, we invite you to repeat the exercise above but with a conservative interpolation algorithm (CWP_SPATIAL_INTERP_FROM_INTERSECTION).
Copy paste the following code to have a field defined on the faces. Adapt the code accordingly to those two changes. Observe the output.

```{prf:algorithm} basic couling algorithm

double *send_field_data = malloc(sizeof(double) * n_elt * n_components);

for (int i = 0; i < n_elt; i++) {
  send_field_data[i] = i;
}
```
