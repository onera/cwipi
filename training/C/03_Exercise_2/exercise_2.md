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

# Exercise 2 : Coupling with a deformed mesh over time

Now that you know how to set up a basic coupling, let's go further by doing several coupling iterations.
At each iteration, the coupling interface mesh of `code 1` is deformed.

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
%reload_ext figure_magics
```

+++

## Initialization

Since the set up is roughly the same as in the previous exercise, it is not split up in as many code cells. Here we ask you to:
  - Initialize MPI
  - Initialize CWIPI for `code 1`
  - Set up the coupling : What value should be chosen for `CWP_Dynamic_mesh_t` since the coupling interface mesh of code 1 is deformed over time?
  - Ask for visualization outputs

```{code-cell}
%%code_block -p exercise_2_code_1 -i 1

#include "cwp.h"
#include "cwp_priv.h"

#include "pdm.h"
#include "pdm_generate_mesh.h"

int
main(int argc, char *argv[]) {

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // Initialize CWIPI :
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

  // Create the coupling :
  int n_part = 1;
  const char  *coupling_name     = "coupling";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  coupled_code_name[0] = "code2";
  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_DEFORMABLE,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  // Set coupling visualisation:
  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");
```

## What changes when the mesh is deformed over time?

Let's have a look again at the pseudo-code of the introduction.

```{prf:algorithm} basic coupling algorithm

**Inputs** Given $code1$ with a mesh $m1$ on which a field that will be sent is defined $sf1$. Given $code2$ with a mesh $m2$ on which a field that will be received is defined $rf2$.

**Output** $rf2$, which is $sf1$ interpolated on $m2$

1. Initialize CWIPI
2. Set coupling between $code1$ and $code2$
3. Describe codes:
   1. $code1$ has a mesh $m1$ on which we define a field $sf1$
   2. $code2$ has a mesh $m2$ and a receive buffer field $rf2$
4. Operate solver iterations:
   1. $code1$ sends $sf1$
   2. $code2$ receives $rf2$
5. Finalize CWIPI
```

*What would you change if `$m2$` was deformed at each `$code2$` iteration?
What about the field that we choose to be the x-coordinate of the coupling interface mesh?*

### Mesh

First we use a simple mesh generation function from ParaDiGM to create our coupling interface mesh : a square (nothing to do).

```{code-cell}
%%code_block -p exercise_2_code_1 -i 2

  // Create mesh :
  int     n_vtx = 0;
  int     n_elt = 0;
  double *coords      = NULL;
  int    *elt_vtx_idx = NULL;
  int    *elt_vtx     = NULL;
  PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[0]),
                                         10,
                                         &n_vtx,
                                         &n_elt,
                                         &coords,
                                         &elt_vtx_idx,
                                         &elt_vtx);
```

The mesh will change at each iteration but only its coordinates since it is but deformed.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 3

  CWP_Mesh_interf_vtx_set(code_name[0],
                              coupling_name,
                              0,
                              n_vtx,
                              coords,
                              NULL);

  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           coupling_name,
                                           CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   coupling_name,
                                   0,
                                   block_id,
                                   n_elt,
                                   elt_vtx_idx,
                                   elt_vtx,
                                   NULL);

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);
```

### Field

Actually, the mesh topology does not change. Thus, at each coupling iteration the number of vertices remains the same.
Thus, it suffices to provide the pointer to the field array and change the values inside it at each iteration.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 4

  const char *field_name      = "a super fancy field";
  int         n_components    = 1;

  CWP_Field_create(code_name[0],
                   coupling_name,
                   field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   CWP_DOF_LOCATION_NODE,
                   CWP_FIELD_EXCH_RECV,
                   CWP_STATUS_ON);

  double *field_data = malloc(sizeof(double) * n_vtx);

  CWP_Field_data_set(code_name[0],
                     coupling_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_TARGET,
                     field_data);
```

In the case the topology of your mesh changes at each iteration, the new field array will be set at each iteration.
It is important to know that the field should still be created before starting the first time step.

## Time iterations

At the beginning of each coupling iteration, we begin a new time step using `CWP_Time_step_beg` which we will terminate
at the end of the iteration with `CWP_Time_step_end`.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 6

  const int    itdeb = 1;
  const int    itend = 10;
  double       ttime = 0.0;
  double       dt    = 0.1;

  double degrad = acos(-1.0)/180.;
  double x = 0.0;
  double y = 0.0;
  double alpha = 2;
  alpha = alpha * degrad;
  double sina = sin(alpha);
  double cosa = cos(alpha);

  for (int it = itdeb; it <= itend; it ++) {

    ttime = (it-itdeb)*dt;

    // Start time step
    CWP_Time_step_beg(code_name[0],
                      ttime);
```

Let's rotate the mesh of `code 1` with respect to `code 2`.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 7

    if (it > itdeb) {
      for (int i = 0; i < n_vtx; i++) {
        x = coords[i * 3    ];
        y = coords[i * 3 + 1];
        coords[i * 3    ] = cosa*x - sina*y;
        coords[i * 3 + 1] = sina*x + cosa*y;
        field_data[i] = coords[3 * i];
      }
    }

```

The aim is to interpolate the field of `code 1` onto the mesh of `code 2`.
*What happens to the interpolation weights when the mesh of `code 1` is deformed?
Thus, what does that induce in your code?*

The chosen tolerance does not change here over time, so we set it before the iteration loop.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 5

  CWP_Spatial_interp_property_set(code_name[0],
                                  coupling_name,
                                  "tolerance",
                                  CWP_DOUBLE,
                                  "0.001");
```

But the weights need to be computed at each iteration after the mesh has been deformed, so that is done in the iteration loop.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 8

    CWP_Spatial_interp_weights_compute(code_name[0],
                                       coupling_name);

```

Now we receive the field send by `code 2`.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 9

    CWP_Field_irecv(code_name[0],
                    coupling_name,
                    field_name);

    CWP_Field_wait_irecv(code_name[0],
                         coupling_name,
                         field_name);
```

Earlier we set a tolerence for the localization algorithm.
To check if that tolerence was large enougth, the function **CWP_N_uncomputed_tgts_get** can be called to retreive the number of unlocated vertices of the coupling interface of `code 1`.
To know which vertices were unlocated the **CWP_Uncomputed_tgts_get** is called.

```{code-cell}
%%code_block -p exercise_1_code_1 -i 13

  int        n_uncomputed_tgts = -1;
  const int *uncomputed_tgts   = NULL;
  n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_name[0],
                                                coupling_name,
                                                field_name,
                                                0);

  uncomputed_tgts = CWP_Uncomputed_tgts_get(code_name[0],
                                            coupling_name,
                                            field_name,
                                            0);
```

Let's have a sneak peek in this algorithm through this animation which will help you understand what we mean by unlocated points.

```{code-cell}
%%localization
unlocated
```

*Spoiler : At the end of the exercise you will see that since the coupling interface mesh of `code 1` moves
there are unlocated points with the tolerance set to 0.001. Increasing it will eventually let all points be located
but at the cost of the time taken by the algorithm. You call play around with the tolerance once you finish the exercise.*

Let's end the iteration.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 10

    CWP_Time_step_end(code_name[0]);

  } // end iterations
```

## Finalize

Let us finish the coupling by freeing the memory allocated for it and ending this program.

```{code-cell}
%%code_block -p exercise_2_code_1 -i 11

  // Delete field :
  CWP_Field_del(code_name[0],
                coupling_name,
                field_name);

  // Delete Mesh :
  CWP_Mesh_interf_del(code_name[0],
                      coupling_name);

  // Delete the coupling :
  CWP_Cpl_del(code_name[0],
              coupling_name);

  // free
  free(intra_comm);
  free(code_name);
  free(coupled_code_name);
  free(coords);
  free(elt_vtx_idx);
  free(elt_vtx);
  free(field_data);

  // Finalize CWIPI :
  CWP_Finalize();

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}

```

## Execution and visualization

Run the following cells to execute to program you just wrote and visualize the basic coupling you implemented.

```{code-cell}
%merge_code_blocks -l c -p exercise_2_code_1 -n 1 -v -c
```

```{code-cell}
%%visualize_dynamic
cwipi_writer/coupling_code1_code2/CHR.case : r_a~super~fancy~field1
cwipi_writer/coupling_code2_code1/CHR.case : s_a~super~fancy~field1
```



