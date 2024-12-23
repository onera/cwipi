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

+++ {"editable": false, "deletable": false}

# Exercise 2 : Coupling with a deformed mesh over time

+++ {"editable": false, "deletable": false}

Now that you know how to set up a basic coupling, let's go further by doing several coupling iterations.
At each iteration, the coupling interface mesh of `code1` is deformed.

+++ {"editable": false, "deletable": false}

## Load magic commands

As usual we start by loading the custom magic commands.

```{code-cell}
---
"editable": false
"deletable": false
---
import os, sys
module_path = os.path.abspath(os.path.join('../../utils'))
if module_path not in sys.path:
    sys.path.append(module_path)
```

```{code-cell}
---
"editable": false
"deletable": false
---
%reload_ext visu_magics
%reload_ext code_magics
%reload_ext figure_magics
```

+++ {"editable": false, "deletable": false}

## Initialization

Since the set up is roughly the same as in the previous exercise, it is not split up in as many code cells. Here we ask you to:
  - Initialize MPI
  - Initialize CWIPI for `code1`
  - Set up the coupling : What value should be chosen for `CWP_Dynamic_mesh_t` since the coupling interface mesh of `code1` is deformed over time?
  - Ask for visualization outputs

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 1

program fortran_new_api_deformable_sol

    use cwp

    implicit none

    include "mpif.h"

  !--------------------------------------------------------------------
  integer                                     :: ierr
  integer                                     :: i_rank, n_rank
  integer(kind=8), parameter                  :: n_vtx_seg = 10

  integer                                     :: n_code
  character(len = 99),                pointer :: code_names(:) => null()
  integer                                     :: is_active_rank = CWP_STATUS_ON
  integer,                            pointer :: intra_comms(:) => null()

  integer                                     :: n_part
  character(len = 99),                pointer :: coupled_code_names(:) => null()
  character(len = 99)                         :: coupling_name

  double precision,   dimension(:,:), pointer :: coords    => null()
  integer(c_long),    dimension(:),   pointer :: vtx_g_num => null()

  integer,            dimension(:),   pointer :: elt_vtx_idx => null()
  integer,            dimension(:),   pointer :: elt_vtx     => null()
  integer(c_long),    dimension(:),   pointer :: elt_g_num   => null()
  integer(c_int)                              :: id_block

  character(len = 99)                         :: field_name
  integer(c_int)                              :: n_components

  double precision                            :: dt
  double precision                            :: degrad
  double precision                            :: x
  double precision                            :: y
  double precision                            :: alpha
  double precision                            :: sina
  double precision                            :: cosa
  double precision                            :: time
  integer                                     :: i, it, itdeb, itend
  integer                                     :: n_vtx, n_elt

  double precision,                   pointer :: field_data(:) => null()

  integer(c_int)                              :: n_uncomputed_tgts
  integer(c_int),                     pointer :: uncomputed_tgts(:) => null()
  !--------------------------------------------------------------------

  ! MPI Initialization :
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! Initialize CWIPI :
  n_code = 1

  allocate(code_names(n_code), &
           intra_comms(n_code))

  code_names(1) = "code1"

  call CWP_Init(mpi_comm_world, &
                n_code,         &
                code_names,     &
                is_active_rank, &
                intra_comms)

  ! Create the coupling :
  ! CWP_DYNAMIC_MESH_DEFORMABLE allows us to take into account the modifications
  ! to the mesh over the coupling steps.
  coupling_name     = "coupling"
  allocate(coupled_code_names(n_code))
  coupled_code_names(1) = "code2"
  n_part = 1;
  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_DEFORMABLE,                           &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  ! Set coupling visualisation:
  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")
```

+++ {"editable": false, "deletable": false}

## What changes when the mesh is deformed over time?

Let's have a look again at the pseudo-code of the introduction.

<!-- ```{prf:algorithm} basic coupling algorithm

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
``` -->
**Inputs** Given `code1` with a mesh `m1` on which a field that will be sent is defined `sf1`. Given `code2` with a mesh `m2` on which a field that will be received is defined `rf2`.

**Output** `rf2`, which is `sf1` interpolated on `m2`

1. Initialize CWIPI
2. Set coupling between `code1` and `code2`
3. Describe codes:
   1. `code1` has a mesh `m1` on which we define a field `sf1`
   2. `code2` has a mesh `m2` and a receive buffer field `rf2`
4. Operate solver iterations:
   1. `code1` sends `sf1`
   2. `code2` receives `rf2`
5. Finalize CWIPI

In this exercise `$code1$` receives a field send by `$code2$`.
Here we decide to deform `$m1$` at each `$code1$` iteration.

*What does that change in our coupling code?
What would happen if `$code1$` would send `sf1`?*

### Mesh

First we use a simple mesh generation function from ParaDiGM to create our coupling interface mesh : a square (nothing to do).
It is composed of triangle elements (i.e. `CWP_BLOCK_FACE_TRIA3`).

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 2

  ! Create mesh :
  call CWPT_generate_mesh_rectangle_simplified(intra_comms(1), &
                                               n_vtx_seg,      &
                                               n_vtx,          &
                                               n_elt,          &
                                               coords,         &
                                               elt_vtx_idx,    &
                                               elt_vtx)
```

+++ {"editable": false, "deletable": false}

The mesh will change at each iteration. Since it is deformed, only its coordinates change.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 3

  call CWP_Mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               coords,        &
                               vtx_g_num)

  id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                       coupling_name,       &
                                       CWP_BLOCK_FACE_POLY)

  call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                        coupling_name, &
                                        0,             &
                                        id_block,      &
                                        n_elt,         &
                                        elt_vtx_idx,   &
                                        elt_vtx,       &
                                        elt_g_num)

  call CWP_Mesh_interf_finalize(code_names(1), &
                                coupling_name)
```

+++ {"editable": false, "deletable": false}

### Field

In this exercise, `code1` receives a field from `code2`.
The mesh changes at each iteration but here we decided that the field `code2` sends wouldn't.
In a real case application, `code2` would send a different field a each iteration.
Since the mesh topology does not change, the coupling code would be the same since a pointer on the field is provided and there are no internal copies of the provided fields inside CWIPI.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 4

  field_name   = "a super fancy field"
  n_components = 1

  call CWP_Field_create(code_names(1),                &
                        coupling_name,                &
                        field_name,                   &
                        CWP_DOUBLE,                   &
                        CWP_FIELD_STORAGE_INTERLACED, &
                        n_components,                 &
                        CWP_DOF_LOCATION_NODE,        &
                        CWP_FIELD_EXCH_RECV,          &
                        CWP_STATUS_ON)

  allocate(field_data(n_vtx*n_components))

  call CWP_Field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_TARGET, &
                          field_data)
```

+++ {"editable": false, "deletable": false}

In the case the topology of your mesh changes at each iteration, the new field array will be set at each iteration.
It is important to know that the field should still be created before starting the first time step.

## Time iterations

At the beginning of each coupling iteration, we begin a new time step using `CWP_Time_step_beg` which we will terminate
at the end of the iteration with `CWP_Time_step_end`.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 6

  ! Interations :
  ! At each iteration the mesh coordinates and the exchanged fields are modified.
  itdeb = 1
  itend = 10
  time  = 0.0d0
  dt    = 0.1d0

  degrad = acos(-1.0)/180.
  x = 0.0
  y = 0.0
  alpha = 2
  alpha = alpha * degrad
  sina = sin(alpha)
  cosa = cos(alpha)

  do it = itdeb, itend

    time = (it-itdeb)*dt

    ! Begin time step :
    call CWP_Time_step_beg(code_names(1), &
                           time)
```

+++ {"editable": false, "deletable": false}

Let's rotate the mesh of `code1` with respect to `code2`.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 7

    if (it > itdeb) then
      do i = 1, n_vtx
        x = coords(1, i)
        y = coords(2, i)
        coords(1, i) = cosa*x - sina*y
        coords(2, i) = sina*x + cosa*y
        field_data(i) = coords(1, i)
      enddo
    endif

```

+++ {"editable": false, "deletable": false}

The aim is to interpolate the field of `code2` onto the mesh of `code1`.
*What happens to the interpolation weights when the mesh of `code1` is deformed?
Thus, what does that induce in your code?*

The chosen tolerance does not change here over time, so we set it before the iteration loop.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 5

  call CWP_Spatial_interp_property_set(code_names(1), &
                                       coupling_name, &
                                       "tolerance",   &
                                       CWP_DOUBLE,    &
                                       "0.001")
```

+++ {"editable": false, "deletable": false}

But the weights need to be computed at each iteration after the mesh has been deformed, so that is done in the iteration loop.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 8

    call CWP_Spatial_interp_weights_compute(code_names(1), &
                                            coupling_name)

```

+++ {"editable": false, "deletable": false}

Now we receive the field send by `code2`.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 9

    call CWP_Field_irecv(code_names(1), &
                         coupling_name, &
                         field_name)

    call CWP_Field_wait_irecv(code_names(1), &
                              coupling_name, &
                              field_name)
```

+++ {"editable": false, "deletable": false}

Earlier we set a tolerance for the localization algorithm.
To check if that tolerance was large enough, the function **CWP_N_uncomputed_tgts_get** can be called to retrieve the number of unlocated vertices of the coupling interface of `code1`.
To know which vertices were unlocated the **CWP_Uncomputed_tgts_get** is called.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1_code_1 -i 13

  n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_names(1),   &
                                                coupling_name,   &
                                                field_name, &
                                                0)

  if (n_uncomputed_tgts /= 0) then
    uncomputed_tgts => CWP_Uncomputed_tgts_get(code_names(1),   &
                                               coupling_name,   &
                                               field_name, &
                                               0)
  endif
```

+++ {"editable": false, "deletable": false}

Let's have a sneak peek in this algorithm through this animation which will help you understand what we mean by unlocated points.

```{code-cell}
---
"editable": false
"deletable": false
---
%%localization
unlocated
```

+++ {"editable": false, "deletable": false}

*Spoiler : At the end of the exercise you will see that since the coupling interface mesh of `code1` moves
there are unlocated points with the tolerance set to 0.001. Increasing it will eventually let all points be located
but at the cost of the time taken by the algorithm. You call play around with the tolerance once you finish the exercise.*

Let's end the iteration.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 10

    call CWP_Time_step_end(code_names(1))

  enddo
```

+++ {"editable": false, "deletable": false}

## Finalize

Let us finish the coupling by freeing the memory allocated for it and ending this program.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 11

  ! Delete field :
  call CWP_Field_Del(code_names(1),   &
                     coupling_name,   &
                     field_name)

  ! Delete Mesh :
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  ! Delete the coupling :
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)

  ! free
  deallocate(code_names)
  deallocate(intra_comms)
  deallocate(coupled_code_names)
  deallocate(field_data)
  deallocate(uncomputed_tgts)

  call cwpt_fortran_free_c(c_loc(coords))
  call cwpt_fortran_free_c(c_loc(elt_vtx_idx))
  call cwpt_fortran_free_c(c_loc(elt_vtx))

  ! Finalize CWIPI :
  call CWP_Finalize()

  ! Finalize MPI :
  call MPI_Finalize(ierr)

end program fortran_new_api_deformable_sol

```

+++ {"editable": false, "deletable": false}

## Execution and visualization

Run the following cells to execute to program you just wrote and visualize the basic coupling you implemented.

```{code-cell}
---
"deletable": false
---
%merge_code_blocks -l fortran -p exercise_2_code_1 -n 1 -v -c
```

```{code-cell}
---
"deletable": false
---
%%visualize_dynamic
cwipi_writer/coupling_code1_code2/CHR.case : r_a~super~fancy~field1
cwipi_writer/coupling_code2_code1/CHR.case : s_a~super~fancy~field1
```

+++ {"editable": false, "deletable": false}

<span style="color:red">*You reached the end of this training. Congratulations ! Feel free to give us feedback.*</span>
