!
! File:   mod_pdm_mesh_nodal.F90
! Author: equemera
!
! Created on July 10, 2017, 1:34 PM
!

MODULE pdm_mesh_nodal

  use pdm
  use iso_c_binding

  implicit none

  contains


!> Create a Mesh nodal structure
!!
!! @param[in]   n_part   Number of partition on the current process
!! @param[out]  mesh     Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine pdm_mesh_nodal_create (n_part, mesh)
    use iso_c_binding

    implicit none

    integer, intent(in)  :: n_part
    type(c_ptr)          :: mesh

    interface
      function pdm_mesh_nodal_create_c (n_part) result(mesh) bind(c, name='PDM_Mesh_nodal_create')
        use iso_c_binding
        implicit none
        integer(c_int), intent (in), value :: n_part
        type(c_ptr)                        :: mesh
      end function pdm_mesh_nodal_create_c
    end interface

    mesh = pdm_mesh_nodal_create_c (n_part)

  end subroutine pdm_mesh_nodal_create


!> \brief Free partially a nodal mesh structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine pdm_mesh_nodal_partial_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value     :: mesh

    interface
      subroutine pdm_mesh_nodal_partial_free_c (mesh) &
        bind(c, name='PDM_Mesh_nodal_partial_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value     :: mesh
      end subroutine pdm_mesh_nodal_partial_free_c
    end interface

    call pdm_mesh_nodal_partial_free_c (mesh)

  end subroutine pdm_mesh_nodal_partial_free


!> Free a nodal mesh structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine pdm_mesh_nodal_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value :: mesh

    interface
      subroutine pdm_mesh_nodal_free_c (mesh) &
        bind(c, name='PDM_Mesh_nodal_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value :: mesh
      end subroutine pdm_mesh_nodal_free_c
    end interface

    call pdm_mesh_nodal_free_c (mesh)

  end subroutine pdm_mesh_nodal_free


!> Define partition vertices
!!
!! @param[in]  mesh     Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part  Partition identifier
!! @param[in]  n_vtx    Number of vertices
!! @param[in]  coords   Interlaced coordinates (size = 3 * \ref n_vtx)
!! @param[in]  numabs   Global numbering
!!

  subroutine pdm_mesh_nodal_coord_set (mesh, id_part, n_vtx, coords, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_vtx
    integer, intent(in)                 :: owner
    double precision, pointer           :: coords(:,:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_coords = C_NULL_PTR
    type(c_ptr) :: c_numabs = C_NULL_PTR

    interface
      subroutine pdm_mesh_nodal_coord_set_c(mesh, id_part, n_vtx, coords, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_coord_set')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_vtx
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: coords
        type (c_ptr),                value :: numabs
      end subroutine pdm_mesh_nodal_coord_set_c
    end interface

    c_coords = c_loc (coords)
    c_numabs = c_loc (numabs)

    call  pdm_mesh_nodal_coord_set_c(mesh, id_part, n_vtx, c_coords, c_numabs, owner)

  end subroutine

END MODULE pdm_mesh_nodal
