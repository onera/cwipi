!-----------------------------------------------------------------------------
! This file is part of the ParaDiGMA library.
!
! Copyright (C) 2020  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE  
  use mpi
#endif  
  use pdm_overlay
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  integer :: code
  integer :: i_rank
  integer :: n_rank

  integer :: id

  integer (c_int), parameter ::  n_partMeshA = 1
  integer (kind = pdm_g_num_s), parameter ::  nGFaceMeshA = 1
  integer (kind = pdm_g_num_s), parameter ::  nGVtxMeshA = 3
  integer (c_int), parameter ::  n_partMeshB = 1
  integer (kind = pdm_g_num_s), parameter ::  nGFaceMeshB = 1
  integer (kind = pdm_g_num_s), parameter ::  nGVtxMeshB = 3
  double precision, parameter :: projectCoeff = 0.d0


  integer (c_int), parameter              :: ipartA = 0
  integer (c_int), parameter              :: nFaceA = 1
  integer (c_int), allocatable                :: faceVtxIdxA(:)
  integer (c_int), allocatable                :: faceVtxA(:)
  integer (kind = pdm_g_num_s), allocatable   :: faceLNToGNA(:)
  integer (c_int), parameter              :: nVtxA = 3
  double precision, allocatable               :: vtxCoordA(:)
  integer (kind = pdm_g_num_s), allocatable   :: vtxLNToGNA(:)

  integer (kind = pdm_g_num_s) :: nGOlFaceA
  integer (kind = pdm_g_num_s) :: nGOlVtxA

  integer (c_int) :: nOlFaceA
  integer (c_int) :: nOlLinkedFaceA
  integer (c_int) :: nOlVtxA
  integer (c_int) :: sOlFaceIniVtxA
  integer (c_int) :: sOlface_vtxA
  integer (c_int) :: sInitToOlFaceA

  integer (c_int), allocatable                :: olFaceIniVtxIdxA(:)
  integer (c_int), allocatable                :: olFaceIniVtxA(:)
  integer (c_int), allocatable                :: olface_vtx_idxA(:)
  integer (c_int), allocatable                :: olface_vtxA(:)
  integer (c_int), allocatable                :: olLinkedface_procIdxA(:)
  integer (c_int), allocatable                :: olLinkedFaceA(:)
  integer (kind = pdm_g_num_s), allocatable   :: olface_ln_to_gnA(:)
  double precision, allocatable               :: olCoordsA(:)
  integer (kind = pdm_g_num_s), allocatable   :: olvtx_ln_to_gnA(:)
  integer (c_int), allocatable                :: initToOlFaceIdxA(:)
  integer (c_int), allocatable                :: initToOlFaceA(:)


  integer (c_int), parameter              :: ipartB = 0
  integer (c_int), parameter              :: nFaceB = 1
  integer (c_int), allocatable                :: faceVtxIdxB(:)
  integer (c_int), allocatable                :: faceVtxB(:)
  integer (kind = pdm_g_num_s), allocatable   :: faceLNToGNB(:)
  integer (c_int), parameter              :: nVtxB = 3
  double precision, allocatable               :: vtxCoordB(:)
  integer (kind = pdm_g_num_s), allocatable   :: vtxLNToGNB(:)

  integer (kind = pdm_g_num_s) :: nGOlFaceB
  integer (kind = pdm_g_num_s) :: nGOlVtxB

  integer (c_int) :: nOlFaceB
  integer (c_int) :: nOlLinkedFaceB
  integer (c_int) :: nOlVtxB
  integer (c_int) :: sOlFaceIniVtxB
  integer (c_int) :: sOlface_vtxB
  integer (c_int) :: sInitToOlFaceB

  integer (c_int), allocatable                :: olFaceIniVtxIdxB(:)
  integer (c_int), allocatable                :: olFaceIniVtxB(:)
  integer (c_int), allocatable                :: olface_vtx_idxB(:)
  integer (c_int), allocatable                :: olface_vtxB(:)
  integer (c_int), allocatable                :: olLinkedface_procIdxB(:)
  integer (c_int), allocatable                :: olLinkedFaceB(:)
  integer (kind = pdm_g_num_s), allocatable   :: olface_ln_to_gnB(:)
  double precision, allocatable               :: olCoordsB(:)
  integer (kind = pdm_g_num_s), allocatable   :: olvtx_ln_to_gnB(:)
  integer (c_int), allocatable                :: initToOlFaceIdxB(:)
  integer (c_int), allocatable                :: initToOlFaceB(:)

  !
  ! Init
  !

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 1 MPI process is mandatory'
    call mpi_finalize(code)
    stop
  end if

  !
  ! Set MeshA
  !

  allocate (faceVtxIdxA(nFaceA+1))
  allocate (faceVtxA(nVtxA))
  allocate (faceLNToGNA(nFaceA))
  allocate (vtxCoordA(3*nVtxA))
  allocate (vtxLNToGNA(nVtxA))

  faceVtxIdxA(1) = 0
  faceVtxIdxA(2) = 3

  faceVtxA(1) = 1
  faceVtxA(2) = 2
  faceVtxA(3) = 3

  faceLNToGNA(1) = 1

  vtxCoordA(1) = 0.d0
  vtxCoordA(2) = 0.d0
  vtxCoordA(3) = 0.d0

  vtxCoordA(4) = 0.5d0
  vtxCoordA(5) = sqrt(3.d0)/2.d0
  vtxCoordA(6) = 0.d0

  vtxCoordA(7) = 1.d0
  vtxCoordA(8) = 0.d0
  vtxCoordA(9) = 0.d0

  vtxLNToGNA(1) = 1
  vtxLNToGNA(2) = 2
  vtxLNToGNA(3) = 3


  !
  ! Set MeshB
  !

  allocate (faceVtxIdxB(nFaceB+1))
  allocate (faceVtxB(nVtxB))
  allocate (faceLNToGNB(nFaceB))
  allocate (vtxCoordB(3*nVtxB))
  allocate (vtxLNToGNB(nVtxB))

  faceVtxIdxB(1) = 0
  faceVtxIdxB(2) = 3

  faceVtxB(1) = 1
  faceVtxB(2) = 2
  faceVtxB(3) = 3

  faceLNToGNB(1) = 1

  vtxCoordB(1) = 0.5d0
  vtxCoordB(2) = -sqrt(3.d0)/4.d0
  vtxCoordB(3) = 0.d0

  vtxCoordB(4) = 1.d0
  vtxCoordB(5) = sqrt(3.d0)/4.d0
  vtxCoordB(6) = 0.d0

  vtxCoordB(7) = 0.d0
  vtxCoordB(8) = sqrt(3.d0)/4.d0
  vtxCoordB(9) = 0.d0

  vtxLNToGNB(1) = 1
  vtxLNToGNB(2) = 2
  vtxLNToGNB(3) = 3

  !
  ! Create a new PDM_mesh_location structure
  !   The MPI communicator and the number of point cloud are setted
  !

  call PDM_ol_create (n_partMeshA, &
                      nGFaceMeshA, &
                      nGVtxMeshA, &
                      n_partMeshB, &
                      nGFaceMeshB, &
                      nGVtxMeshB, &
                      projectCoeff, &
                      MPI_COMM_WORLD, &
                      id)

  call PDM_ol_parameter_set (id, &
                             PDM_OL_CAR_LENGTH_TOL, &
                             1.d-4)

  call PDM_ol_parameter_set (id, &
                             PDM_OL_EXTENTS_TOL, &
                             1.d-4)

  call PDM_ol_input_mesh_set (id, &
                              PDM_OL_MESH_A, &
                              ipartA, &
                              nFaceA, &
                              faceVtxIdxA, &
                              faceVtxA, &
                              faceLNToGNA, &
                              nVtxA, &
                              vtxCoordA, &
                              vtxLNToGNA)

  call PDM_ol_input_mesh_set (id, &
                              PDM_OL_MESH_B, &
                              ipartB, &
                              nFaceB, &
                              faceVtxIdxB, &
                              faceVtxB, &
                              faceLNToGNB, &
                              nVtxB, &
                              vtxCoordB, &
                              vtxLNToGNB)

  call PDM_ol_compute (id);

  call PDM_ol_dump_times (id);

  call PDM_ol_mesh_dim_get (id, &
                            PDM_OL_MESH_A, &
                            nGOlFaceA, &
                            nGOlVtxA)

  call PDM_ol_part_mesh_dim_get (id, &
                                 PDM_OL_MESH_A, &
                                 ipartA, &
                                 nOlFaceA, &
                                 nOlLinkedFaceA, &
                                 nOlVtxA,&
                                 sOlFaceIniVtxA,&
                                 sOlface_vtxA, &
                                 sInitToOlFaceA);

  allocate(olFaceIniVtxIdxA(nFaceA+1))
  allocate(olFaceIniVtxA(sOlFaceIniVtxA))
  allocate(olface_vtx_idxA(nOlFaceA+1))
  allocate(olface_vtxA(sOlface_vtxA))
  allocate(olLinkedface_procIdxA(n_rank+1))
  allocate(olLinkedFaceA(4*nOlLinkedFaceA))
  allocate(olface_ln_to_gnA(nOlFaceA))
  allocate(olCoordsA(3*nOlVtxA))
  allocate(olvtx_ln_to_gnA(nOlVtxA))
  allocate(initToOlFaceIdxA(nFaceA+1))
  allocate(initToOlFaceA(sInitToOlFaceA))

  call PDM_ol_mesh_entities_get (id, &
                                 PDM_OL_MESH_A, &
                                 ipartA, &
                                 olFaceIniVtxIdxA, &
                                 olFaceIniVtxA, &
                                 olface_vtx_idxA, &
                                 olface_vtxA, &
                                 olLinkedface_procIdxA, &
                                 olLinkedFaceA, &
                                 olface_ln_to_gnA, &
                                 olCoordsA, &
                                 olvtx_ln_to_gnA, &
                                 initToOlFaceIdxA, &
                                 initToOlFaceA)


  call PDM_ol_mesh_dim_get (id, &
                            PDM_OL_MESH_B, &
                            nGOlFaceB, &
                            nGOlVtxB)

  call PDM_ol_part_mesh_dim_get (id, &
                                 PDM_OL_MESH_B, &
                                 ipartB, &
                                 nOlFaceB, &
                                 nOlLinkedFaceB, &
                                 nOlVtxB,&
                                 sOlFaceIniVtxB,&
                                 sOlface_vtxB, &
                                 sInitToOlFaceB);


  allocate(olFaceIniVtxIdxB(nFaceB+1))
  allocate(olFaceIniVtxB(sOlFaceIniVtxB))
  allocate(olface_vtx_idxB(nOlFaceB+1))
  allocate(olface_vtxB(sOlface_vtxB))
  allocate(olLinkedface_procIdxB(n_rank+1))
  allocate(olLinkedFaceB(4*nOlLinkedFaceB))
  allocate(olface_ln_to_gnB(nOlFaceB))
  allocate(olCoordsB(3*nOlVtxB))
  allocate(olvtx_ln_to_gnB(nOlVtxB))
  allocate(initToOlFaceIdxB(nFaceB+1))
  allocate(initToOlFaceB(sInitToOlFaceB))

  call PDM_ol_mesh_entities_get (id, &
                                 PDM_OL_MESH_B, &
                                 ipartB, &
                                 olFaceIniVtxIdxB, &
                                 olFaceIniVtxB, &
                                 olface_vtx_idxB, &
                                 olface_vtxB, &
                                 olLinkedface_procIdxB, &
                                 olLinkedFaceB, &
                                 olface_ln_to_gnB, &
                                 olCoordsB, &
                                 olvtx_ln_to_gnB, &
                                 initToOlFaceIdxB, &
                                 initToOlFaceB)


  call PDM_ol_del (id)

  deallocate(olFaceIniVtxIdxA)
  deallocate(olFaceIniVtxA)
  deallocate(olface_vtx_idxA)
  deallocate(olface_vtxA)
  deallocate(olLinkedface_procIdxA)
  deallocate(olLinkedFaceA)
  deallocate(olface_ln_to_gnA)
  deallocate(olCoordsA)
  deallocate(olvtx_ln_to_gnA)
  deallocate(initToOlFaceIdxA)
  deallocate(initToOlFaceA)

  deallocate(olFaceIniVtxIdxB)
  deallocate(olFaceIniVtxB)
  deallocate(olface_vtx_idxB)
  deallocate(olface_vtxB)
  deallocate(olLinkedface_procIdxB)
  deallocate(olLinkedFaceB)
  deallocate(olface_ln_to_gnB)
  deallocate(olCoordsB)
  deallocate(olvtx_ln_to_gnB)
  deallocate(initToOlFaceIdxB)
  deallocate(initToOlFaceB)

  deallocate (faceVtxIdxA)
  deallocate (faceVtxA)
  deallocate (faceLNToGNA)
  deallocate (vtxCoordA)
  deallocate (vtxLNToGNA)

  deallocate (faceVtxIdxB)
  deallocate (faceVtxB)
  deallocate (faceLNToGNB)
  deallocate (vtxCoordB)
  deallocate (vtxLNToGNB)

  call mpi_finalize(code)


end program testf
