!-----------------------------------------------------------------------------
! This file is part of the CWIPI library. 
!
! Copyright (C) 2011  ONERA
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

subroutine printStatus(iiunit, status)
  use cwipi
  
  implicit none
  integer :: status
  integer :: iiunit

  select case (status)
     case(cwipi_exchange_ok)
        write(iiunit,*) "Exchange ok"
     case(cwipi_exchange_bad_receiving)
        write(iiunit,*) "no or bad receiving"
     case default
        write(iiunit,*) "Unknown receiving status"
        stop
  end select

end subroutine printStatus

!
! ------------------------------------------------------------------------------
! Fonction d'interpolation bidon, juste pour voir si c'est bien pris en compte
! ------------------------------------------------------------------------------
!

subroutine  userInterpolation(entitiesDim, &
                              nLocalVertex, &
                              nLocalElement, &
                              nLocalPolyhedra, &
                              nDistantPoint, &
                              localCoordinates, &
                              localConnectivityIndex, &
                              localConnectivity, &
                              localPolyFaceIndex, &
                              localPolyCellToFaceConnec, &
                              localPolyFaceConnecIdx, &
                              localPolyFaceConnec, &
                              disPtsCoordinates, &
                              disPtsLocation, &
                              disPtsDistance, &
                              disPtsBaryCoordIdx, &
                              disPtsBaryCoord, &
                              stride, &
                              solverType, &
                              localField, &
                              distantField)

  use cwipi

  implicit none

  integer :: entitiesDim
  integer :: nLocalVertex
  integer :: nLocalElement
  integer :: nLocalPolyhedra
  integer :: nDistantPoint
  real(8) :: localCoordinates(*)
  integer :: localConnectivityIndex(*)
  integer :: localConnectivity(*)
  integer :: localPolyFaceIndex(*)
  integer :: localPolyCellToFaceConnec(*)
  integer :: localPolyFaceConnecIdx(*)
  integer :: localPolyFaceConnec(*)
  real(8) :: disPtsCoordinates(*)
  integer :: disPtsLocation(*)
  real(4) :: disPtsDistance(*)
  integer :: disPtsBaryCoordIdx(*)
  real(8) :: disPtsBaryCoord(*)
  integer :: stride
  integer :: solverType
  real(8) :: localField(*)
  real(8) :: distantField(*)

  integer :: i

  if (solverType .eq. cwipi_solver_cell_center) then
     do i = 1, nDistantPoint
        distantField(i) = localField(disPtsLocation(i))
     enddo
  else
     print*, 'Error in _userInterpolation : bad solver_type'
     stop
  endif

end subroutine userInterpolation

!
! -----------------------------------------------------------------------------
! Programme de tests
! -----------------------------------------------------------------------------
!

program testf

  use mpi
  
  use cwipi
  
  use modDeterminant
  use baseSimplex2D
  use baseSimplex3D

  implicit none

  interface
       subroutine  userInterpolation(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer :: entitiesDim
         integer :: nLocalVertex
         integer :: nLocalElement
         integer :: nLocalPolyhedra
         integer :: nDistantPoint
         real(8) :: localCoordinates(*)
         integer :: localConnectivityIndex(*)
         integer :: localConnectivity(*)
         integer :: localPolyFaceIndex(*)
         integer :: localPolyCellToFaceConnec(*)
         integer :: localPolyFaceConnecIdx(*)
         integer :: localPolyFaceConnec(*)
         real(8) :: disPtsCoordinates(*)
         integer :: disPtsLocation(*)
         real(4) :: disPtsDistance(*)
         integer :: disPtsBaryCoordIdx(*)
         real(8) :: disPtsBaryCoord(*)
         integer :: stride
         integer :: solverType
         real(8) :: localField(*)
         real(8) :: distantField(*)
       end subroutine userInterpolation
    end interface

!  integer, pointer :: location(:)
!  integer, pointer :: baryCooIdx(:)
!  real(8), pointer :: baryCoo(:)
!  real(8), pointer :: tmpDbl(:)
  integer :: nLocatedPoints
  integer :: nNotLocatedPoints
  integer :: nDistantPoints

  integer :: localcom, localGroup, p1Group, p1Comm
  integer :: iRank, currentRank, localcommsize
  character (len = 4) :: proc
  character (len = 5) :: codeName, codeCoupledName
  integer :: code
  integer :: iiunit
  integer :: ivalue
  real(8) :: dvalue

  real(8) :: xmin = -100.d0
  real(8) :: xmax =  100.d0
  real(8) :: ymin = -100.d0
  real(8) :: ymax =  100.d0
  integer :: nx   = 24
  integer :: ny   = 28
  integer :: initrandom = 2

  integer nVert, nCell, lconnecindex
  character(10)      :: name
  
  integer, parameter :: nVertm = 4000
  integer, parameter :: nCellm = 4000
  integer, parameter :: lconnecindexm = 12000
  integer, parameter :: nptstolocate = 21

  real(8), pointer :: coordstolocate(:)

  real(8), pointer :: vertx(:)
  integer, pointer :: connecindex(:)
  integer, pointer :: connec     (:)

  real(8), pointer :: values(:)
  real(8), pointer :: localvalues(:)

  real(4), pointer :: distLocPts(:)

  integer status

  integer i, order, k

  integer :: vpar = 10
  character (len =  6) :: cpar = "niterf"
  character (len = 30) :: disstr = ""

  integer :: stride = 1
  integer :: rl(1)
  integer :: dislocalcommsize
  integer :: commWorldSize

  integer :: n_partition, n2, codeId
  
  integer :: nVertSeg

  integer :: nLocatedPts

  real(8) :: randLevel

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, iRank, code)
  call mpi_comm_size(mpi_comm_world, commWorldSize, code)

  write(proc,'(i4.4)') iRank
  iiunit = 9
  
  open(unit=iiunit, file='fortran_surf_PiPj_'//proc//'.txt', &
       form='formatted', status='unknown')

  if (iRank == 0) then
     print '(/"START: fortran_surf_PiPj")'
  endif

  n_partition = 1
  do while ((2 * n_partition**2) < commWorldSize)
     n_partition = n_partition + 1
  enddo

  n2 = 2 * n_partition**2
  
  if (n2 /= commWorldSize) then
     if (iRank == 0) then
        print *, '      Not executed : only available if the number of processus in the form of 2 * n_partition**2'
     endif
     call mpi_finalize(code)
     stop;
  endif


  ! -----------------------------------------
  ! Initialisation des maillages
  ! -----------------------------------------

  !>  Vertices
  !>  5
  !>  0.   0.   0.    1
  !>  1.   0.   0.    1
  !>  0.   1.   0.    1
  !>  0.   0.   1.    1
  !>  0.73 0.73 0.73  1
  
  !>  Tetrahedra
  !>  2
  !>  1 2 3 4  1
  !>  5 2 4 3  1


  if( iRank==0 )then
    
    nVert=4
    nCell=1
    allocate( vertx(1:3*nVert) )  !> 4 sommets
    allocate( connecindex(1:2) )  !> 1 tetra
    allocate( connec     (1:4) )  !> 1 tetra
    
    coord(01:03)=[0d0,0d0,0d0]    
    coord(04:06)=[1d0,0d0,0d0]    
    coord(07:09)=[0d0,1d0,0d0]    
    coord(10:12)=[0d0,0d0,1d0]    
    
    connec(1:4)=[1,2,3,4]
    connecindex(1:2)=[0,4]
    
  else( iRank==1 )then
  
    nVert=1
    nCell=1
    allocate( vertx(1:4*3) )  !> 4 sommets
    allocate( connecindex(1:2  ) )  !> 1 tetra
    allocate( connec     (1:4  ) )  !> 1 tetra
    
    coord(01:03)=[0.73d0,0.73d0,0.73d0]    
    coord(04:06)=[1.00d0,0.00d0,0.00d0]    
    coord(07:09)=[0.00d0,0.00d0,1.00d0]    
    coord(10:12)=[0.00d0,1.00d0,0.00d0]    
    
    connec(1:4)=[1,2,3,4]
    connecindex(1:2)=[0,4]
    
  endif
  
  !> Ecriture des maillages au format mesh de l'inria
  write(name,'("Tetra",i1,".mesh")')iRank
  open(unit=100,file=trim(name),action='write',status='unknown')
  write(100,'("MeshVersionFormatted 1"/)')
  write(100,'("Dimension 3"/)')
  write(100,'("Vertices")')
  write(100,*)nVert
  j=0
  do iVert=1,nVert
    write(100,'(3(e22.15,1x),i2)')coord(j+1:j+3),0  ; j=j+3
  enddo
  write(100,'(/"Tetrahedra")')
  write(100,*)nCell
  do iCell=1,nCell
    write(100,'(*(i6,1x))')connec(connecindex(iCell)+1:connecindex(iCell+1)),0
  enddo
  write(100,'(/"End")')
  close(100)



!
! -----------------------------------------
! Initialisation de l'interface de couplage
! -----------------------------------------
!
 
  call cwipi_set_output_listing_f(iiunit)

  if (iRank < commWorldSize / 2) then
     codeName = 'code1'
     codeId = 1
     codeCoupledName = 'code2'
  else 
     codeName = 'code2'
     codeId = 2
     codeCoupledName = 'code1'
  endif

  call cwipi_init_f (mpi_comm_world, &
                     codeName, &
                     localcom)

!
! ------------------------------------------------
! Creation du fichier de sortie listing
! (en parallele, un fichier listing par processus)
! ------------------------------------------------
!

  call mpi_comm_rank(localcom, currentRank, code)
  call mpi_comm_size(localcom, localcommsize, code)

  write(iiunit,*)
  write(iiunit,*) "dump apres initialisation"
  write(iiunit,*) "-------------------------"
  write(iiunit,*)

  call cwipi_dump_appli_properties_f

!
! -------------------------------------
! Test de definition des points
! a interpoler
! -------------------------------------
!
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 3"
  write(iiunit, *)

  if (iRank == 0) then
     print*, '       Create coupling'  
  endif

  call cwipi_create_coupling_f("test2D_3", &
                               cwipi_cpl_parallel_with_part,&
                               codeCoupledName, &
                               2,     & ! Dimension des entites geometriques
                               0.1d0, & ! Tolerance geometrique
                               cwipi_static_mesh, &
                               cwipi_solver_cell_center, &
                               1, &
                               "Ensight Gold",&
                               "text")

!
! Construction du maillage

  nVertSeg   = 10
  randLevel    = 0.1d0
  nVert      = nVertSeg * nVertSeg
  nCell        = (nVertSeg - 1) * (nVertSeg - 1)

  allocate(vertx(3 * nVert))
  allocate(coordstolocate(3 * nVert))

  allocate(connecindex(nCell + 1))
  allocate(connec(4 * nCell))

  allocate(values(nCell))
  allocate(localvalues(nVert))

  if (iRank == 0) then
     print*, '       Create mesh'
  endif

  call grid_mesh_f(xmin, &
                   xmax, &
                   ymin, &
                   ymax, &
                   randLevel, &
                   nVertSeg, &
                   n_partition, & 
                   coords,  &
                   connecindex,&
                   connec,&
                   localcom)

  call cwipi_define_mesh_f("test2D_3", &
                           nVert, &
                           nCell, &
                           coords, &
                           connecindex, &
                           connec)

!
! Definition des points a localiser

  do i = 1, 3 * nVert
     coordstolocate(i) = 0.75 * vertx(i)
  enddo

  if (iRank == 0) then
     print*, '       Set points to locate'
  endif

  call cwipi_set_points_to_locate_f("test2D_3", &
                                     nptstolocate, &
                                     coordstolocate)


!
! Envoi de la coordonnee Y a codeC
! Reception de la coordonnee Y provenant de codec

  do i = 1, nCell
     if (codeId == 1) then
        values(i) = vertx(3*(i-1) + 1)
     else
        values(i) = vertx(3*(i-1) + 2)
     endif
  enddo

  stride = 1

  if (iRank == 0) then
     print*, '       Exchange'
  endif

  call cwipi_exchange_f ("test2D_3", &
                         "echange1", &
                         stride, &
                         1, &
                         0.1d0, &
                         "cooy", &
                         values, &
                         "coox", &
                         localvalues, &
                         userInterpolation, &
                         nNotLocatedPoints, &
                         status)
  
  call cwipi_get_n_located_pts_f("test2D_3", nLocatedPts)

  allocate(distLocPts(nLocatedPts))

  call cwipi_dist_located_pts_get_f("test2D_3", distLocPts)

  call printStatus(iiunit, status)
  write(iiunit,*) "valeurs recues test2D_3"
  write(iiunit,*) (localvalues(i),i=1,nptstolocate)

!
! Suppression de l'objet couplage

  if (iRank == 0) then
     print*, '       Delete coupling'
  endif

  call cwipi_delete_coupling_f("test2D_3");
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

  call cwipi_finalize_f()

  deallocate(coords)
  deallocate(coordstolocate)
  deallocate(connecindex)
  deallocate(connec)
  deallocate(values)
  deallocate(localvalues)

  call mpi_finalize(code)

end program testf
