/*============================================================================
 * Functions about high order meshes
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2018       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_error.h>
#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_ho.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_triangulate.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_ho.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _fvmc_ho_user_fcts_t {
  fvmc_ho_location_fct_t location_tetra;
  fvmc_ho_location_fct_t location_prism;
  fvmc_ho_location_fct_t location_pyramid;
  fvmc_ho_location_fct_t location_hexa;
  fvmc_ho_location_fct_t location_tria;
  fvmc_ho_location_fct_t location_quad;
  fvmc_ho_location_fct_t location_edge;
  
  fvmc_ho_weight_fct_t weight_tetra;
  fvmc_ho_weight_fct_t weight_prism;
  fvmc_ho_weight_fct_t weight_pyramid;
  fvmc_ho_weight_fct_t weight_hexa;
  fvmc_ho_weight_fct_t weight_tria;
  fvmc_ho_weight_fct_t weight_quad;
  fvmc_ho_weight_fct_t weight_edge;


  fvmc_ho_interp_fct_t interp_tetra;
  fvmc_ho_interp_fct_t interp_prism;
  fvmc_ho_interp_fct_t interp_pyramid;
  fvmc_ho_interp_fct_t interp_hexa;
  fvmc_ho_interp_fct_t interp_tria;
  fvmc_ho_interp_fct_t interp_quad;
  fvmc_ho_interp_fct_t interp_edge;

} fvmc_ho_user_fcts_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static fvmc_ho_user_fcts_t *_user_fcts = NULL;


/*============================================================================
 * Private function definitions
 *============================================================================*/


/*   subroutine outsideTrian2(trian2,xyz, uv,max_dist2) */
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* #define outsideTrian2 0 */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     ! Calcul distance point xyz à trian2 (1,2,3,4,5,6) */
/*     ! 3 */
/*     ! 6 5 */
/*     ! 1 4 2 */
/*     ! 4 triangles : */
/*     ! 6    5    3    6 5 */
/*     ! 1 4  4 2  6 5  4  */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     use space_cell_vertices, only: uvTriangle */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     real(8), intent(in)    :: trian2(1:3,1:6) */
/*     real(8), intent(in)    :: xyz  (1:3) */
/*     real(8), intent(out)   :: uv   (1:2) */
/*     real(8), intent(out)   :: max_dist2 */
/*     !> */
/*     integer                :: iVert,nVert */
/*     integer                :: nTria,iTria,tria1 (1:3,1:4) */
/*     integer                :: nLine,iLine,lines (1:3,1:9) */
/*     real(8)                :: dist2    , dist2_ (    1:9) */
/*     logical                :: inside   , inside_(    1:9) */
/*     real(8)                :: uv1 (1:2), uv1_   (1:2,1:9) */
/*     real(8)                :: trian1(1:3,1:3) */
/*     real(8)                :: xyzP(1:3) */
/*     real(8)                :: ai  (1:6) */
/*     real(8)                :: pt0(1:3) */
/*     real(8)                :: pt1(1:3),pt2(1:3),pt3(1:3) */
/*     real(8)                :: u */
/*     real(8)                :: dxyz     (    1:3) */
/*     real(8)                :: dist2Vert(    1:6) */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* #if outsideTrian2!=0 */
/*     print '(">>> outsideTrian2")' */
/*     do iVert=1,6 */
/*       print '(4x,"outsideTrian2 trian",i1,"=",3(f12.5,1x))',iVert,trian2(1:3,iVert) */
/*     enddo */
/*     print '(4x,"outsideTrian2    xyz=",3(f12.5,1x))',xyz(1:3) */
/* #endif */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     max_dist2=hugeD2 */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     !> Triangle P2 */
/*     !>  3 */
/*     !>  6 5 */
/*     !>  1 4 2 */
/*     !> */
/*     !> -> 4 Triangles P1 */
/*     !> 6    5    3    6 5 */
/*     !> 1 4  4 2  6 5  4 */
    
/*     tria1(1:3,1)=[1,4,6] */
/*     tria1(1:3,2)=[4,2,5] */
/*     tria1(1:3,3)=[6,5,3] */
/*     tria1(1:3,4)=[4,5,6] */
/*     nTria=4 */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     !> xyz is projected on Triangles P1 => xyzP */
/*     !> if xyzP \in TriangleP1 dist2=|xyz,xyzP|^2 */
/*     inside_(1:nTria)=.false. */
/*     loopTriangles: do iTria=1,nTria */
/*       pt1(1:3)=trian2(1:3,tria1(1,iTria)) */
/*       pt2(1:3)=trian2(1:3,tria1(2,iTria)) */
/*       pt3(1:3)=trian2(1:3,tria1(3,iTria)) */

/*                         // Projection du point sur le plan du triangle                   */
/*       call projectPlan3D(pt1=pt1,pt2=pt2,pt3=pt3,xyz=xyz(1:3), xyzP=xyzP(1:3),dist2=dist2) */
      
/*       trian1(1:3,1)=pt1(1:3) */
/*       trian1(1:3,2)=pt2(1:3) */
/*       trian1(1:3,3)=pt3(1:3) */
/*                         // Dans triangle P1 courant */
/*       call trian1_xyz_uv(trian1=trian1,xyz=xyzP, uv=uv1) */
      
/*       call baseTrianP1(u=uv1(1),v=uv1(2),ai=ai(1:3)) ! print '("outsideTrian2 xyzP=",3(e12.5,1x)," ai=",3(e12.5,1x),"dist2=",e12.5)',xyzP(1:3),ai(1:3),dist2 */
/*       if( 0d0<=ai(1) .and. 0d0<=ai(2) .and. 0d0<=ai(3) )then */
/*         dist2_  (iTria)=dist2 */
/*         inside_ (iTria)=.true. */

/*           // uv1_ coeff interpolés dans trian                                                    */
/*         uv1_(1:2,iTria)= uvTriangle(1:2,tria1(1,iTria))*ai(1) & */
/*         &               +uvTriangle(1:2,tria1(2,iTria))*ai(2) & */
/*         &               +uvTriangle(1:2,tria1(3,iTria))*ai(3) */
       
/* #if outsideTrian2=!0 */
/*         call baseTrianP2(u=uv1_(1,iTria),v=uv1_(2,iTria),ai=ai(1:6)) */
/*         pt0(1:3)=matmul(trian2(1:3,1:6),ai(1:6)) */
/*         !if( 1d-15<norm2(pt0(1:3)-xyzP(1:3)) )then */
/*           print '(/4x,"outsideTrian2 2D Lagrange")' */
/*           print '( 4x,"outsideTrian2 uv  =",3(e12.5,1x))',uv(1:2) */
/*           print '( 4x,"outsideTrian2 xyzP=",3(e12.5,1x))',xyzP(1:3) */
/*           print '( 4x,"outsideTrian2 pt0 =",3(e12.5,1x))',pt0(1:3) */
/*           print '( 4x,"outsideTrian2 dist(pt0,xyzP)=",e12.5)',norm2(pt0(1:3)-xyzP(1:3)) */
/*          !print '(4x,"Stop File: ",a," line: ",i6)',__FILE__,__LINE__+1 ; stop */
/*         !endif */
/* #endif */
        
/*       else */
/*         inside_(iTria)=.false. */
/*       endif */
      
/* #if outsideTrian2!=0 */
/*       if( inside_(iTria) )then */
/*         print '(4x,"outsideTrian2:projectPlan3D in: ",l," xyz=",3(e12.5,1x),"xyzP=",3(e12.5,1x),"dist2=",e12.5)',inside_(iTria),xyz(1:3),xyzP(1:3),dist2 */
/*       endif */
/* #endif */
/*     enddo loopTriangles */
    
/*     !> max_dist2 = min( |xyz,xyzP|^2 ) where inside_(:)==.true. */
/*     do iTria=1,nTria */
/*       if( inside_(iTria) )max_dist2=min(max_dist2,dist2_(iTria)) */
/*     enddo */
    
/*     !> Keeping iEdge where inside_=true and dist2=max_dist2 */
/*     do iTria=1,nTria */
/*       if( inside_(iTria) .and. max_dist2==dist2_(iTria) )then */
/*         uv(1:2)=uv1_(1:2,iTria) */
/*       endif */
/*     enddo */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     !> Connectivity of lines */
/*     lines(1:2,1)=[1,4] */
/*     lines(1:2,2)=[2,4] */
/*     lines(1:2,3)=[2,5] */
/*     lines(1:2,4)=[5,3] */
/*     lines(1:2,5)=[3,6] */
/*     lines(1:2,6)=[6,1] */
/*     lines(1:2,7)=[4,5] !> internal line */
/*     lines(1:2,8)=[5,6] !> internal line */
/*     lines(1:2,9)=[6,4] !> internal line */
/*     nLine=9 */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     !> point is projected on trian2 lines */
    
/*     inside_(1:nLine)=.false. */
/*     loopLines: do iLine=1,nLine */
/*       pt1(1:3)=trian2(1:3,lines(1,iLine)) */
/*       pt2(1:3)=trian2(1:3,lines(2,iLine)) */
      
/*       call projectLine3D(pt1=pt1(1:3),pt2=pt2(1:3),xyz=xyz(1:3) ,xyzP=xyzP(1:3),dist2=dist2) */
      
/*       pt2(1:3)=pt2 (1:3)-pt1(1:3) */
/*       pt0(1:3)=xyzP(1:3)-pt1(1:3) */
/*       u=dot_product(pt0,pt2)/dot_product(pt2,pt2) ; u=2d0*u-1d0 ! u \in [-1,+1] */
/*       call baseEdgeP1(u=u,ai=ai(1:2)) */
/*       if( 0d0<=ai(1) .and. 0d0<=ai(2) )then */
/*         dist2_ (    iLine)=dist2 */
/*         uv1_   (1:2,iLine)= uvTriangle(1:2,lines(1,iLine))*ai(1) & */
/*         &                  +uvTriangle(1:2,lines(2,iLine))*ai(2) */
/*         inside_(    iLine)=.true. */
        
/* #if outsideTrian2=!0 */
/*         call baseTrianP2(u=uv1_(1,iLine),v=uv1_(2,iLine),ai=ai) */
/*         pt0(1:3)=matmul(trian2(1:3,1:6),ai(1:6)) */
        
/*         if( 1d-15<norm2(pt0(1:3)-xyzP(1:3)) )then */
/*           pt1(1:3)=trian2(1:3,lines(1,iLine)) */
/*           pt2(1:3)=trian2(1:3,lines(2,iLine)) */
          
/*           pt0(1:3)=5d-1*((1d0-u)*pt1(1:3)+(1+u)*pt2(1:3)) */
/*           print '(/4x,"outsideTrian2 1D Lagrange")' */
/*           print '( 4x,"outsideTrian2 u   =",1(e12.5,1x))',u */
/*           print '( 4x,"outsideTrian2 pt1 =",3(e12.5,1x))',pt1(1:3) */
/*           print '( 4x,"outsideTrian2 pt2 =",3(e12.5,1x))',pt2(1:3) */
/*           print '( 4x,"outsideTrian2 xyzP=",3(e12.5,1x))',xyzP(1:3) */
/*           print '( 4x,"outsideTrian2 pt0 =",3(e12.5,1x))',pt0(1:3) */
/*           print '( 4x,"outsideTrian2 dist(pt0,xyzP)=",e12.5)',norm2(pt0(1:3)-xyzP(1:3)) */
          
/*           pt0(1:3)=matmul(trian2(1:3,1:6),ai(1:6)) */
/*           print '(/4x,"outsideTrian2 2D Lagrange")' */
/*           print '( 4x,"outsideTrian2 uv  =",3(e12.5,1x))',uv1_(1:2,iLine) */
/*           print '( 4x,"outsideTrian2 pt1 =",3(e12.5,1x))',trian2(1:3,1) */
/*           print '( 4x,"outsideTrian2 pt2 =",3(e12.5,1x))',trian2(1:3,2) */
/*           print '( 4x,"outsideTrian2 pt3 =",3(e12.5,1x))',trian2(1:3,3) */
/*           print '( 4x,"outsideTrian2 pt4 =",3(e12.5,1x))',trian2(1:3,4) */
/*           print '( 4x,"outsideTrian2 pt5 =",3(e12.5,1x))',trian2(1:3,5) */
/*           print '( 4x,"outsideTrian2 pt6 =",3(e12.5,1x))',trian2(1:3,6) */
/*           print '( 4x,"outsideTrian2 xyzP=",3(e12.5,1x))',xyzP(1:3) */
/*           print '( 4x,"outsideTrian2 pt0 =",3(e12.5,1x))',pt0(1:3) */
/*           print '( 4x,"outsideTrian2 dist(pt0,xyzP)=",e12.5)',norm2(pt0(1:3)-xyzP(1:3)) */
          
/*          !print '(4x,"Stop File: ",a," line: ",i6)',__FILE__,__LINE__+1 ; stop */
/*         endif */
/* #endif */
        
/*       endif */
      
/* #if outsideTrian2!=0 */
/*       if( inside_(iLine) )then */
/*         print '(4x,"outsideTrian2:projectLine3D in: ",l," xyz=",3(e12.5,1x),"xyzP=",3(e12.5,1x),"dist2=",e12.5)',inside_(iLine),xyz(1:3),xyzP(1:3),dist2 */
/*       endif */
/* #endif */
/*     enddo loopLines */
    
/*     !> Computing max_dist2 */
/*     do iLine=1,nLine */
/*       if( inside_(iLine) )max_dist2=min(max_dist2,dist2_(iLine)) */
/*     enddo */
    
/*     !> Keeping iEdge where inside_=true and dist2=max_dist2 */
/*     do iLine=1,nLine */
/*       if( inside_(iLine) .and. max_dist2==dist2_(iLine) )then */
/*         uv (1:2)=uv1_(1:2,iLine) */
/*       endif */
/*     enddo */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     !> min distance(point, tetra vertices) */
    
/*     nVert=size(trian2,2) */
/*     do iVert=1,nVert */
/*       dxyz(1:3)=xyz(1:3)-trian2(1:3,iVert) */
/*       dist2Vert(iVert)=dot_product(dxyz,dxyz) */
/*      !print '("iVert=",i1," dist2Vert=",e12.5)',iVert,dist2Vert(iVert) */
/*     enddo */
    
/*     dist2=minval(dist2Vert(1:nVert)) */
/*     do iVert=1,nVert */
/*       if( dist2Vert(iVert)==dist2 )then */
/*         uv1 (1:2)=uvTriangle(1:2,iVert) */
/*         xyzP(1:3)=trian2    (1:3,iVert) */
/*       endif */
/*     enddo */
    
/* #if outsideTrian2!=0 */
/*     print '(/4x,"outsideTrian2:projectVert3D in: ",l," xyz=",3(e12.5,1x),"xyzP=",3(e12.5,1x),"dist2=",e12.5)',.true.,xyz(1:3),xyzP(1:3),dist2 */
/* #endif */
     
/*     if( dist2<max_dist2 )then */
/*       max_dist2=dist2 */
/*       uv(1:2)=uv1(1:2) */
      
/* #if outsideTrian2!=0 */
/*       print '(4x,"outsideTrian2:projectVert3D xyz=",3(e12.5,1x),"xyzP=",3(e12.5,1x),"dist2=",e12.5)',xyz(1:3),xyzP(1:3),max_dist2 */
/* #endif */
      
/*     endif */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* #if outsideTrian2!=0 */
/*    !print '(4x,"outsideTrian2 max_dist2=",e12.5," uvw=",3(e12.5,1x))',max_dist2,uvw(1:3) */
/*     print '(4x,"outsideTrian2 max_dist2=",e12.5)',max_dist2 */
/*     print '("<<< outsideTrian2")' */
/* #endif */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* #undef outsideTrian2 */
/*     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*     return */
/*   end subroutine outsideTrian2 */




/*----------------------------------------------------------------------------
 * 
 * Default point location in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates if outside (or NULL)
 * 
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_in_cell_3d (const fvmc_element_t type,
                             const int order,
                             const int n_node,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_in_cell_3d : Not implemented yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_on_cell_2d (const fvmc_element_t type,
                             const int order,
                             const int n_node,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_on_cell_2d : Not implemented yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_on_cell_1d (const fvmc_element_t type,
                             const int order,
                             const int n_node,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_on_cell_1d : Not implemented yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Default weight computataion in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point inside cell  
 *   weight            --> weights
 * 
 *----------------------------------------------------------------------------*/

static void
_default_weight_in_cell_3d (const fvmc_element_t type,
                           const int order,
                           const int n_node,
                           const int *ho_vertex_num,
                           const double *vertex_coords,
                           const double *point_coords,
                           double *weight)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_weight_in_cell_3d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default weight computataion on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point on cell  
 *   weight            --> weights
 * 
 *----------------------------------------------------------------------------*/

static void 
 _default_weight_on_cell_2d (const fvmc_element_t type,
                            const int order,
                            const int n_node,
                            const int *ho_vertex_num,
                            const double *vertex_coords,
                            const double *point_coords,
                            double *weight)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_weight_on_cell_2d : Not implement yet\n"));
}

/*----------------------------------------------------------------------------
 * 
 * Default weight computataion on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point on cell  
 *   weight            --> weights
 * 
 *----------------------------------------------------------------------------*/

static void 
_default_weight_on_cell_1d (const fvmc_element_t type,
                           const int order,
                           const int n_node,
                           const int *ho_vertex_num,
                           const double *vertex_coords,
                           const double *point_coords,
                           double *weight)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_weight_on_cell_1d : Not implemented yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords in a 3D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

static void 
_default_interp_in_cell_3d (const fvmc_element_t type,
                            const int order,
                            const int n_node,
                            const int *ho_vertex_num,
                            const int *local_to_user,
                            const double *vertex_coords,
                            const double *point_coords,
                            const float *distance,
                           const double *point_proj_coords,
                            const double *weight,
                            const int stride_field,
                            const double *src_field,
                            double *target_field)
                           
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_interp_in_cell_3d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

static void 
_default_interp_on_cell_2d (const fvmc_element_t type,
                            const int order,
                            const int n_node,
                            const int *ho_vertex_num,
                            const int *local_to_user,
                            const double *vertex_coords,
                            const double *point_coords,
                            const float *distance,
                           const double *point_proj_coords,
                            const double *weight,
                            const int stride_field,
                            const double *src_field,
                            double *target_field)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_interp_on_cell_2d : Not implemented yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords on a 1D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

static void 
_default_interp_on_cell_1d (const fvmc_element_t type,
                            const int order,
                            const int n_node,
                            const int *ho_vertex_num,
                            const int *local_to_user,
                            const double *vertex_coords,
                            const double *point_coords,
                            const float *distance,
                            const double *point_proj_coords,
                            const double *weight,
                            const int stride_field,
                            const double *src_field,
                            double *target_field)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_weight_on_cell_1d : Not implemented yet\n"));
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 * 
 * Unset elementary functions
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_user_elementary_functions_unset (void)
{
  if (_user_fcts != NULL) {
    free (_user_fcts);
    _user_fcts = NULL;
  }
}


/*----------------------------------------------------------------------------
 * 
 * Set elementary functions
 * 
 * parameters:
 *   location_tetra    <-- Location in a tetrahedron
 *   location_prism    <-- Location in a prism
 *   location_pyramid  <-- Location in a pyramid
 *   location_hexa     <-- Location in a hexaedron
 *   location_tria     <-- Location on a triangle
 *   location_quad     <-- Location on a quandragle
 *   location_edge     <-- Location on a edge
 *   weight_tetra       <-- Weight computation in a tetrahedron
 *   weight_prism       <-- Weight computation in a prism
 *   weight_pyramid     <-- Weight computation in a pyramid
 *   weight_hexa        <-- Weight computation in a hexaedron
 *   weight_tria        <-- Weight computation on a triangle
 *   weight_quad        <-- Weight computation on a quandragle
 *   weight_edge        <-- Weight computation on a edge
 *   interp_tetra       <-- Interpolation in a tetrahedron
 *   interp_prism       <-- Interpolation in a prism
 *   interp_pyramid     <-- Interpolation in a pyramid
 *   interp_hexa        <-- Interpolation in a hexaedron
 *   interp_tria        <-- Interpolation on a triangle
 *   interp_quad        <-- Interpolation on a quandragle
 *   interp_edge        <-- Interpolation on a edge
 *
 *----------------------------------------------------------------------------*/

void
fvmc_ho_user_elementary_functions_set (fvmc_ho_location_fct_t location_tetra,
                                       fvmc_ho_location_fct_t location_prism,
                                       fvmc_ho_location_fct_t location_pyramid,
                                       fvmc_ho_location_fct_t location_hexa,
                                       fvmc_ho_location_fct_t location_tria,
                                       fvmc_ho_location_fct_t location_quad,
                                       fvmc_ho_location_fct_t location_edge,
                                       fvmc_ho_weight_fct_t weight_tetra,
                                       fvmc_ho_weight_fct_t weight_prism,
                                       fvmc_ho_weight_fct_t weight_pyramid,
                                       fvmc_ho_weight_fct_t weight_hexa,
                                       fvmc_ho_weight_fct_t weight_tria,
                                       fvmc_ho_weight_fct_t weight_quad,
                                       fvmc_ho_weight_fct_t weight_edge,
                                       fvmc_ho_interp_fct_t interp_tetra,
                                       fvmc_ho_interp_fct_t interp_prism,
                                       fvmc_ho_interp_fct_t interp_pyramid,
                                       fvmc_ho_interp_fct_t interp_hexa,
                                       fvmc_ho_interp_fct_t interp_tria,
                                       fvmc_ho_interp_fct_t interp_quad,
                                       fvmc_ho_interp_fct_t interp_edge)
{

  if (_user_fcts == NULL) {
    _user_fcts = (fvmc_ho_user_fcts_t *) malloc (sizeof(fvmc_ho_user_fcts_t));
  }

  _user_fcts->location_tetra   = location_tetra;
  _user_fcts->location_prism   = location_prism;
  _user_fcts->location_pyramid = location_pyramid;
  _user_fcts->location_hexa    = location_hexa;
  _user_fcts->location_tria    = location_tria;
  _user_fcts->location_quad    = location_quad;
  _user_fcts->location_edge    = location_edge;

  _user_fcts->weight_tetra   = weight_tetra;
  _user_fcts->weight_prism   = weight_prism;
  _user_fcts->weight_pyramid = weight_pyramid;
  _user_fcts->weight_hexa    = weight_hexa;
  _user_fcts->weight_tria    = weight_tria;
  _user_fcts->weight_quad    = weight_quad;
  _user_fcts->weight_edge    = weight_edge;

  _user_fcts->interp_tetra   = interp_tetra;
  _user_fcts->interp_prism   = interp_prism;
  _user_fcts->interp_pyramid = interp_pyramid;
  _user_fcts->interp_hexa    = interp_hexa;
  _user_fcts->interp_tria    = interp_tria;
  _user_fcts->interp_quad    = interp_quad;
  _user_fcts->interp_edge    = interp_edge;
}


/*----------------------------------------------------------------------------
 * 
 * Point location in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates if outside (or NULL)
 * 
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_in_cell_3d (const fvmc_element_t type,
                             const int order,
                             const int n_node,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{


  printf("-- fvmc_ho_location_in_cell_3d -- deb\n"); 

  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      return (_user_fcts->location_tetra) (order,
                                           n_node,
                                           ho_vertex_num,
                                           vertex_coords,
                                           point_coords,
                                           projected_coords);
      break;
      
    case FVMC_CELL_PRISM:
      return (_user_fcts->location_prism) (order,
                                           n_node,
                                           ho_vertex_num,
                                           vertex_coords,
                                           point_coords,
                                           projected_coords);
      break;
      
    case FVMC_CELL_PYRAM:
      return (_user_fcts->location_pyramid) (order,
                                             n_node,
                                             ho_vertex_num,
                                             vertex_coords,
                                             point_coords,
                                             projected_coords);
      break;

    case FVMC_CELL_HEXA:
      return (_user_fcts->location_hexa) (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;

    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_in_cell_3d : Not a high order 3D element type\n"));
    } 
  }

  else {

    return _default_location_in_cell_3d (type,
                                         order,
                                         n_node,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords);
    
  }

  printf("-- fvmc_ho_location_in_cell_3d -- fin\n"); 
}

/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_on_cell_2d (const fvmc_element_t type,
                             const int order,
                             const int n_node,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{

  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      return (_user_fcts->location_tria) (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;
      
    case FVMC_FACE_QUAD:
      return (_user_fcts->location_quad) (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_on_cell_2d : Not a high order 2D element type\n"));
    }
  }

  else {

    return _default_location_on_cell_2d (type,
                                         order,
                                         n_node,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords);
    
  }
}

/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_on_cell_1d (const fvmc_element_t type,
                             const int order,
                             const int n_node,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{

  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      return (_user_fcts->location_edge) (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_on_cell_1d : Not a high order 1D element type\n"));
    }
  }

  else {

    return _default_location_on_cell_1d (type,
                                         order,
                                         n_node,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords);

  }
  
}


/*----------------------------------------------------------------------------
 * 
 * Compute weight in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point inside cell  
 *   weight            --> weights
 * 
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_weight_in_cell_3d (const fvmc_element_t type,
                          const int order,
                          const int n_node,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *weight)
{
  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      (_user_fcts->weight_tetra) (order,
                                 n_node,
                                 ho_vertex_num,
                                 vertex_coords,
                                 point_coords,
                                 weight);
      break;
      
    case FVMC_CELL_PRISM:
      (_user_fcts->weight_prism) (order,
                                 n_node,
                                 ho_vertex_num,
                                 vertex_coords,
                                 point_coords,
                                 weight);
      break;
      
    case FVMC_CELL_PYRAM:
      (_user_fcts->weight_pyramid) (order,
                                   n_node,
                                   ho_vertex_num,
                                   vertex_coords,
                                   point_coords,
                                   weight);
      break;

    case FVMC_CELL_HEXA:
      (_user_fcts->weight_hexa) (order,
                                n_node,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                weight);
      break;

    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_weight_in_cell_3d : Not a high order 3D element type\n"));
    } 

  }

  else {

    _default_weight_in_cell_3d (type,
                               order,
                               n_node,
                               ho_vertex_num,
                               vertex_coords,
                               point_coords,
                               weight); 

  }
}

/*----------------------------------------------------------------------------
 * 
 * Compute weight on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point on cell  
 *   weight            --> weights
 * 
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_weight_on_cell_2d (const fvmc_element_t type,
                          const int order,
                          const int n_node,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *weight)
{
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      (_user_fcts->weight_tria) (order,
                                n_node,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                weight);
      break;
      
    case FVMC_FACE_QUAD:
      (_user_fcts->weight_quad) (order,
                                n_node,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                weight);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_weight_on_cell_2d : Not a high order 2D element type\n"));
    }
  }

  else {

    _default_weight_on_cell_2d (type,
                               order,
                               n_node,
                               ho_vertex_num,
                               vertex_coords,
                               point_coords,
                               weight); 

  }
}

/*----------------------------------------------------------------------------
 * 
 * Compute weight on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point on cell  
 *   weight            --> weights
 * 
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_weight_on_cell_1d (const fvmc_element_t type,
                          const int order,
                          const int n_node,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *weight)
{
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      (_user_fcts->weight_edge) (order,
                                n_node,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                weight);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_weight_on_cell_1d : Not a high order 1D element type\n"));
    }

  }

  else {

    _default_weight_on_cell_1d (type,
                               order,
                               n_node,
                               ho_vertex_num,
                               vertex_coords,
                               point_coords,
                               weight); 

  }
}


/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords in a 3D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_in_cell_3d (const fvmc_element_t type,
                           const int order,
                           const int n_node,
                           const int *ho_vertex_num,
                           const int *local_to_user,
                           const double *vertex_coords,
                           const double *point_coords,
                            const float *distance,
                           const double *point_proj_coords,
                          const double *weight,
                           const int stride_field,
                           const double *src_field,
                           double *target_field)
                           
{
  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      (_user_fcts->interp_tetra) (order,
                                  n_node,
                                  ho_vertex_num,
                                  local_to_user, 
                                  vertex_coords,
                                  point_coords,
                                  distance,
                                  point_proj_coords,
                                  weight,
                                  stride_field,
                                  src_field,
                                  target_field);

      break;
      
    case FVMC_CELL_PRISM:

      (_user_fcts->interp_prism) (order,
                                  n_node,
                                  ho_vertex_num,
                                  local_to_user, 
                                  vertex_coords,
                                  point_coords,
                                  distance,
                                  point_proj_coords,
                                  weight,
                                  stride_field,
                                  src_field,
                                  target_field);

      break;
      
    case FVMC_CELL_PYRAM:
      (_user_fcts->interp_pyramid ) (order,
                                     n_node,
                                     ho_vertex_num,
                                     local_to_user, 
                                     vertex_coords,
                                     point_coords,
                                     distance,
                                     point_proj_coords,
                                     weight,
                                     stride_field,
                                     src_field,
                                     target_field);

      break;

    case FVMC_CELL_HEXA:
      (_user_fcts->interp_hexa) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;

    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_interp_in_cell_3d : Not a high order 3D element type\n"));
    } 

  }

  else {

    _default_interp_in_cell_3d (type,
                                order,
                                n_node,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                distance,
                                point_proj_coords,
                                weight,
                                stride_field,
                                src_field,
                                target_field);

  }
}


/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_on_cell_2d (const fvmc_element_t type,
                           const int order,
                           const int n_node,
                           const int *ho_vertex_num,
                           const int *local_to_user,
                           const double *vertex_coords,
                           const double *point_coords,
                           const float *distance,
                           const double *point_proj_coords,
                           const double *weight,
                           const int stride_field,
                           const double *src_field,
                           double *target_field)
{
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      (_user_fcts->interp_tria) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;
      
    case FVMC_FACE_QUAD:
      (_user_fcts->interp_quad) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_interp_on_cell_2d : Not a high order 2D element type\n"));
    }

  }

  else {

    _default_interp_on_cell_2d (type,
                                order,
                                n_node,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                distance,
                                point_proj_coords,
                                weight,
                                stride_field,
                                src_field,
                                target_field);
  }
}


/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_on_cell_1d (const fvmc_element_t type,
                           const int order,
                           const int n_node,
                           const int *ho_vertex_num,
                           const int *local_to_user,
                           const double *vertex_coords,
                           const double *point_coords,
                           const float *distance,
                           const double *point_proj_coords,
                           const double *weight,
                           const int stride_field,
                           const double *src_field,
                           double *target_field)
{
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      (_user_fcts->interp_edge) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_interp_on_cell_1d : Not a high order 1D element type\n"));
    }

  }

  else {

    _default_interp_on_cell_1d (type,
                                order,
                                n_node,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                distance,
                                point_proj_coords,
                                weight,
                                stride_field,
                                src_field,
                                target_field);

  }
}
