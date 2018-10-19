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

#include "bftc_error.h"
#include "bftc_mem.h"
#include "bftc_printf.h"
#include "fvmc_point_location.h"
#include "cwipi_config.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_triangulate.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_ho_basis.h"

#if defined (HAVE_SPACE_BASIS) 
#include "spacebasis.h"
#endif

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

typedef struct _fvmc_ho_basis_user_elt_t {
  fvmc_ho_basis_fct_t elt_basis;
} fvmc_ho_basis_user_elt_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static int                 _n_ijk_tria_space = 0;
static int               **_ijk_tria_space = NULL;

static int                 _n_ijk_quad_space = 0;
static int               **_ijk_quad_space = NULL;

static int                 _n_ijk_1D_space = 0;
static int               **_ijk_1D_space = NULL;

static fvmc_ho_basis_user_elt_t *_user_edge = NULL;

static fvmc_ho_basis_user_elt_t *_user_tria = NULL;
static fvmc_ho_basis_user_elt_t *_user_quad = NULL;

static fvmc_ho_basis_user_elt_t *_user_tetra = NULL;
static fvmc_ho_basis_user_elt_t *_user_hexa = NULL;
static fvmc_ho_basis_user_elt_t *_user_prism = NULL;
static fvmc_ho_basis_user_elt_t *_user_pyra = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * 
 * Monomonial product
 * Procedure utilisee pour calculer (rapidement) les fonctions de base des simplex
 * dont les noeuds d'interpolation sont equidistants.
 * 
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   u               <-- u (size =  n_pts)
 *   fn              --> fn (size =  n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_monomialProduct
(
 const int order,
 const int n_pts,
 const double *u,
       double *fn
)
{  

  double constant;

  for (int i = 0; i < n_pts; i++) {
    fn[i]= 1.;
  }

  for (int i = 0; i < n_pts; i++) {

    constant = 1. / (double) (i - n_pts);
    
    for (int j = 0; i < n_pts; i++) {
      fn[j] *= ((double) i - (double) order * u[j]) * constant;
    }
  }
}

/*----------------------------------------------------------------------------
 * 
 * Triangle Pn basis
 * 
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uv              <-- u (size = 2 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_tria_pn
(
 const int order,
 const int n_pts,
 const double *uv,
 double *weights
)
{

  if (order == 1) {

    if (n_pts != 1) {

      double *u  = malloc (sizeof(double) *n_pts);
      double *v  = malloc (sizeof(double) *n_pts);
      
      for (int i = 0; i < n_pts; i++) {
        u[i]  = uv[2*i];
        v[i]  = uv[2*i+1];
      }
      
      for (int i = 0; i < n_pts; i++) {
        weights[3*i]   = 1. - u[i] - v[i];
        weights[3*i+1] = u[i];
        weights[3*i+2] = v[i];
      }
      
      free (u);
      free (v);
    }

    else {

      double u = uv[0];
      double v = uv[1];
      
      weights[0] = 1. - u - v;
      weights[1] = u;
      weights[2] = v;

    }
  }
  
  else if (order == 2) {

    if (n_pts != 1) {
      double *u  = malloc (sizeof(double) *n_pts);
      double *v  = malloc (sizeof(double) *n_pts);
      double *w  = malloc (sizeof(double) *n_pts);
      double *u2 = malloc (sizeof(double) *n_pts);
      double *v2 = malloc (sizeof(double) *n_pts);
      double *w2 = malloc (sizeof(double) *n_pts);
    

      for (int i = 0; i < n_pts; i++) {
        u[i]  = uv[2*i];
        v[i]  = uv[2*i+1];
        w[i]  = 1. - u[i] - v[i];
        u2[i] = 2. * u[i];
        v2[i] = 2. * v[i];
        w2[i] = 2. * w[i];
      }
      
      for (int i = 0; i < n_pts; i++) {
        
        weights[6*i+0] = w[i] * (-1. + w2[i]);  /* (i,j,k)=(0,0,2) */
        weights[6*i+1] = u2[i] * w2[i];         /* (i,j,k)=(1,0,1) */
        weights[6*i+2] = u[i] * (-1. + u2[i]);  /* (i,j,k)=(2,0,0) */
        weights[6*i+3] = v2[i] * w2[i];         /* (i,j,k)=(0,1,1) */
        weights[6*i+4] = u2[i] * v2[i];         /* (i,j,k)=(1,1,0) */
        weights[6*i+5] = v[i] * (-1. + v2[i]);  /* (i,j,k)=(0,2,0) */
      }
    
      free (u);
      free (v);
      free (w);
      free (u2);
      free (v2);
      free (w2);
    }
    else {

      double u = uv[0];
      double v = uv[1];
      
      double w  = 1. - u - v;
      double u2 = 2. * u;
      double v2 = 2. * v;
      double w2 = 2. * w;
      
      weights[0] = w * (-1. + w2);  /* (i,j,k)=(0,0,2) */
      weights[1] = u2 * w2;         /* (i,j,k)=(1,0,1) */
      weights[2] = u * (-1. + u2);  /* (i,j,k)=(2,0,0) */
      weights[3] = v2 * w2;         /* (i,j,k)=(0,1,1) */
      weights[4] = u2 * v2;         /* (i,j,k)=(1,1,0) */
      weights[5] = v * (-1. + v2);  /* (i,j,k)=(0,2,0) */

    }
  }

  else if (order == 3) {

    for (int i = 0;  i < n_pts; i++) {
      
      double u = uv[2*i];
      double u3 = 3.*u;
      double u3m1 = (u3-1.)*5e-1;
      
      double v = uv[2*i+1];
      double v3 = 3.*v;
      double v3m1 = (v3-1.)*5e-1;
      
      double w = 1. - u - v;
      double w3 = 3.*w;
      double w3m1 = (w3-1.)*5e-1;
      
      weights[2*i+0] = w*w3m1*(w3-2.);   // (i,j,k)=(003)
      weights[2*i+3] = u*u3m1*(u3-2.);   // (i,j,k)=(300)
      weights[2*i+9] = v*v3m1*(v3-2.);   // (i,j,k)=(030)
      weights[2*i+1] = u3*w3*w3m1;       // (i,j,k)=(102)
      
      double coef = u3*u3m1;
      weights[2*i+2] = coef*w3;           //(i,j,k)=(201)
      weights[2*i+6] = coef*v3;           // (i,j,k)=(210)
      
      coef=v3*v3m1;
      weights[2*i+8] = coef*u3;           // (i,j,k)=(120)
      weights[2*i+7] = coef*w3;           // (i,j,k)=(021)
      
      coef=v3*w3;
      weights[2*i+4] = coef*w3m1;         // (i,j,k)=(012)
      weights[2*i+5] = coef*u3;           // (i,j,k)=(111)
      
    }
  }
  
  else {

    double *u  = malloc (sizeof(double) *n_pts);
    double *v  = malloc (sizeof(double) *n_pts);
    double *w  = malloc (sizeof(double) *n_pts);
    double *fu = malloc (sizeof(double) *n_pts);
    double *fv = malloc (sizeof(double) *n_pts);
    double *fw = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = uv[2*i];
      v[i] = uv[2*i+1];
      w[i] = 1. - u[i] - v[i];
    }
    
    //print '("size(ai)=",i2,"x",i2)',size(ai,1),size(ai,2)

    /* for (int j = 0; j < n_nodes; j++) { */


    /*   // Attention monomonial product a revoir */

    /*   do iMod=1,nMod */
    /*   iu=ijk(1,iMod) */
    /*   iv=ijk(2,iMod) */
    /*   iw=ord-iu-iv */
      
    /*   call monomialProduct(ord=ord,n=iu,u=u, fn=fu) */
    /*   call monomialProduct(ord=ord,n=iv,u=v, fn=fv) */
    /*   call monomialProduct(ord=ord,n=iw,u=w, fn=fw) */
      
    /*   !print '(3x,"iMod=",i3,3x,"iu=",i3," iv=",i3," iw=",i3)',iMod,iu,iv,iw */
    /*   do iNod=1,nNod */
    /*     ai(iMod,iNod)=fu(iNod)*fv(iNod)*fw(iNod) */
    /*   enddo */
      
    /* enddo */
    
    /* deallocate(u,v,w,fu,fv,fw) */
    /* !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    



    
#if defined (HAVE_SPACE_BASIS)
    
    if (_ijk_tria_space == NULL) {
      _n_ijk_tria_space = FVMC_MAX (order-1, 10);
      _ijk_tria_space = malloc (sizeof(int *) * _n_ijk_tria_space);
      for (int i = 0; i < _n_ijk_tria_space; i++) {
        _ijk_tria_space[i] = NULL;
      }
    }

    if (order > _n_ijk_tria_space) {
      _ijk_tria_space = realloc (_ijk_tria_space, sizeof(int *) * _n_ijk_tria_space);
      for (int i = _n_ijk_tria_space; i < order; i++) {
        _ijk_tria_space[i] = NULL;
      }
    }
    
    int *__ijk_tria_space = _ijk_tria_space[order-1];

    if (__ijk_tria_space == NULL) {
      const int n_nodes = (order+1)*(order+2)/2;
      __ijk_tria_space = malloc (sizeof(int) * 2 * n_nodes);
      int k = 0;
      for (int j = 0; j < order+1; j++) {
        for (int i = 0; i < order+1-j; i++) {
          __ijk_tria_space[2*k]   = i;
          __ijk_tria_space[2*k+1] = j;
          k += 1;
        }
      }
    }

    SNB_setT3BasisEqui_uv (order, n_pts, __ijk_tria_space, uv, weights);

#else
    bftc_error(__FILE__, __LINE__, 0,
               _("_basis_tria_pn not yet implemented for order > 2 without space basis \n"));
#endif

  }
}




/*----------------------------------------------------------------------------
 * 
 * Quadrangle Qn basis
 * 
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uv              <-- u (size = 2 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/


static void
_basis_quad_qn
(
 const int order,
 const int n_pts,
 const double *uv,
 double *weights
)
{
  
  if (order == 1) {

    double *u = malloc (sizeof(double) *n_pts);
    double *v = malloc (sizeof(double) *n_pts);
    double *u1 = malloc (sizeof(double) *n_pts);
    double *v1 = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = uv[2*i];
      v[i] = uv[2*i+1];
      u1[i] = (1 - u[i]);
      v1[i] = (1 - v[i]);
    }
    
    for (int i = 0; i < n_pts; i++) {
      weights[4*i+0] = u1[i] * v1[i];
      weights[4*i+1] = u[i] * v1[i];
      weights[4*i+2] = u1[i] * v[i];
      weights[4*i+3] = u[i] * v[i];
    }

    free(u1);
    free(v1);
    free(u);
    free(v);

  }

  else if (order == 2) {
    
    double *u = malloc (sizeof(double) *n_pts);
    double *v = malloc (sizeof(double) *n_pts);

    double *uM = malloc (sizeof(double) *n_pts);
    double *uP = malloc (sizeof(double) *n_pts);
    double *u0 = malloc (sizeof(double) *n_pts);

    double *au1 = malloc (sizeof(double) *n_pts);
    double *au2 = malloc (sizeof(double) *n_pts);
    double *au3 = malloc (sizeof(double) *n_pts);
    
    double *vM = malloc (sizeof(double) *n_pts);
    double *vP = malloc (sizeof(double) *n_pts);
    double *v0 = malloc (sizeof(double) *n_pts);

    double *av1 = malloc (sizeof(double) *n_pts);
    double *av2 = malloc (sizeof(double) *n_pts);
    double *av3 = malloc (sizeof(double) *n_pts);

    
    for (int i = 0; i < n_pts; i++) {
      u[i] = uv[2*i];
      v[i] = uv[2*i+1];

      uM[i] = 2*(1-u[i]);
      uP[i] = 2*u[i];
      u0[i] = u[i]-0.5;
      
      au1[i] = -uM[i] * u0[i]; 
      au2[i] =  uM[i] * uP[i];
      au3[i] =  u0[i] * uP[i];
    
      vM[i] = 2*(1-v[i]);
      vP[i] = 2*v[i];
      v0[i] = v[i]-0.5;

      av1[i] = -vM[i] * v0[i]; 
      av2[i] =  vM[i] * vP[i];
      av3[i] =  v0[i] * vP[i];
    }
    
    for (int i = 0; i < n_pts; i++) {
      weights[9*i+0]=au1[i]*av1[i];
      weights[9*i+1]=au2[i]*av1[i];
      weights[9*i+2]=au3[i]*av1[i];
      weights[9*i+3]=au1[i]*av2[i];
      weights[9*i+4]=au2[i]*av2[i];
      weights[9*i+5]=au3[i]*av2[i];
      weights[9*i+6]=au1[i]*av3[i];
      weights[9*i+7]=au2[i]*av3[i];
      weights[9*i+8]=au3[i]*av3[i];
    }
    
    free(u);
    free(v);

    free(uM);
    free(uP);
    free(u0);

    free(au1);
    free(au2);
    free(au3);
    
    free(vM);
    free(vP);
    free(v0);

    free(av1);
    free(av2);
    free(av3);

  }

  else {

#if defined (HAVE_SPACE_BASIS)
    
    if (_ijk_quad_space == NULL) {
      _n_ijk_quad_space = FVMC_MAX (order-1, 10);
      _ijk_quad_space = malloc (sizeof(int *) * _n_ijk_quad_space);
      for (int i = 0; i < _n_ijk_quad_space; i++) {
        _ijk_quad_space[i] = NULL;
      }
    }

    if (order > _n_ijk_quad_space) {
      _ijk_quad_space = realloc (_ijk_quad_space, sizeof(int *) * _n_ijk_quad_space);
      for (int i = _n_ijk_quad_space; i < order; i++) {
        _ijk_quad_space[i] = NULL;
      }
    }
    
    int *__ijk_quad_space = _ijk_quad_space[order-1];

    if (__ijk_quad_space == NULL) {
      const int n_nodes = (order+1)*(order+1);
      __ijk_quad_space = malloc (sizeof(int) * 2 * n_nodes);
      int k = 0;
      for (int j = 0; j < order+1; j++) {
        for (int i = 0; i < order+1; i++) {
          __ijk_quad_space[k++] = i;
          __ijk_quad_space[k++] = j;
        }
      }
    }

    double *u = malloc (sizeof(double) *n_pts);
    double *v = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = 2 * uv[2*i]   - 1;
      v[i] = 2 * uv[2*i+1] - 1;
    }
    
    SNB_setQ4BasisEqui_uv (order, n_pts, __ijk_quad_space,
                           u, v,
                           weights);

    free (u);
    free (v);
    
#else
    bftc_error(__FILE__, __LINE__, 0,
               _("_basis_quad_pn not yet implemented for order > 2 without space basis \n"));
#endif
  }
}
  



/*----------------------------------------------------------------------------
 * 
 * high order basis
 * 
 * parameters:
 *   type            <-- element type
 *
 * return:
 *
 *----------------------------------------------------------------------------*/

static fvmc_ho_basis_user_elt_t *
_get_user_elt (fvmc_element_t elt_type)
{

  fvmc_ho_basis_user_elt_t *user_elt = NULL;
  
  switch(elt_type) {

  case FVMC_EDGE:
    user_elt = _user_edge;
    break;
    
  case FVMC_FACE_TRIA:
    user_elt = _user_tria;
    break;
    
  case FVMC_FACE_QUAD:
    user_elt = _user_quad;
    break;
    
  case FVMC_CELL_TETRA:
    user_elt = _user_tetra;
    break;
    
  case FVMC_CELL_PYRAM:
    user_elt = _user_pyra;
    break;
    
  case FVMC_CELL_PRISM:
    user_elt = _user_prism;
    break;
    
  case FVMC_CELL_HEXA:
    user_elt = _user_hexa;
    break;

  default:
    bftc_error(__FILE__, __LINE__, 0,
               _("fvmc_ho_user_elt_unset : Unvailable element type\n"));
  }

  return user_elt;
}

/*----------------------------------------------------------------------------
 * 
 * default high order basis
 * 
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_pts           <-- number of points 
 *   uvw             <-- uvw (size = elt_dim * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *----------------------------------------------------------------------------*/

static void
_default_elt_basis
(
const fvmc_element_t type,
const int order,
const int n_pts,
const double *uvw,
      double *weights 
)
{
  switch (type) {

  case FVMC_FACE_TRIA:
    _basis_tria_pn (order, n_pts, uvw, weights);
    break;

  case FVMC_FACE_QUAD:
    _basis_quad_qn (order, n_pts, uvw, weights);
    break;

  case FVMC_EDGE:
  case FVMC_CELL_TETRA:
  case FVMC_CELL_PRISM:
  case FVMC_CELL_PYRAM:
  case FVMC_CELL_HEXA:
  default: 
    bftc_error(__FILE__, __LINE__, 0,
               _("_default_elts_basis : '%d' element type not yet implemented\n"),
               type);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * 
 * Unset a user element
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_basis_user_elt_unset (fvmc_element_t elt_type)
{

  fvmc_ho_basis_user_elt_t *_user_elt = _get_user_elt (elt_type);
  
  if (_user_elt != NULL) {
    free (_user_elt);
    _user_elt = NULL;
  }
  
}


void
fvmc_ho_basis_user_elts_unset (void)
{
  fvmc_ho_basis_user_elt_unset (FVMC_EDGE);
  fvmc_ho_basis_user_elt_unset (FVMC_FACE_TRIA);
  fvmc_ho_basis_user_elt_unset (FVMC_FACE_QUAD);
  fvmc_ho_basis_user_elt_unset (FVMC_CELL_TETRA);
  fvmc_ho_basis_user_elt_unset (FVMC_CELL_HEXA);
  fvmc_ho_basis_user_elt_unset (FVMC_CELL_PRISM);
  fvmc_ho_basis_user_elt_unset (FVMC_CELL_PYRAM);
}

/*----------------------------------------------------------------------------
 * 
 * Unset a user element
 * 
 *----------------------------------------------------------------------------*/


void
fvmc_ho_basis_user_elt_set (fvmc_element_t elt_type,
                            fvmc_ho_basis_fct_t elt_basis)
{
  fvmc_ho_basis_user_elt_t *user_elt = _get_user_elt (elt_type);
  
  if (user_elt == NULL) {
    user_elt = (fvmc_ho_basis_user_elt_t *) malloc (sizeof(fvmc_ho_basis_user_elt_t));
  }

  user_elt->elt_basis = elt_basis;

}


/*----------------------------------------------------------------------------
 * 
 * high order basis
 * 
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_nodes         <-- number of nodes
 *   n_pts           <-- number of points 
 *   uvw             <-- uvw
 *   uvw             <-- uvw (size = n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *----------------------------------------------------------------------------*/

void
fvmc_ho_basis
(
const fvmc_element_t type,
const int order,
const int n_nodes,
const int n_pts,
const double *uvw,
      double *weights 
)
{
  fvmc_ho_basis_user_elt_t *user_elt = _get_user_elt (type);
  
  if (user_elt != NULL) {
    if (user_elt->elt_basis != NULL) {
      (user_elt->elt_basis) (order,
                             n_nodes,
                             n_pts,
                             uvw,
                             weights);
    }
    else {
      _default_elt_basis (type,
                          order,
                          n_pts,
                          uvw,
                          weights);
    }
  }

  else {
    _default_elt_basis (type,
                        order,
                        n_pts,
                        uvw,
                        weights);

  }
}

/*----------------------------------------------------------------------------
 * 
 * Free static variables
 * 
 *----------------------------------------------------------------------------*/


void
fvmc_ho_basis_free
(
 void
)
{
#if defined (HAVE_SPACE_BASIS) 
  if (_ijk_tria_space != NULL) {

    for (int i = 0; i < _n_ijk_tria_space; i++) {
      if (_ijk_tria_space[i] != NULL) {
        free(_ijk_tria_space[i]);
        _ijk_tria_space[i] = NULL;
      }
    }
    free (_ijk_tria_space);
    _ijk_tria_space = NULL;
  
  }

  if (_ijk_quad_space != NULL) {

    for (int i = 0; i < _n_ijk_quad_space; i++) {
      if (_ijk_quad_space[i] != NULL) {
        free(_ijk_quad_space[i]);
        _ijk_quad_space[i] = NULL;
      }
    }
    free (_ijk_quad_space);
    _ijk_quad_space = NULL;
    
  }
  
  if (_ijk_1D_space != NULL) {

    for (int i = 0; i < _n_ijk_1D_space; i++) {
      if (_ijk_1D_space[i] != NULL) {
        free(_ijk_1D_space[i]);
        _ijk_1D_space[i] = NULL;
      }
    }
    
    free (_ijk_1D_space);
    _ijk_1D_space = NULL;
  }

#endif
  
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
