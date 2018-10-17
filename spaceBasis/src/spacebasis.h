#ifndef __SPACEBASIS_H__
#define __SPACEBASIS_H__

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/**
 *
 * \brief Free an array allocated in fortran
 *
 * \param [in]     s_array  Size of array 
 * \param [inout]  array    Array
 * 
 */

void
SNB_free_double
(
 double    *array,
 const int  s_array
);


/**
 *
 * \brief Compute interpolation weights for nVtx (u,v) coordinates in a high order triangle
 *
 * \param [in]     order    Order 
 * \param [in]     nVtx     Number of vertices
 * \param [in]     ijk      Node location in ijk grid (size = 2 * nMod)
 * \param [in]     uv       uv coordinates (size = 2 * nVtx)
 * \param [in]     weights  Interpolation weights (size = nMod * nVtx)
 * 
 */

void
SNB_setT3BasisEqui_uv
(
 const int     order,
 const int     nVtx,
 const int    *ijk,
 const double *uv,
       double *weights
);


/**
 *
 * \brief Compute interpolation weights for nVtx (u,v) coordinates in a high order triangle
 *
 * \param [in]     order    Order 
 * \param [in]     nVtx     Number of vertices
 * \param [in]     ijk      Node location in ijk grid (size = 2 * nMod)
 * \param [in]     uv       u coordinates (size = nVtx)
 * \param [in]     v        v coordinates (size = nVtx)
 * \param [in]     ai       Interpolation weights (size = nMod * nVtx)
 * 
 */

void
SNB_setT3BasisEqui_u_v
(
 const int     order,
 const int     nVtx,
 const int    *ijk,
 const double *u,
 const double *v,
       double *weights
);


/**
 *
 * \brief Compute interpolation weights for nVtx (u,v) coordinates in a P3 triangle
 *
 * \param [in]     nVtx     Number of vertices
 * \param [in]     u        u coordinates (size = nVtx)
 * \param [in]     v        v coordinates (size = nVtx)
 * \param [in]     ai       Interpolation weights (size = nMod * nVtx)
 * 
 */

void
SNB_setT3Basis_u_v_P3
(
 const int    nVtx,
 const double *u,
 const double *v,
       double *weights
);


/**
 *
 * \brief Compute interpolation weights for nVtx (u,v) coordinates in a P3 triangle
 *
 * \param [in]     n_pts    Number of points
 * \param [in]     uv       uv coordinates (size = 2 * n_pts)
 * \param [in]     ai       Interpolation weights (size = nMod * nVtx)
 * 
 */

void
SNB_setT3Basis_P3
(
 const int    n_pts,
 const double *uv,
       double *weights
);


/**
 *
 * \brief Compute weights from u and v
 *
 * \param [in]     order    Order 
 * \param [in]     nVtx     Number of vertices 
 * \param [in]     uvw      uvw coordinates
 * \param [in]     display  Flag to dispay results
 * 
 */

void
SNB_setQ4BasisEqui_uv
(
 const int     order,
 const int     nVtx,
 const int    *ijk,
 const double *u,
 const double *v,
       double *weights
);


/**
 *
 * \brief Compute uvw of nodes for the current order
 *
 * \param [in]     order    Order 
 * \param [inout]  uvw      uvw coordinates
 * \param [in]     display  Flag to dispay results
 * 
 */

void
SNB_nodes2D
(
 const int     order,
       double **uvw,
 const int     display
);


/**
 *
 * \brief Compute ab coordinates from uv coordinates
 *
 * \param [in]     nVtx     Number of vertices 
 * \param [inout]  uv       uv coordinates
 * \param [out]    a        A coordinates 
 * \param [out]    b        B coordinates
 * \param [in]     display  Flag to dispay results
 * 
 */

void
SNB_nodes2Duv2ab
(
 const int     order,
 const double *uvw,
       double **a,
       double **b,
 const int     display
);


/**
 *
 * \brief Compute the Vandermonde matrix
 *
 * \param [in]     order    Order 
 * \param [inout]  a        A coordinates (size = n_nodes) 
 * \param [inout]  b        B coordinates (size = n_nodes)
 * \param [out]    vand     Vandermonde Matrix
 * 
 */

void
SNB_vandermonde2D
(
 const int     order,
 const double *a,
 const double *b,
       double **vand
);


/**
 *
 * \brief Compute lagrange polynome
 *
 * \param [in]      order     Order 
 * \param [in]      nVtx      Number of vertices
 * \param [in]      a         A coordinates (size = n_Vtx) 
 * \param [in]      b         B coordinates (size = n_Vtx)
 * \param [in]      vand      Vandermonde matrix (size = n_nodes X n_nodes)
 * \param [intout]  lx        Lagrange Coefficient (size = n_nVtx X n_nodes if transpose,
 *                                                  n_nodes X n_nVtx otherwise)
 * \param [in]      transpose Flag to transpose lx
 * 
 */

void
SNB_lagrange2Dv
(
 const int     order,
 const int     nVtx,
 const double *vand,
 const double *a,
 const double *b,
       double *lx
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SPACEBASIS_H__ */
