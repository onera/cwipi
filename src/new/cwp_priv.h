#ifndef __CWP_PRIV_H__
#define __CWP_PRIV_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-20  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


#include "cwp.h" 

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Suppress warning
 *============================================================================*/

#if defined(__INTEL_COMPILER)
#define CWP_PRAGMA_TO_STR(x) _Pragma(#x)
#define CWP_INTEL_SUPPRESS_WARNING_PUSH _Pragma("warning(push)")
#define CWP_INTEL_SUPPRESS_WARNING(w) CWP_PRAGMA_TO_STR(warning(disable:w))
#define CWP_INTEL_SUPPRESS_WARNING_POP _Pragma("warning(pop)")
#define CWP_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    CWP_INTEL_SUPPRESS_WARNING_PUSH CWP_INTEL_SUPPRESS_WARNING(w)
#else // CWP_INTEL
#define CWP_INTEL_SUPPRESS_WARNING_PUSH
#define CWP_INTEL_SUPPRESS_WARNING(w)
#define CWP_INTEL_SUPPRESS_WARNING_POP
#define CWP_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // CWP_INTEL

#if defined(__clang__)
#define CWP_PRAGMA_TO_STR(x) _Pragma(#x)
#define CWP_CLANG_SUPPRESS_WARNING_PUSH _Pragma("clang diagnostic push")
#define CWP_CLANG_SUPPRESS_WARNING(w) CWP_PRAGMA_TO_STR(clang diagnostic ignored w)
#define CWP_CLANG_SUPPRESS_WARNING_POP _Pragma("clang diagnostic pop")
#define CWP_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    CWP_CLANG_SUPPRESS_WARNING_PUSH CWP_CLANG_SUPPRESS_WARNING(w)
#else // CWP_CLANG
#define CWP_CLANG_SUPPRESS_WARNING_PUSH
#define CWP_CLANG_SUPPRESS_WARNING(w)
#define CWP_CLANG_SUPPRESS_WARNING_POP
#define CWP_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // CWP_CLANG

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#define CWP_PRAGMA_TO_STR(x) _Pragma(#x)
#define CWP_GCC_SUPPRESS_WARNING_PUSH _Pragma("GCC diagnostic push")
#define CWP_GCC_SUPPRESS_WARNING(w) CWP_PRAGMA_TO_STR(GCC diagnostic ignored w)
#define CWP_GCC_SUPPRESS_WARNING_POP _Pragma("GCC diagnostic pop")
#define CWP_GCC_SUPPRESS_WARNING_WITH_PUSH(w)                                                  \
    CWP_GCC_SUPPRESS_WARNING_PUSH CWP_GCC_SUPPRESS_WARNING(w)
#else // CWP_GCC
#define CWP_GCC_SUPPRESS_WARNING_PUSH
#define CWP_GCC_SUPPRESS_WARNING(w)
#define CWP_GCC_SUPPRESS_WARNING_POP
#define CWP_GCC_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // DOCTEST_GCC

#define CWP_UNUSED(x) (void)(x)  

CWP_GCC_SUPPRESS_WARNING("-Wcast-qual")
CWP_GCC_SUPPRESS_WARNING("-Wunknown-pragmas")

/*============================================================================
 * Type
 *============================================================================*/


/* TODO: CWP_Interpolation_t doit etre bascule dans un header prive (inutilise dans cwp.h) */

/**
 * \enum CWP_Interpolation_t
 * \brief Interpolation type
 *
 * CWP_Interpolation_t gives the different ways to interpolate
 */

typedef enum {

  CWP_INTERPOLATION_DEFAULT,  /*!< Default interpolation */
  CWP_INTERPOLATION_USER      /*!< User interpolation */

} CWP_Interpolation_t;

/*============================================================================
 * User interpolation type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public functions prototypes
 *============================================================================*/


/**
 *
 * \brief Setting of a FORTRAN user interpolation from location.
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_loc_set_f
(
 const char *local_code_name,
 const char *cpl_id,
 void       *fct
);


/**
 *
 * \brief Setting of a FORTRAN user interpolation from intersection.
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN .
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_inter_set_f
(
 const char *local_code_name,
 const char *cpl_id,
 void       *fct
);


/**
 *
 * \brief Setting of a FORTRAN user interpolation from closest points
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN .
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_closest_set_f
(
 const char *local_code_name,
 const char *cpl_id,
 void       *fct
);





void
CWP_surf_gen_init
(char* genName,
  int nx, int ny, int nPart, MPI_Comm* comm, double prop, double width, double randomVar
);

void
CWP_surf_gen_compute
(char* genName
);



void
CWP_surf_gen_by_block_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum, int* nElts,
  int* nTri , int** eltsConnecTri , CWP_g_num_t** eltsGnumTri,
  int* nQuad, int** eltsConnecQuad, CWP_g_num_t** eltsGnumQuad,
  int* nPoly, int** eltsConnecPolyIndex, int** eltsConnecPoly, CWP_g_num_t** eltsGnumPoly
);

void
CWP_surf_gen_one_connectivity_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum,
  int* nElts, int** eltsConnecIndex, int** eltsConnec, CWP_g_num_t** eltsGnum
);

void
CWP_surf_face_edge_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum,
  int* nFace, int** faceEdgeIdx, int** faceEdge,
  int* nEdge, int** edgeVtxIdx, int** edgeVtx,
  CWP_g_num_t** faceLNToGN
);


void
CWP_surf_gen_tri_field_get
( char* genName, int i_part,
  double** field
);

void
CWP_surf_gen_quad_field_get
( char* genName, int i_part,
  double** field
);

void
CWP_surf_gen_poly_field_get
( char* genName, int i_part,
  double** field
);



CWP_g_num_t*
CWP_GlobalNumGet
(
 const char  *local_code_name,
 const char  *cpl_id,
 const int    id_block,
 const int    i_part
);

MPI_Comm
CWP_Connectable_comm_get
(
 char* local_code_name
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWP_H__ */
