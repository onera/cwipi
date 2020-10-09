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


/*=============================================================================
 * Macro definitions
 *============================================================================*/

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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWP_H__ */
