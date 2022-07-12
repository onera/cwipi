/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"

#include "pdm_part_domain_interface.h"
#include "pdm_part_domain_interface_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/


/*============================================================================
 * Public function definitions
 *============================================================================*/
PDM_part_domain_interface_t*
PDM_part_domain_interface_create
(
const int                    n_interface,
const int                    n_domain,
const int                   *n_part,
PDM_domain_interface_mult_t  multidomain_interface,
PDM_ownership_t              ownership,
PDM_MPI_Comm                 comm
)
{
  PDM_part_domain_interface_t *dom_intrf = (PDM_part_domain_interface_t *) malloc (sizeof(PDM_part_domain_interface_t));
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_domain          = n_domain;
  dom_intrf->n_part            = malloc(n_domain * sizeof(int));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    dom_intrf->n_part[i_domain] = n_part[i_domain];
  }
  dom_intrf->multidomain_intrf = multidomain_interface;
  dom_intrf->ownership         = ownership;
  dom_intrf->comm              = comm;

  dom_intrf->interface_pn_vtx        = (int          ***) malloc(n_domain * sizeof(int          **));
  dom_intrf->interface_vtx_ln_to_gn  = (PDM_g_num_t ****) malloc(n_domain * sizeof(PDM_g_num_t ***));
  dom_intrf->interface_sgn_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_vtx_idx   = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_dom_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));

  dom_intrf->interface_pn_face       = (int          ***) malloc(n_domain * sizeof(int          **));
  dom_intrf->interface_face_ln_to_gn = (PDM_g_num_t ****) malloc(n_domain * sizeof(PDM_g_num_t ***));
  dom_intrf->interface_sgn_face      = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_face      = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_face_idx  = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_dom_face      = (int         ****) malloc(n_domain * sizeof(int         ***));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    dom_intrf->interface_pn_vtx       [i_domain] = (int          **) malloc(n_part[i_domain] * sizeof(int          *));
    dom_intrf->interface_vtx_ln_to_gn [i_domain] = (PDM_g_num_t ***) malloc(n_part[i_domain] * sizeof(PDM_g_num_t **));
    dom_intrf->interface_sgn_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_vtx_idx  [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_dom_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));

    dom_intrf->interface_pn_face      [i_domain] = (int          **) malloc(n_part[i_domain] * sizeof(int          *));
    dom_intrf->interface_face_ln_to_gn[i_domain] = (PDM_g_num_t ***) malloc(n_part[i_domain] * sizeof(PDM_g_num_t **));
    dom_intrf->interface_sgn_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_face_idx [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_dom_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));

    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      dom_intrf->interface_pn_vtx       [i_domain][i_part] = NULL;
      dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] = NULL;
      dom_intrf->interface_sgn_vtx      [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_vtx      [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] = NULL;
      dom_intrf->interface_dom_vtx      [i_domain][i_part] = NULL;
      dom_intrf->interface_pn_face      [i_domain][i_part] = NULL;
      dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = NULL;
      dom_intrf->interface_sgn_face     [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_face     [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_face_idx [i_domain][i_part] = NULL;
      dom_intrf->interface_dom_face     [i_domain][i_part] = NULL;
    }
  }

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    dom_intrf->is_result[i] = 0;
  }

  dom_intrf->translation_vect   = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_direction = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_center    = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_angle     = (double  *) malloc(n_interface * sizeof(double  ));

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    dom_intrf->translation_vect  [i_interface] = NULL;
    dom_intrf->rotation_direction[i_interface] = NULL;
    dom_intrf->rotation_center   [i_interface] = NULL;
    dom_intrf->rotation_angle    [i_interface] = 0.;
  }

  return dom_intrf;
}

void
PDM_part_domain_interface_set
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind,
 int                           i_domain,
 int                           i_part,
 int                          *interface_pn,
 PDM_g_num_t                 **interface_ln_to_gn,
 int                         **interface_sgn,
 int                         **interface_ids,
 int                         **interface_ids_idx,
 int                         **interface_dom
)
{
  assert(i_part < dom_intrf->n_part[i_domain]);

  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    dom_intrf->interface_pn_vtx      [i_domain][i_part] = interface_pn;
    dom_intrf->interface_vtx_ln_to_gn[i_domain][i_part] = interface_ln_to_gn;
    dom_intrf->interface_sgn_vtx     [i_domain][i_part] = interface_sgn;
    dom_intrf->interface_ids_vtx     [i_domain][i_part] = interface_ids;
    dom_intrf->interface_ids_vtx_idx [i_domain][i_part] = interface_ids_idx;
    dom_intrf->interface_dom_vtx     [i_domain][i_part] = interface_dom;
  }
  else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    dom_intrf->interface_pn_face      [i_domain][i_part] = interface_pn;
    dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = interface_ln_to_gn;
    dom_intrf->interface_sgn_face     [i_domain][i_part] = interface_sgn;
    dom_intrf->interface_ids_face     [i_domain][i_part] = interface_ids;
    dom_intrf->interface_ids_face_idx [i_domain][i_part] = interface_ids_idx;
    dom_intrf->interface_dom_face     [i_domain][i_part] = interface_dom;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }
}

void
PDM_part_domain_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind,
 int                            i_domain,
 int                            i_part,
 int                          **interface_pn,
 PDM_g_num_t                 ***interface_ln_to_gn,
 int                         ***interface_sgn,
 int                         ***interface_ids,
 int                         ***interface_ids_idx,
 int                         ***interface_dom
)
{
  assert(i_part < dom_intrf->n_part[i_domain]);

  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    *interface_pn       = dom_intrf->interface_pn_vtx      [i_domain][i_part];
    *interface_ln_to_gn = dom_intrf->interface_vtx_ln_to_gn[i_domain][i_part];
    *interface_sgn      = dom_intrf->interface_sgn_vtx     [i_domain][i_part];
    *interface_ids      = dom_intrf->interface_ids_vtx     [i_domain][i_part];
    *interface_ids_idx  = dom_intrf->interface_ids_vtx_idx [i_domain][i_part];
    *interface_dom      = dom_intrf->interface_dom_vtx     [i_domain][i_part];
  }
  else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    *interface_pn       = dom_intrf->interface_pn_face      [i_domain][i_part];
    *interface_ln_to_gn = dom_intrf->interface_face_ln_to_gn[i_domain][i_part];
    *interface_sgn      = dom_intrf->interface_sgn_face     [i_domain][i_part];
    *interface_ids      = dom_intrf->interface_ids_face     [i_domain][i_part];
    *interface_ids_idx  = dom_intrf->interface_ids_face_idx [i_domain][i_part];
    *interface_dom      = dom_intrf->interface_dom_face     [i_domain][i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }
}

int
PDM_part_domain_interface_n_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf
)
{
  return dom_intrf->n_interface;
}

void
PDM_part_domain_interface_free
(
 PDM_part_domain_interface_t  *dom_intrf
)
{


  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain ) {

    for(int i_part = 0; i_part < dom_intrf->n_part[i_domain]; ++i_part) {

      if(dom_intrf->ownership == PDM_OWNERSHIP_KEEP) {
        if(dom_intrf->interface_pn_vtx       [i_domain][i_part] != NULL) {
          free(dom_intrf->interface_pn_vtx       [i_domain][i_part]);
          dom_intrf->interface_pn_vtx       [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part]);
          dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_sgn_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_sgn_vtx      [i_domain][i_part]);
          dom_intrf->interface_sgn_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_vtx      [i_domain][i_part]);
          dom_intrf->interface_ids_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part]);
          dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_dom_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_dom_vtx      [i_domain][i_part]);
          dom_intrf->interface_dom_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_pn_face      [i_domain][i_part] != NULL) {
          free(dom_intrf->interface_pn_face      [i_domain][i_part]);
          dom_intrf->interface_pn_face      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_face_ln_to_gn[i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part]);
          dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_sgn_face     [i_domain][i_part]);
          dom_intrf->interface_sgn_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_face     [i_domain][i_part]);
          dom_intrf->interface_ids_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face_idx [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_face_idx [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_face_idx [i_domain][i_part]);
          dom_intrf->interface_ids_face_idx [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_dom_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_dom_face     [i_domain][i_part]);
          dom_intrf->interface_dom_face     [i_domain][i_part] = NULL;
        };
      }
    }

    free(dom_intrf->interface_pn_vtx       [i_domain]);
    free(dom_intrf->interface_vtx_ln_to_gn [i_domain]);
    free(dom_intrf->interface_sgn_vtx      [i_domain]);
    free(dom_intrf->interface_ids_vtx      [i_domain]);
    free(dom_intrf->interface_ids_vtx_idx  [i_domain]);
    free(dom_intrf->interface_dom_vtx      [i_domain]);

    free(dom_intrf->interface_pn_face      [i_domain]);
    free(dom_intrf->interface_face_ln_to_gn[i_domain]);
    free(dom_intrf->interface_sgn_face     [i_domain]);
    free(dom_intrf->interface_ids_face     [i_domain]);
    free(dom_intrf->interface_ids_face_idx [i_domain]);
    free(dom_intrf->interface_dom_face     [i_domain]);

  }
  free(dom_intrf->n_part);

  free(dom_intrf->interface_pn_vtx       );
  free(dom_intrf->interface_vtx_ln_to_gn );
  free(dom_intrf->interface_sgn_vtx      );
  free(dom_intrf->interface_ids_vtx      );
  free(dom_intrf->interface_ids_vtx_idx  );
  free(dom_intrf->interface_dom_vtx      );

  free(dom_intrf->interface_pn_face      );
  free(dom_intrf->interface_face_ln_to_gn);
  free(dom_intrf->interface_sgn_face     );
  free(dom_intrf->interface_ids_face     );
  free(dom_intrf->interface_ids_face_idx );
  free(dom_intrf->interface_dom_face     );

  for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface) {
    if(dom_intrf->translation_vect[i_interface]   != NULL) {
      free(dom_intrf->translation_vect[i_interface]);
      dom_intrf->translation_vect[i_interface] = NULL;
    }
    if(dom_intrf->rotation_direction[i_interface]   != NULL) {
      free(dom_intrf->rotation_direction[i_interface]);
      dom_intrf->rotation_direction[i_interface] = NULL;
    }
    if(dom_intrf->rotation_center[i_interface]   != NULL) {
      free(dom_intrf->rotation_center[i_interface]);
      dom_intrf->rotation_center[i_interface] = NULL;
    }
  }

  free(dom_intrf->translation_vect  );
  free(dom_intrf->rotation_direction);
  free(dom_intrf->rotation_center   );
  free(dom_intrf->rotation_angle    );

  free(dom_intrf);
}




void
PDM_part_domain_interface_translation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
  const double                       *vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->translation_vect[i_interface] == NULL);

  dom_intrf->translation_vect[i_interface] = (double *) malloc( 3 * sizeof(double));

  for(int i = 0; i < 3; ++i) {
    dom_intrf->translation_vect[i_interface][i] = vect[i];
  }

}

void
PDM_part_domain_interface_rotation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
  const double                       *direction,
  const double                       *center,
  const double                        angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->rotation_direction[i_interface] == NULL);
  assert(dom_intrf->rotation_center   [i_interface] == NULL);

  dom_intrf->rotation_direction[i_interface] = (double *) malloc( 3 * sizeof(double));
  dom_intrf->rotation_center   [i_interface] = (double *) malloc( 3 * sizeof(double));

  for(int i = 0; i < 3; ++i) {
    dom_intrf->rotation_direction[i_interface][i] = direction[i];
    dom_intrf->rotation_center   [i_interface][i] = center   [i];
  }
  dom_intrf->rotation_angle[i_interface] = angle;
}



void
PDM_part_domain_interface_translation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
        double                      **vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->translation_vect[i_interface] != NULL){

    *vect = (double *) malloc( 3 * sizeof(double));
    double* _vect = *vect;

    for(int i = 0; i < 3; ++i) {
      _vect[i] = dom_intrf->translation_vect[i_interface][i];
    }
  } else {
    *vect = NULL;
  }
}

void
PDM_part_domain_interface_rotation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
        double                      **direction,
        double                      **center,
        double                       *angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->rotation_direction[i_interface] != NULL) {
    assert(dom_intrf->rotation_center   [i_interface] != NULL);

    *direction = (double *) malloc( 3 * sizeof(double));
    *center    = (double *) malloc( 3 * sizeof(double));
    double *_direction = *direction;
    double *_center    = *center   ;

    for(int i = 0; i < 3; ++i) {
      _direction[i] = dom_intrf->rotation_direction[i_interface][i];
      _center   [i] = dom_intrf->rotation_center   [i_interface][i];
    }
    *angle = dom_intrf->rotation_angle[i_interface];
  } else {
    *direction = NULL;
    *center    = NULL;
    *angle     = 0;
  }
}
