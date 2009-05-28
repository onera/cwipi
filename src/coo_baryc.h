#ifndef __COO_BARYC_H__
#define __COO_BARYC_H__

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <fvm_defs.h>
#include <fvm_nodal.h>
#include <fvm_nodal_append.h>
#include <fvm_nodal_order.h>
#include <fvm_nodal_triangulate.h>
#include <fvm_writer.h>
#include <fvm_locator.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include <fvm_parall.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

void coo_baryc(const fvm_locator_t* locator, 
               const int   nMeshCoords, 
               const double *meshCoords,
               const int   nElts,                
               const int  *nMeshElts, 
               const int  *meshElts, 
               int  *n_dist_points, 
               int  **nDistBarCoords, 
               double **distBarCoords);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BAR_COORDS_H__ */
