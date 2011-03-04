#ifndef __COO_BARYC_H__
#define __COO_BARYC_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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
namespace cwipi {

void coo_baryc(const fvm::fvm_locator_t* locator,
               const int   nMeshCoords,
               const double *meshCoords,
               const int   nElts,
               const int  *nMeshElts,
               const int  *meshElts,
               int  *n_dist_points,
               int  **nDistBarCoords,
               double **distBarCoords);

}

#endif /* __BAR_COORDS_H__ */
