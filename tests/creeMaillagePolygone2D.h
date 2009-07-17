#ifndef __CREE_MAIL_POLY_H__
#define __CREE_MAIL_POLY_H__

#include <fvm_defs.h>
#include <fvm_nodal.h>
#include <mpi.h>

#include "couplings_config.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

void creeMaillagePolygone2D(int order, 
                            MPI_Comm localComm,
                            double xmin, 
                            double xmax, 
                            double ymin, 
                            double ymax,
                            int initRandom,
                            int nx,
                            int ny,
                            int *nVertex,
                            double **meshCoords,
                            int *nElts,
                            int **eltsConnecPointer,
                            int **eltsConnec
                            );

void PROCF(creeMaillagePolygone2D_f, CREEMAILLAGEPOLYGONE2D_F)(int *order,
                                                               MPI_Comm *localComm,
                                                               double *xmin, 
                                                               double *xmax, 
                                                               double *ymin, 
                                                               double *ymax,
                                                               int *initRandom,
                                                               int *nx,
                                                               int *ny,
                                                               int *nVertex,
                                                               double *meshCoords_f,
                                                               int *nElts,
                                                               int *lEltsConnecPointer_f,
                                                               int *eltsConnecPointer_f,
                                                               int *eltsConnec_f
                                                               );
#endif
