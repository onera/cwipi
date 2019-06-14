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
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>

#include "pdm_poly_surf_gen.h"
#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "grid_mesh_pdm.h"

#define ABS(a)     ((a) <  0  ? -(a) : (a))


/*----------------------------------------------------------------------
 *                                                                     
 * Return a random number in [-1, 1]                                    
 *                                                                     
 * parameters:
 *   coupling_id         <-- Coupling id               
 *   nNotLocatedPoints   <-- Number of not located points
 *---------------------------------------------------------------------*/

static double _random01()
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / ABS(rsigna - rsignb);
  double resultat =   sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}








int
_create_split_mesh
(
 int           imesh,
 PDM_MPI_Comm  pdm_mpi_comm,
 PDM_g_num_t  nVtxSeg,
 double        length,
 int           nPart,
PDM_part_split_t           method,
 int           haveRandom,
 PDM_g_num_t   *nGFace,
 PDM_g_num_t   *nGVtx,
 PDM_g_num_t   *nGEdge,
 int           *nTPart,
 int           *nEdgeGroup,
 int           **nVtx,
 int           **nElts,
 int           ***connecIdx,
 int           ***connec,
 double        ***coords
)
{
  struct timeval t_elaps_debut;

  int myRank;
  int numProcs;

  PDM_MPI_Comm_rank (pdm_mpi_comm, &myRank);
  PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

  double        xmin = 0.;
  double        xmax = length;
  double        ymin = 0.;
  double        ymax = length;
  PDM_g_num_t  nx = nVtxSeg;
  PDM_g_num_t  ny = nVtxSeg;
  int           dNFace;
  int           dNVtx;
  int           dNEdge;
  int          *dFaceVtxIdx;
  PDM_g_num_t *dFaceVtx;
  double       *dVtxCoord;
  PDM_g_num_t *dfaceEdge;
  PDM_g_num_t *dedgeVtx;
  PDM_g_num_t *dEdgeFace;
  int          *dEdgeGroupIdx;
  PDM_g_num_t   *dEdgeGroup;

  int           initRandom = 0;

  /*
   *  Create mesh i
   */

  if (imesh == 1) {
    nx *= 2;
    ny *= 2;
  }

  ++initRandom;

  gettimeofday(&t_elaps_debut, NULL);

  PDM_poly_surf_gen (pdm_mpi_comm,
                     xmin,
                     xmax,
                     ymin,
                     ymax,
                     haveRandom,
                     initRandom,
                     nx,       
                     ny,
                     nGFace,
                     nGVtx,  
                     nGEdge,
                     &dNVtx,
                     &dVtxCoord,
                     &dNFace,
                     &dFaceVtxIdx,
                     &dFaceVtx,
                     &dfaceEdge,
                     &dNEdge,
                     &dedgeVtx,
                     &dEdgeFace,
                     nEdgeGroup,
                     &dEdgeGroupIdx,
                     &dEdgeGroup);
  

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc (dNFace*sizeof(int));
  int *dedgeVtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

  dedgeVtxIdx[0] = 0;
  for (int i = 0; i < dNEdge; i++) {
    dedgeVtxIdx[i+1] = 2 + dedgeVtxIdx[i];
  }

  /*
   *  Split mesh i
   */

  int ppartId;

  int nPropertyCell = 0;
  int *renum_properties_cell = NULL;
  int nPropertyFace = 0;
  int *renum_properties_face = NULL;
  
  
  int dNCell=dNFace;
  
  int* dCellFaceIdx = (int*)malloc(sizeof(int)*(1+dNCell));
  dCellFaceIdx[0]=0;
  for(int i=0;i<dNCell;i++){
    dCellFaceIdx[i+1]=dCellFaceIdx[i]+1;
  }
  
  int* dCellFace = (int*)malloc( sizeof(int) * dCellFaceIdx[dNCell] );
  for(int i=0;i<dNCell;i++){
    for(int j=dCellFaceIdx[i]; j < dCellFaceIdx[i+1]; j++){
      dCellFace[j]=j;
    }
  }

  PDM_part_create (&ppartId,
                   pdm_mpi_comm,
                   method,
                   "PDM_PART_RENUM_CELL_NONE",
                   "PDM_PART_RENUM_FACE_NONE",
                   nPropertyCell,
                   renum_properties_cell,
                   nPropertyFace,
                   renum_properties_face,
                   nPart,
                   dNCell,
                   dNFace,
                   dNVtx,
                   0,
                   dCellFaceIdx,
                   dCellFace   ,
                   NULL,
                   NULL,
                   0,
                   NULL,
                   NULL,
                   dFaceVtxIdx,
                   dFaceVtx,
                   NULL,
                   dVtxCoord,
                   NULL,
                   NULL,
                   NULL);
while(1==1){}
 /* free (dCellPart);

  free (dVtxCoord);
  free (dFaceVtxIdx);
  free (dFaceVtx);
  free (dfaceEdge);
  free (dedgeVtxIdx);
  free (dedgeVtx);
  free (dEdgeFace);
  free (dEdgeGroupIdx);
  free (dEdgeGroup);
*/

  double **tmpcoords            = (double**)malloc(sizeof(double*)*nPart);
  int    **tmpconnecIdx         = (int**)malloc(sizeof(int*)*nPart);
  int    **tmpconnec            = (int**)malloc(sizeof(int*)*nPart);
  /* Bandwidth */
  int* tmpnVtx          = (int*)malloc(sizeof(int)*nPart);
  int* tmpnElts            = (int*)malloc(sizeof(int)*nPart);
  
  for (int ipart = 0; ipart < nPart; ipart++) {

    int nFace;
    int nEdge;
    int nEdgePartBound;
    int nProc;
    int sfaceEdge;
    int sedgeVtx;
    int sEdgeGroup;
    int nEdgeGroup2;
 
    int too1,too2;
    PDM_part_part_dim_get (ppartId,
                           ipart,
                           &(tmpnElts[ipart]),
                           &nEdge,
                           &nEdgePartBound,
                           &(tmpnVtx[ipart]),
                           &nProc,
                           nTPart,
                           &sfaceEdge,
                           &sedgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);
    
    int          *cellTag;
    int          *cellFaceIdx;
    int          *faceEdgeIdx;
    int          *faceEdge;
    
    int          *edgeVtxIdx;
    int          *edgeVtx;
    int          *cellFace;
    PDM_g_num_t *cellLNToGN;
    int          *faceTag;
    int          *faceCell;
    int          *faceVtxIdx;
    int          *faceVtx;
    PDM_g_num_t *faceLNToGN;
    int          *facePartBoundProcIdx;
    int          *facePartBoundPartIdx;
    int          *facePartBound;
    int          *vtxTag;
    double       *vtx;
    PDM_g_num_t *vtxLNToGN;
    int          *faceGroupIdx;
    int          *faceGroup;
    PDM_g_num_t *faceGroupLNToGN;

   // (*coords)           [i_part] ;//= (double*)malloc(sizeof(double)*nVtx );
   // (*eltsConnecPointer)[i_part] ;//= (int*)   malloc(sizeof(int)   *nFace);
   // (*eltsConnec)       [i_part] ;//= (int*)   malloc(sizeof(int)   *nFace);

    PDM_part_part_val_get(ppartId,
                          ipart,        //Partition ID
                          &cellTag,     //3D
                          &cellFaceIdx, //3D
                          &cellFace,    //3D
                          &cellLNToGN,  //3D
                          &faceTag,     //Uselesss
                          &faceCell,    // ??
                          &faceVtxIdx,  //Connectivity IDX
                          &faceVtx,     //Connectivity
                          &faceLNToGN,  //Element global Numbering 
                          &facePartBoundProcIdx,
                          &facePartBoundPartIdx,
                          &facePartBound,
                          &vtxTag,    //UseLess
                          &(tmpcoords[ipart]),       //Vertices coordinates
                          &vtxLNToGN, //Vertices global numbering
                          &faceGroupIdx,
                          &faceGroup,
                          &faceGroupLNToGN); 

  tmpconnecIdx[ipart] = faceVtxIdx ;
  tmpconnec   [ipart] = faceVtx    ;

/*
   tmpconnecIdx[ipart] = (int*) malloc(sizeof(int) *(1+tmpnElts[ipart]));

   int idx =0;
   tmpconnecIdx[ipart][0]=0;
   for(int i=0;i<tmpnElts[ipart];i++){
     for(int j=faceEdgeIdx[i]; j<faceEdgeIdx[i+1]; j++){
       int edge = faceEdge[j];
       for(int k=edgeVtxIdx[edge]; k<edgeVtxIdx[edge+1]; k++){
         idx++; 
       }
     }
     tmpconnecIdx[ipart][i+1] = idx;
   }
   
   tmpconnec   [ipart] = (int*) malloc(sizeof(int) * tmpconnecIdx[ipart][tmpnElts[ipart]]);

   idx =0;
   for(int i=0;i<tmpnElts[ipart];i++){
     for(int j=faceEdgeIdx[i]; j<faceEdgeIdx[i+1]; j++){
       int edge = faceEdge[j];
       for(int k=edgeVtxIdx[edge]; k<edgeVtxIdx[edge+1]; k++){
         tmpconnec[ipart][idx] = edgeVtx[k];
         idx++; 
       }
     }
   }
*/
     
  }
  *nVtx  = tmpnVtx;
  *nElts = tmpnElts;
  
  *coords = tmpcoords;
  *connecIdx = tmpconnecIdx; 
  *connec    = tmpconnec; 

  return ppartId;
}




/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   xmin                <-- grid xmin (global) 
 *   xmax                <-- grid xmax (global)
 *   ymin                <-- grid ymin (global)
 *   ymax                <-- grid ymax (global)
 *   randLevel           <-- random level
 *   nVertexSeg          <-- number of vertices in X and Y
 *   nPartitionSeg       <-- number of partitions in X and Y
 *   nVertex             --> number of vertices of the local partition
 *   coords              --> coordinates of vertices of the local partition
 *   nElts               --> number of elements of the local partition
 *   eltsConnecPointer   --> connectivity index of the local partition 
 *   eltsConnec          --> connectivity of the local partition
 *   localComm           <-- MPI comm of the global grid
 *---------------------------------------------------------------------*/


void  grid_mesh_pdm(double xmin, 
                double xmax, 
                double ymin, 
                double ymax, 
                double randLevel,
                int nVertexSeg,
                int nPart, 
                int** nVtx,
                double ***coords,
                int** nElts ,
                int ***eltsConnecPointer,
                int ***eltsConnec,
                MPI_Comm localComm)
{

  /* Compute local partition bounds with random level */

  int           dNEdge;
  int           haveRandom = 0;
  int           initRandom = 0;
  PDM_g_num_t   nGFace;
  PDM_g_num_t   nGVtx; 
  PDM_g_num_t   nGEdge;
  int           *nEdgeGroup;  

   ++initRandom;
     
   PDM_MPI_Comm pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(&localComm);


  /*
   *  Split mesh i
   */



  int nPropertyCell = 0;
  int *renum_properties_cell = NULL;
  int nPropertyFace = 0;
  int *renum_properties_face = NULL;

#ifdef PDM_HAVE_PARMETIS  
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH  
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif  
  int nFaceGroup;


  /*
   *  Create a partitioned mesh
   */
  
  
  
  int imesh = 0;
  int nTPart;
  double length = (xmax - xmin)/(double)nVertexSeg ;
  int ppartId =_create_split_mesh (imesh,
                      pdm_localComm,
                      nVertexSeg,
                      length,
                      nPart,
                      method,
                      haveRandom,
                      &nGFace,
                      &nGVtx,
                      &nGEdge,
                      &nTPart,
                      &nEdgeGroup,
                      nVtx,
                      nElts,
                      eltsConnecPointer,
                      eltsConnec,
                      coords
                      );
                      
                      
   printf("nPart %i nGFace %i\n",nPart,nGFace);   

}



/* rotation d'un angle alpha en degre */
void mesh_rotate_pdm(double *coords, int npt, double alpha) {
  double degrad = acos(-1.0)/180.;
  double sina, cosa, x, y;
  
  alpha = alpha * degrad;
  sina = sin(alpha);
  cosa = cos(alpha);

  for (int i = 0; i < npt; i++) {
    x = coords[i*3];
    y = coords[i*3 +1];
    coords[i*3]   = cosa*x - sina*y;
    coords[i*3+1] = sina*x + cosa*y;
  }
}

void carre2rond_pdm(double xmin, double xmax, double ymin, double ymax, double *coords, int nVertex)  {
  double xc, yc , x , y;
  double hyp, l;
  // tolerance pour eviter division par ZERO
  const double toler = 0.00000001;

  // remarque : ne fonctionne que dans un plan XY
  // coord centre
  xc = (xmax + xmin)/2.;
  yc = (ymax + ymin)/2.;

  for (int i = 0; i < nVertex ; i++) {
    x = coords[i*3] - xc;
    y = coords[i*3+1] - yc;

    hyp = sqrt(x*x + y*y);
    if (hyp > toler) {
      if (fabs(x) > fabs(y)) {
    l = fabs(cos(acos(x / hyp)));
      } else {
    l = fabs(sin(asin(y / hyp)));
      }
      coords[i*3] = xc + x*l;
      coords[i*3+1] = yc + y*l;
    }
  }
}

