#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include <fvm_writer.h>
#include <fvm_nodal_append.h>
#include <fvm_nodal_order.h>
#include <fvm_nodal_triangulate.h>
#include <mpi.h>

//#include "quickSort.h"
#include "creeMaillagePolygone2D.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
  Attention :
  -----------
  Les fichiers inclus par <metis.h> posent problème à minima sous Linux du à
  des redéclarations de fonctions de <stdlib.h>. On recopie dans ce cas le
  minimum nécessaire issu de <metis.h> (correspondant à METIS 4.0)
*/

typedef int idxtype;
void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *,
                              int *, int *, int *, int *, int *, idxtype *);
void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *,
                         int *, int *, int *, int *, int *, idxtype *);

#ifdef __cplusplus
}
#endif

#define ABS(a)     ((a) <  0  ? -(a) : (a)) 
#define MIN(a,b)   ((a) > (b) ?  (b) : (a)) 

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


static double random01() 
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / ABS(rsigna - rsignb);
  double resultat =   sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}


static int _partitionner(int tableau[], int p, int r) 
{
  int pivot = tableau[p], i = p-1, j = r+1;
  int temp;
  while(1) {
    do
      j--;
    while(tableau[j] > pivot);
    do
      i++;
    while(tableau[i] < pivot);
    if(i<j) {
      temp = tableau[i];
      tableau[i] = tableau[j];
      tableau[j] = temp;
    }
    else
      return j;
  }
  return j;
}


/* Fonctions publiques */

static void _quickSort(int tableau[], int p, int r) 
{
  int q;
  if(p<r) {
    q = _partitionner(tableau, p, r);
    _quickSort(tableau, p, q);
    _quickSort(tableau, q+1, r);
  }
  return;
}

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
                            int **eltsConnec) 
{
  MPI_Status status;

  int nRank;
  MPI_Comm_size(localComm, &nRank);

  int localRank;
  MPI_Comm_rank(localComm, &localRank);

  int nEdges = 0;
  int *globalVertexNum = NULL;
  int *globalEltNum = NULL;
  int *downConnectivity = NULL;
  int *edges = NULL;
  int *edgeToFace = NULL;
  int *neighbourPointer = NULL;
  int *neighbour = NULL;

  if (localRank == 0) {
    const double coefRand = .3;
    //const double coefRand = 0.45;
    srand(initRandom);
    
    /* Construction des coordonnes */
    /* --------------------------- */
    
    /* nx et ny doivent etre pair et superieur a 4 */

    int nx1 = nx;
    int ny1 = ny;
    
    if (nx1 < 4)
      nx1 = 4;
    
    if (ny1 < 4)
      ny1 = 4;
    
    if (nx1 % 2 != 0)
      nx1 += 1;
    
    if (ny1 % 2 != 0)
      ny1 += 1;
    
    /* cotes d un octogone */
    
    double cote1;
    if (nx1 == 4) 
      cote1 = (xmax - xmin) / (sqrt(2) + 1);
    else 
      cote1 = (xmax - xmin) / ((nx1 - 2)/2 + ((nx1 - 2)/2 - 1) * sqrt(2) + sqrt(2));
    double cote2;
    if (ny1 == 4) 
      cote2 = (ymax - ymin) / (sqrt(2) + 1);
    else 
      cote2 = (ymax - ymin) / ((ny1 - 2)/2 + ((ny1 - 2)/2 - 1) * sqrt(2) + sqrt(2));
    
    /* Sommets */
    
    *nVertex = 0;
    
    /* Allocation temporaire (surdimensionnee) */
    
    BFT_MALLOC(*meshCoords, 3*(nx1*ny1), double );
    
    int cpt = 0;
    int cptMax;
    if (ny1 > 4)
      cptMax = ny1 + (ny1-4)/2;
    else
      cptMax = ny1;
    
    double ycourant = ymin;
    int iy = 0;
    const double eps = 1e-5;
    while(cpt < cptMax) {
      if (cpt % 3 == 1 || cpt % 3 == 2) {
        for (int ix = 0; ix < nx1/2; ix++) {
          (*meshCoords)[3*(*nVertex)]   = xmin + ix * (1+sqrt(2)) * cote1;
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }
      }
      else if ((cpt % 3) == 0) {
        if ((cpt == (cptMax-1)) || (cpt == 0)) {
          (*meshCoords)[3*(*nVertex)]   = xmin;
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }
        
        (*meshCoords)[3*(*nVertex)]   = xmin + cote1/sqrt(2);
        (*meshCoords)[3*(*nVertex)+1] = ycourant;
        (*meshCoords)[3*(*nVertex)+2] = 0.;
        (*nVertex)++;
        
        for (int ix = 2; ix < nx1-1; ix++) {
          if (ix % 2 == 0) 
            (*meshCoords)[3*(*nVertex)] = (*meshCoords)[3*((*nVertex)-1)]+cote1;
          else
            (*meshCoords)[3*(*nVertex)] = (*meshCoords)[3*((*nVertex)-1)]+cote1*sqrt(2);
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }
        if ((cpt == (cptMax-1)) || (cpt == 0)) {
          (*meshCoords)[3*(*nVertex)]   = xmax;
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }
      }
      cpt++;
      if ((cpt % 3 == 1) || (cpt % 3 == 0))
        ycourant += cote2/sqrt(2);
      else 
        ycourant += cote2;
    }
    
    BFT_REALLOC(*meshCoords, 3*(*nVertex), double );
    for (int ix = 0; ix <(*nVertex) ; ix++) {
      if (ABS(xmin-(*meshCoords)[3*ix]) > eps &&
        ABS(xmax-(*meshCoords)[3*ix]) > eps &&
          ABS(ymax-(*meshCoords)[3*ix+1]) > eps &&
          ABS(ymin-(*meshCoords)[3*ix+1]) > eps) {
        (*meshCoords)[3*ix] += random01() * coefRand * cote1;
        (*meshCoords)[3*ix+1] += random01() * coefRand * cote2;
        
      }
      //(*meshCoords)[3*ix+2]=21*sin((*meshCoords)[3*ix+1]/10.)*sin((*meshCoords)[3*ix]/10.);
      (*meshCoords)[3*ix+2]=0.;
    }
    
    
    /* Construction des blocs d'elements */
    /* --------------------------------- */
    
    *nElts = 0;
    
    /* Allocation avec surestimation de la taille */
    /* Optimisation : Peut-etre reajuste */
    
    BFT_MALLOC(*eltsConnecPointer, nx1*ny1+1, int);
    BFT_MALLOC(*eltsConnec, 8*nx1*ny1, int);
    (*eltsConnecPointer)[0] = 0; 

    for(int itype = 0; itype < 3; itype++) {
      
      int itype1 = itype;
      
      if (order != 1) 
        itype1 = 2 - itype; 
 
      int n1;
      int n2;
      int n3;
      
      if (itype1 == 0) {
        /* Triangles */
        
        /* -- Premiere ligne */
        n1 = 1;
        n2 = n1+nx1;
        
        for (int ix = 0; ix < nx1/2; ix++) {
          int ix1 = 2 * ix;
          int ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n1+ix1;
          (*eltsConnec)[ideb+1] = n1+ix1+1;
          (*eltsConnec)[ideb+2] = n2+ix;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;  
        }
        
        /* -- Autres triangles (un a gauche un a droite */
        int nbLi = (ny1-4)/2;
        n1   = 1 + nx1 + nx1/2;
        n2   = n1 + nx1/2;
        n3   = n2 + nx1-2;
        
        for (int itri = 0; itri < nbLi; itri++) {
          int n4 = n1 + nx1/2 - 1;
          int n5 = n2 + nx1-2 - 1;
          int n6 = n3 + nx1/2 - 1;
          int ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n1;
          (*eltsConnec)[ideb+1] = n2;
          (*eltsConnec)[ideb+2] = n3;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;  
          
          ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n4;
          (*eltsConnec)[ideb+1] = n6;
          (*eltsConnec)[ideb+2] = n5;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;  
          n1 = n3 + nx1/2;
          n2 = n1 + nx1/2;
          n3 = n2 + nx1-2;
        }
        
        /* -- Derniere ligne */ 
        n2 = n1+nx1/2;
        
        for (int ix = 0; ix < nx1/2; ix++) {
          int ix1 = 2 * ix;
          int ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n1+ix;
          (*eltsConnec)[ideb+1] = n2+ix1+1;
          (*eltsConnec)[ideb+2] = n2+ix1;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;  
        }
        int nbTriangle = (*nElts);
      }
      
      else if (itype1 == 1) {
        
        /* Quadrangles */
        int nxQuad = (nx1-4)/2;
        int nyQuad = (ny1-4)/2;
        int nbQuadrangle = nxQuad*nyQuad;
        
        for (int iy = 0; iy < nyQuad; iy++) {
          for (int ix = 0; ix < nxQuad; ix++) {
            n1 = iy*(2*nx1-2) + nx1 + nx1/2 + 1 + ix + 1;
            n2 = iy*(2*nx1-2) + 2*nx1 + 1 + 2*ix + 1;
            n3 = n2 + 1;
            int n4 = iy*(2*nx1-2) + 3*nx1 - 2 + 1 + ix + 1 ;
            int ideb = (*eltsConnecPointer)[*nElts] ;
            (*eltsConnec)[ideb]   = n1;
            (*eltsConnec)[ideb+1] = n3;
            (*eltsConnec)[ideb+2] = n4;
            (*eltsConnec)[ideb+3] = n2;
            *nElts += 1;
            (*eltsConnecPointer)[*nElts] = ideb + 4;  
          }
        }
      }
      
      else if (itype1 == 2) {
        
        /* Polygones */
        int nxPoly = (nx1-2)/2;
        int nyPoly = (ny1-2)/2;
        int nbPoly = nxPoly * nyPoly; 
        
        int delta = 0;
        
        for (int iy = 0; iy < nyPoly; iy++) {
          
          if (iy == 1)
            delta += 2*nx1;
          else if (iy !=0)
            delta += 2*nx1-2;      
          for (int ix = 0; ix < nxPoly; ix++) {
            if (iy == 0)
              n1 = delta + 1 + 2*ix +1;
            else
              n1 = delta + 1 + 2*ix ;
            n2 = n1 + 1;
            n3 = iy*(2*nx1-2) + 1 + nx1 + ix;
            int n4 = n3 + 1;
            int n5 = iy*(2*nx1-2) + 1 + nx1 + nx1/2 + ix;
            int n6 = n5 + 1;
            int n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix;
            if (iy == (nyPoly - 1))
              n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix + 1;
            int n8 = n7 + 1;
            int ideb = (*eltsConnecPointer)[*nElts] ;
            (*eltsConnec)[ideb]   = n1;
            (*eltsConnec)[ideb+1] = n2;
            (*eltsConnec)[ideb+2] = n4;
            (*eltsConnec)[ideb+3] = n6;
            (*eltsConnec)[ideb+4] = n8;
            (*eltsConnec)[ideb+5] = n7;
            (*eltsConnec)[ideb+6] = n5;
            (*eltsConnec)[ideb+7] = n3;
            *nElts += 1;
            (*eltsConnecPointer)[*nElts] = ideb + 8;  
          }
        }
      }
    }
  }  
  
  /* Decoupage du maillage par Metis */
  /* ------------------------------- */
  
  if (nRank > 1) {
    if (localRank == 0) {
      
      int localNVertex = 0;
      int localNElts = 0;
      double *localCoords = NULL;

      /* Connectivite descendante */
      /* ------------------------ */
      
      /* Construction de la connectivite descendante et de la table des voisins 
         A sortir dans une fonction propore*/
      
      /* - Construction des aretes + connectivite descendante */
      /* Numerotation de 1 a n */
      /* tri des aretes par le min des sommets */
      
      int *aretes = NULL;
      BFT_MALLOC(aretes, (*eltsConnecPointer)[*nElts]*2, int) ;
      
      downConnectivity = NULL;
      BFT_MALLOC(downConnectivity, (*eltsConnecPointer)[*nElts], int) ;
      
      int **triAre = NULL;
      BFT_MALLOC(triAre, *nVertex, int*) ;
      
      int *nAreVertex = NULL;
      int nDefaultAreVertex = 8;
      BFT_MALLOC(nAreVertex, *nVertex, int) ;
      
      for (int i = 0; i < *nVertex; i++) {
        nAreVertex[i] = nDefaultAreVertex;
        BFT_MALLOC(triAre[i], nAreVertex[i], int) ;
        for (int j = 0; j < nAreVertex[i]; j++) {
          triAre[i][j] = -1;
        }
      }
      
      /* - Construction des aretes - Boucle sur les elements  */
      
      int iare = 0;
      for (int i = 0; i < *nElts; i++) {
        for (int j = (*eltsConnecPointer)[i]; j < (*eltsConnecPointer)[i+1]; j++) {
          aretes[2*iare]  = (*eltsConnec)[j];
          if (j != ((*eltsConnecPointer)[i+1] - 1))
            aretes[2*iare+1] = (*eltsConnec)[j+1];
          else 
            aretes[2*iare+1] = (*eltsConnec)[(*eltsConnecPointer)[i]];
          downConnectivity[j] = iare+1;
          
          int minVertex = MIN(aretes[2*iare], aretes[2*iare+1]) - 1;
          int k = 0;
          while ((k++ < nAreVertex[minVertex]) && (triAre[minVertex][k-1] > -1));      
          k--;
          if (k ==  nAreVertex[minVertex]){
            nAreVertex[minVertex] *= 2;
            BFT_REALLOC(triAre[minVertex], nAreVertex[minVertex], int) ;
            triAre[minVertex][k] = iare+1;
          }
          else if (triAre[minVertex][k] == -1) 
            triAre[minVertex][k] = iare+1;
          
          iare++;
        }
      }
      
      /* - Elimination des doublons - Boucle sur les sommets  */
      
      int *ptAretesCompactees = NULL;
      BFT_MALLOC(ptAretesCompactees, (*eltsConnecPointer)[*nElts], int) ;
      for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++)
        ptAretesCompactees[i] = -1;
      
      int iareCompactee = 1;
      for (int i = 0; i < *nVertex; i++) {
        int j = 0;
        while(triAre[i][j] > -1) {
          int k = j+1;
          int iare1 = triAre[i][j] - 1;
          if (ptAretesCompactees[iare1] == -1) {
            ptAretesCompactees[iare1] = iareCompactee++;
            if (k < nAreVertex[i]) {
              while(triAre[i][k] > - 1) {
                int iare2 = triAre[i][k] - 1;
                if (ptAretesCompactees[iare2] == -1)
                  if (((aretes[2*iare1] == aretes[2*iare2]) && 
                       (aretes[2*iare1+1] == aretes[2*iare2+1])) || 
                      ((aretes[2*iare1] == aretes[2*iare2+1]) && 
                       (aretes[2*iare1+1] == aretes[2*iare2])))
                    ptAretesCompactees[iare2] = ptAretesCompactees[iare1]; 
                k += 1;
              }
            }
          }
          j +=1;
        }
        if (triAre[i] != NULL)
          BFT_FREE(triAre[i]);
      }
      
      nEdges = iareCompactee-1;
      if (triAre != NULL)
        BFT_FREE(triAre);
      
      if (nAreVertex != NULL)
        BFT_FREE(nAreVertex);


      /* - Renumerotation connectivite descendante - Boucle sur la connectivite des éléments */
      
      for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++) 
        downConnectivity[i] = ptAretesCompactees[downConnectivity[i]-1];
      
      
      /* - Compression du tableau de description des aretes - Boucle sur la connectivitédes éléments */
      
      //edges = NULL;
      //BFT_MALLOC(edges, nEdges*2, int) ;
      
      //for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++) {
      //  edges[2*(ptAretesCompactees[i]-1)]   = aretes[2*i];
      //  edges[2*(ptAretesCompactees[i]-1)+1] = aretes[2*i+1];
      // }
      
      if (aretes != NULL)
        BFT_FREE(aretes);
      
      if (ptAretesCompactees != NULL)
        BFT_FREE(ptAretesCompactees);
      
      /* - Rangement des elements par couple suivant l'arete commune */
      
      edgeToFace = NULL;
      BFT_MALLOC(edgeToFace, nEdges*2, int);
      
      for (int i = 0; i < nEdges*2; i++)
        edgeToFace[i] = -1;
      
      for (int i = 0; i < *nElts; i++) {
        for (int j = (*eltsConnecPointer)[i]; j < (*eltsConnecPointer)[i+1]; j++) {
          if (edgeToFace[2*(downConnectivity[j]-1)] == -1 )
            edgeToFace[2*(downConnectivity[j]-1)] = i;
          else if (edgeToFace[2*(downConnectivity[j]-1)+1] == -1 )
            edgeToFace[2*(downConnectivity[j]-1)+1] = i;
          else
            bft_error(__FILE__, __LINE__, 0, "Arete a plus de 2 facettes !!!!\n");
        }
      }
      
      if (downConnectivity != NULL)
        BFT_FREE(downConnectivity);

      /* - Creation de la table des voisins (Numerotation de 1 a n)
         Le voisin d'une face de bord est temporairement marque a -1 */
      
      int *tmpVoisins = NULL;
      BFT_MALLOC(tmpVoisins, (*eltsConnecPointer)[*nElts], int) ;
      
      for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++)
        tmpVoisins[i] = -2;
      
      for (int i = 0; i < nEdges; i++) {
        int elt1 = edgeToFace[2*i];
        int elt2 = edgeToFace[2*i+1];
        if (elt1 != -1) {
          int j = 0;
          while(tmpVoisins[(*eltsConnecPointer)[elt1] + j++] != -2);
          tmpVoisins[(*eltsConnecPointer)[elt1] + --j] = elt2;
          }
        if (elt2 != -1) {
          int j = 0;
          while(tmpVoisins[(*eltsConnecPointer)[elt2] + j++] != -2);
          tmpVoisins[(*eltsConnecPointer)[elt2] + --j] = elt1;
        }
      }
      
      BFT_FREE(edgeToFace);

      /* - Filtrage des faces de bords */
      
      neighbour = NULL;
      BFT_MALLOC(neighbour, (*eltsConnecPointer)[*nElts], int) ;
      
      neighbourPointer = NULL;
      BFT_MALLOC(neighbourPointer, *nElts+1, int) ;
      
      int nVoisin = 0;
      neighbourPointer[0] = nVoisin;
      for (int i = 0; i < *nElts; i++) {
        for (int j = (*eltsConnecPointer)[i]; j < (*eltsConnecPointer)[i+1]; j++) {
          if (tmpVoisins[j] != -1)
            neighbour[nVoisin++] = tmpVoisins[j];
        }
        neighbourPointer[i+1] = nVoisin;
      }
      BFT_REALLOC(neighbour, nVoisin, int) ;
      
      if (tmpVoisins != NULL)
        BFT_FREE(tmpVoisins);
      
      /* Decoupage du maillage par Metis */
      /* ------------------------------- */
        
      int     wgtflag    = 0 ; /* Pas de pondération pour les faces ou cellules */
      int     numflag    = 0 ; /* Numérotation de 1 à n (type Fortran) */
      int     options[5] = {0, 3, 1, 1, 0} ; /* Par défaut si options[0] = 0 */
      int     edgecut    = 0 ; /* <-- nombre de faces sur la partition */
      
      int *numDomElt = NULL;
      BFT_MALLOC(numDomElt, (*nElts), int);
      assert(sizeof(idxtype) == sizeof(int));
      if (nRank < 8)
        
        METIS_PartGraphRecursive(nElts,
                                 (idxtype*)neighbourPointer,
                                 (idxtype*)neighbour,
                                 NULL,       /* vwgt   : poids des cellules */
                                 NULL,       /* adjwgt : poids des faces    */
                                 &wgtflag,
                                 &numflag,
                                 &nRank,
                                 options,
                                 &edgecut,
                                 numDomElt) ;
      
      else
        
        METIS_PartGraphKway(nElts,
                            (idxtype*)neighbourPointer,
                            (idxtype*)neighbour,
                            NULL,       /* vwgt   : poids des cellules */
                            NULL,       /* adjwgt : poids des faces    */
                            &wgtflag,
                            &numflag,
                            &nRank,
                            options,
                            &edgecut,
                            numDomElt) ;

      if (neighbour != NULL)
        BFT_FREE(neighbour);
      
      if (neighbourPointer != NULL)
        BFT_FREE(neighbourPointer);

      int  localEltsSize              = (*nElts)/nRank;   /* Estimation du nombre d'elements locaux 
                                                             pour allocation memoire */
      int  localEltsConnecSize        = 6*(*nElts)/nRank; /* Estimation de la taille 
                                                             de la connectivite locale 
                                                             pour allocation memoire*/
      
      
      int *localEltsConnecPointer = NULL;
        
      BFT_MALLOC(localEltsConnecPointer, localEltsSize+1, int);
      BFT_MALLOC(globalEltNum, localEltsSize, int);
      
      int *localEltsConnec = NULL;
      
      BFT_MALLOC(localEltsConnec, localEltsConnecSize, int);
      BFT_MALLOC(globalVertexNum, localEltsConnecSize, int);
      
      /* Pour chaque proc construction du maillage local envoi des donnees 
         On finit par le proc 0 */
      
      for (int irank = nRank-1; irank >= 0; irank--) {
        int idom = irank;
        localNElts = 0;
        for (int ielt = 0; ielt < (*nElts); ielt++) {
          if (numDomElt[ielt] == idom) {
            if (localEltsSize <= localNElts) {
              localEltsSize *= 2;
              BFT_REALLOC(localEltsConnecPointer, localEltsSize+1, int );
              BFT_REALLOC(globalEltNum, localEltsSize, int );
            }
            globalEltNum[localNElts++] = ielt+1;
          }
        }
        
        int tmpSize = 0;
        localEltsConnecPointer[0] = 0;
        for (int ielt = 0; ielt < localNElts; ielt++) {
          for (int i = (*eltsConnecPointer)[globalEltNum[ielt]-1]; i < (*eltsConnecPointer)[globalEltNum[ielt]]; i++) {
            if (localEltsConnecSize <= tmpSize) {
              localEltsConnecSize *= 2;
              BFT_REALLOC(localEltsConnec, localEltsConnecSize, int );
              BFT_REALLOC(globalVertexNum, localEltsConnecSize, int );
            }
            globalVertexNum[tmpSize]   = (*eltsConnec)[i];
            localEltsConnec[tmpSize++] = (*eltsConnec)[i];
          }
          localEltsConnecPointer[ielt+1] = tmpSize;
        }
        
        /* Tri et elimination des doublons */
        
        localNVertex = 0;
        _quickSort(globalVertexNum, 0, tmpSize-1);
        int i = 0;
        int val = globalVertexNum[i];
        localNVertex  = 1;

        bft_printf("Avant tri : \n");
        for (int i = 0; i < tmpSize; i++)
          bft_printf(" %i",globalVertexNum[i]);
        bft_printf("\n");

        do {
          i++;
          while (( i < tmpSize && (val == globalVertexNum[i])))
            i++;
          if (i < tmpSize) {
            val = globalVertexNum[i];
            globalVertexNum[localNVertex++] = val;
          } 
        } while (i < tmpSize);
        

        bft_printf("Apres tri : \n");
        for (int i = 0; i <localNVertex ; i++)
          bft_printf(" %i",globalVertexNum[i]);
        bft_printf("\n");

        /* Renumerotation de la connectivite */
        
        int *tmpRenum = NULL;
        BFT_MALLOC(tmpRenum, val, int); 
        
        for (int i = 0; i < val; i++)
          tmpRenum[i] = -1;
        
        for (int i = 0; i < localNVertex; i++)
          tmpRenum[globalVertexNum[i]-1] = i;
        
        for (int i = 0; i < localEltsConnecPointer[localNElts]; i++)
          localEltsConnec[i] = tmpRenum[localEltsConnec[i]-1]+1; 
        
        if (tmpRenum != NULL)
          BFT_FREE(tmpRenum);
        
        /* Coordonnees des sommets locaux */
        
        if (irank == nRank-1)
          BFT_MALLOC(localCoords, 3*localNVertex, double);
        else
          BFT_REALLOC(localCoords, 3*localNVertex, double);
        
        for (int i = 0; i < localNVertex; i++) {
          localCoords[3*i]   = (*meshCoords)[3*(globalVertexNum[i]-1)];
          localCoords[3*i+1] = (*meshCoords)[3*(globalVertexNum[i]-1)+1];
          localCoords[3*i+2] = (*meshCoords)[3*(globalVertexNum[i]-1)+2];
        }
        
        /* Envoi des donnees maillages si different du proc 0 */
        
        if (irank != 0) {
          
          bft_printf_flush();
          /* Envoi les infos concernant les sommets */
          int ierror = MPI_SUCCESS;
          ierror = MPI_Send(&localNVertex, 1, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS)
            bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
          
          ierror = MPI_Send(localCoords, 3*localNVertex, MPI_DOUBLE, irank, 0, localComm);
          if (ierror != MPI_SUCCESS)
            bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
          
          ierror = MPI_Send(globalVertexNum, localNVertex, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS)
            bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
          
          /* Envoi les infos concernant les elements */
          
          ierror = MPI_Send(&localNElts, 1, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS)
            bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
          
          ierror = MPI_Send(globalEltNum, localNElts, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS)
            bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
          
          ierror = MPI_Send(localEltsConnecPointer, localNElts+1, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS)
            bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
          
          ierror = MPI_Send(localEltsConnec, localEltsConnecPointer[localNElts], MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS)
            bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
        }
      }
      
      /* Suppression des donnees globales du maillage */
      
      if (numDomElt != NULL)
        BFT_FREE(numDomElt);

      if ((*meshCoords) != NULL)
        BFT_FREE(*meshCoords);
      
      if ((*eltsConnecPointer) != NULL)
        BFT_FREE(*eltsConnecPointer);
      
      if ((*eltsConnec) != NULL)
        BFT_FREE(*eltsConnec);
      
      /* Contenu du maillage local sur proc 0 */
      
      BFT_REALLOC(localCoords, 3*localNVertex, double);
      BFT_REALLOC(localEltsConnec, localEltsConnecPointer[localNElts], int);
      BFT_REALLOC(localEltsConnecPointer, localNElts+1, int);
      BFT_REALLOC(globalVertexNum, localNVertex, int);
      BFT_REALLOC(globalEltNum, localNElts, int);
      
      *nVertex = localNVertex;
      *nElts = localNElts;
      *meshCoords = localCoords;
      *eltsConnec = localEltsConnec;
      *eltsConnecPointer = localEltsConnecPointer;
      
    } 
  
    /* Reception des donnees maillage si different du proc 0 */
    
    if (localRank != 0) {
    
      /* Reception des infos concernant les sommets */
    
      int ierror = MPI_SUCCESS;
      ierror = MPI_Recv(nVertex, 1, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS)
        bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
      
      BFT_MALLOC(*meshCoords, 3*(*nVertex), double);
      
      ierror = MPI_Recv(*meshCoords, 3*(*nVertex), MPI_DOUBLE, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS)
        bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
      
      BFT_MALLOC(globalVertexNum, *nVertex, int);
      
      ierror = MPI_Recv(globalVertexNum, *nVertex, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS)
        bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
      
      /* Reception des infos concernant les elements */
      
      ierror = MPI_Recv(nElts, 1, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS)
        bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
    
      BFT_MALLOC(globalEltNum, *nElts, int);
      
      ierror = MPI_Recv(globalEltNum, *nElts, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS)
        bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
      
      BFT_MALLOC(*eltsConnecPointer, (*nElts)+1, int);
      
      ierror = MPI_Recv(*eltsConnecPointer, (*nElts)+1, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS)
        bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
    
      BFT_MALLOC(*eltsConnec, (*eltsConnecPointer)[(*nElts)], int);
      
      ierror = MPI_Recv(*eltsConnec, (*eltsConnecPointer)[(*nElts)], MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS)
        bft_error(__FILE__, __LINE__, 0, "Erreur MPI\n");
    }
  }

  bft_printf("connectivite : \n");
  for (int i = 0; i < *nElts ; i++) {
    for (int j = (*eltsConnecPointer)[i]; j < (*eltsConnecPointer)[i+1] ; j++)
      bft_printf(" %i", (*eltsConnec)[j]);
    bft_printf("\n");
  }



  if (globalEltNum != NULL)
    BFT_FREE(globalEltNum);

  if (globalVertexNum != NULL)
    BFT_FREE(globalVertexNum);

  /* Construction de la structure nodale */
  /* ----------------------------------- */
  
/*   int nbTriangle   = 0; */
/*   int nbQuadrangle = 0; */
/*   int nbPoly       = 0; */

/*   for (int i = 0; i < *nElts; i++) { */
/*     int nCurrentEltVertex = (*eltsConnecPointer)[i+1] - (*eltsConnecPointer)[i]; */
/*     if (nCurrentEltVertex == 3) { */
/*       if ((nbQuadrangle != 0) || (nbPoly != 0))  */
/*         bft_error(__FILE__, __LINE__, 0, "Les elements ne sont pas ordonnes\n"); */
/*       ++nbTriangle; */
/*     } */
    
/*     else if (nCurrentEltVertex == 4) { */
/*       if (nbPoly != 0)  */
/*         bft_printf("Les elements ne sont pas ordonnes\n"); */
/*       ++nbQuadrangle; */
/*     } */
    
/*     else if (nCurrentEltVertex > 4) { */
/*       ++nbPoly; */
/*     } */
    
/*     else */
/*       bft_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n"); */
    
/*   } */
}


void PROCF(creemaillagepolygone2d_f, CREEMAILLAGEPOLYGONE2D_F)(int *order,
							       MPI_Fint* localFComm,
                                                               /*MPI_Comm *localComm,*/
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
                                                               ) 
{
  MPI_Comm localComm = MPI_Comm_f2c(*localFComm);
  int nVertex_f = *nVertex;
  int nElts_f = *nElts;

  double *meshCoords = NULL;
  int    *eltsConnecPointer = NULL;
  int    *eltsConnec = NULL;

  creeMaillagePolygone2D(*order,
                         localComm,
                         *xmin, 
                         *xmax, 
                         *ymin, 
                         *ymax,
                         *initRandom,
                         *nx,
                         *ny,
                         nVertex,
                         &meshCoords,
                         nElts,
                         &eltsConnecPointer,
                         &eltsConnec);

  if (nVertex_f < *nVertex)
    bft_error(__FILE__, __LINE__, 0, "Augmenter le nombre de sommets Fortran a : %i \n", *nVertex);

  if (nElts_f < *nElts)
    bft_error(__FILE__, __LINE__, 0, "Augmenter le nombre d'elements a : %i \n", *nElts);

  if (*lEltsConnecPointer_f < eltsConnecPointer[*nElts])
    bft_error(__FILE__, __LINE__, 0, "Augmenter la taille du tableau de connectivite a : %i \n", eltsConnecPointer[*nElts]);

  for(int i = 0; i < 3*(*nVertex); i++)
    meshCoords_f[i] = meshCoords[i];

  for(int i = 0; i < *nElts + 1; i++)
    eltsConnecPointer_f[i] = eltsConnecPointer[i];

  for(int i = 0; i < eltsConnecPointer[*nElts]; i++)
    eltsConnec_f[i] = eltsConnec[i];

  BFT_FREE(meshCoords);
  BFT_FREE(eltsConnecPointer);
  BFT_FREE(eltsConnec);
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
