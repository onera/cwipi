#include "coo_baryc.h"
#include <assert.h>
#include <stdio.h>
#define ABS(a)     ((a) <  0  ? -(a) : (a))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum {
  X,
  Y,
  Z
} ;


#define PRODUIT_VECTORIEL(prod_vect, vect1, vect2) ( \
prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z],   \
prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z],   \
prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y]   )


#define PRODUIT_SCALAIRE(vect1, vect2)                        ( \
  vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z] )


#define MODULE(vect) \
     sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])


#define DETERMINANT(vect1, vect2, vect3) ( \
   ((vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z]) * vect3[X]) \
 + ((vect2[X] * vect1[Z] - vect1[X] * vect2[Z]) * vect3[Y]) \
 + ((vect1[X] * vect2[Y] - vect2[X] * vect1[Y]) * vect3[Z]) )

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*----------------------------------------------------------------------------
 *  Fonction qui projette une face sur un plan parallèle à la face.
 *  Ce plan est ensuite assimilé au plan (Oxy).
 *----------------------------------------------------------------------------*/

static void projection_plan_moyen
(
 const int     nbr_som_fac,
 double *const coo_som_fac,
 double *const coo_point_dist
)
{

  int    icoo ;
  int    isom ;
  int    itri ;

  double   cost ;
  double   sint ;

  double   coo_tmp ;

  double   vect1[3] ;
  double   vect2[3] ;

  double  prod_vect[3] ;

  double  barycentre_fac[3] ;
  double  normale_fac[3] ;

  double *coo_som_fac_tmp = NULL;
  double  coo_point_dist_tmp[3] ;


  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Calcul des coordonnées du barycentre B du polygone P */
  /*======================================================*/

  for (icoo = 0 ; icoo < 3 ; icoo++) {

    barycentre_fac[icoo] = 0. ;

    for (isom = 0 ; isom < nbr_som_fac ; isom++)
      barycentre_fac[icoo] += coo_som_fac[3*isom+icoo] ;

    barycentre_fac[icoo] /= nbr_som_fac ;

  }

  for (icoo = 0 ; icoo < 3 ; icoo++)
    normale_fac[icoo] = 0. ;

  /* Calcul de la normale */
  /*======================*/

  for (itri = 0 ; itri < nbr_som_fac ; itri++) {

    for (icoo = 0 ; icoo < 3 ; icoo++) {

      vect1[icoo] = coo_som_fac[3*itri+icoo] - barycentre_fac[icoo] ;

      if (itri < nbr_som_fac - 1)
        vect2[icoo] = coo_som_fac[3*(itri+1)+icoo] - barycentre_fac[icoo] ;
      else
        vect2[icoo] = coo_som_fac[icoo] - barycentre_fac[icoo] ;

    }

    normale_fac[0] += vect1[1] * vect2[2] - vect2[1] * vect1[2] ;
    normale_fac[1] += vect2[0] * vect1[2] - vect1[0] * vect2[2] ;
    normale_fac[2] += vect1[0] * vect2[1] - vect2[0] * vect1[1] ;

  }


  /* Projection dans un plan parallèle à la face */
  /*=============================================*/

  /* On ramène l'origine au centre de gravité de la fac */

  for (isom = 0 ; isom < nbr_som_fac ; isom++)
    for (icoo = 0 ; icoo < 3 ; icoo++)
      coo_som_fac[3*isom+icoo] -= barycentre_fac[icoo] ;

  for (icoo = 0 ; icoo < 3 ; icoo++)
    coo_point_dist[icoo] -= barycentre_fac[icoo] ;

  if (ABS(normale_fac[0]) > 1.e-12 || ABS(normale_fac[1]) > 1.e-12) {

    /* Première rotation d'axe (Oz) et d'angle (Ox, proj normale sur Oxy) */

    BFT_MALLOC(coo_som_fac_tmp, 3*nbr_som_fac, double) ;

    vect1[0] = 1. ;
    vect1[1] = 0. ;
    vect1[2] = 0. ;

    vect2[0] = normale_fac[0] ;
    vect2[1] = normale_fac[1] ;
    vect2[2] = 0. ;

    PRODUIT_VECTORIEL(prod_vect, vect1, vect2) ;

    cost = PRODUIT_SCALAIRE(vect1, vect2) / MODULE(vect2) ;

    if (prod_vect[2] > 0.)
      sint =  MODULE(prod_vect) / MODULE(vect2) ;
    else
      sint = -MODULE(prod_vect) / MODULE(vect2) ;

    for (isom = 0 ; isom < nbr_som_fac ; isom++) {

      coo_som_fac_tmp[3*isom] =
         cost*coo_som_fac[3*isom] + sint*coo_som_fac[3*isom+1] ;
      coo_som_fac_tmp[3*isom+1] =
        -sint*coo_som_fac[3*isom] + cost*coo_som_fac[3*isom+1] ;
      coo_som_fac_tmp[3*isom+2] = coo_som_fac[3*isom+2] ;

    }

    coo_point_dist_tmp[0] = cost*coo_point_dist[0] + sint*coo_point_dist[1] ;
    coo_point_dist_tmp[1] = -sint*coo_point_dist[0] + cost*coo_point_dist[1] ;
    coo_point_dist_tmp[2] = coo_point_dist[2] ;

    /* Deuxième rotation d'axe (Oy) et d'angle (Oz', proj normale sur Ox'z) */

    vect1[0] =  0. ;
    vect1[1] =  0. ;
    vect1[2] =  1. ;

    vect2[0] =
      sqrt(normale_fac[0]*normale_fac[0] + normale_fac[1]*normale_fac[1]) ;
    vect2[1] = 0. ;
    vect2[2] = normale_fac[2] ;

    PRODUIT_VECTORIEL(prod_vect, vect1, vect2) ;

    cost = PRODUIT_SCALAIRE(vect1, vect2) / MODULE(vect2) ;

    if (prod_vect[2] > 0.)
      sint =  MODULE(prod_vect) / MODULE(vect2) ;
    else
      sint = -MODULE(prod_vect) / MODULE(vect2) ;


    for (isom = 0 ; isom < nbr_som_fac ; isom++) {

      coo_som_fac[3*isom] =
         cost*coo_som_fac_tmp[3*isom] + sint*coo_som_fac_tmp[3*isom + 2] ;
      coo_som_fac[3*isom+1] = coo_som_fac_tmp[3*isom+1] ;
      coo_som_fac[3*isom+2] = 0. ;

    }

    coo_point_dist[0] =
      cost*coo_point_dist_tmp[0] + sint*coo_point_dist_tmp[2] ;
    coo_point_dist[1] = coo_point_dist_tmp[1] ;
    coo_point_dist[2] = 0. ;


    BFT_FREE(coo_som_fac_tmp) ;

  }
  else {

    /* On écrase seulement la coordonnée z du sommet, en intervertissant
       éventuellement les coordonnées dans le plan de projection (Oxy).  */

    if (normale_fac[2] > 0.) {
      for (isom = 0 ; isom < nbr_som_fac ; isom++)
        coo_som_fac[3*isom+2] = 0. ;

      coo_point_dist[2] = 0. ;
    }
    else {
      for (isom = 0 ; isom < nbr_som_fac ; isom++) {
        coo_tmp = coo_som_fac[3*isom] ;
        coo_som_fac[3*isom] = coo_som_fac[3*isom+1] ;
        coo_som_fac[3*isom+1] = coo_tmp ;
        coo_som_fac[3*isom+2] = 0. ;
      }
      coo_tmp = coo_point_dist[0] ;
      coo_point_dist[0] = coo_point_dist[1] ;
      coo_point_dist[1] = coo_tmp ;
      coo_point_dist[2] = 0. ;
    }
  }
}


void coo_baryc(const fvm_locator_t* locator,
               const int   nMeshCoords,
               const double *meshCoords,
               const int   nElts,
               const int  *nMeshElts,
               const int  *meshElts,
               int  *n_dist_points,
               int    **nDistBarCoords,
               double **distBarCoords)
{
  /* Boucle sur les points distants */
  *n_dist_points = fvm_locator_get_n_dist_points(locator);
  const fvm_lnum_t *dist_locations = fvm_locator_get_dist_locations(locator);
  const fvm_coord_t *dist_coords = fvm_locator_get_dist_coords(locator);
  fvm_coord_t coo_point_dist[3];

  /* Tableaux locaux */

  const double eps = 1e-5;
  double *coo_som_fac = NULL;
  double *s           = NULL;
  double *dist        = NULL;
  double *aire        = NULL;
  double *proScal     = NULL;
  int tailleDistBarCoords;
  int taille_coo_som_fac;

  BFT_MALLOC(*nDistBarCoords, (*n_dist_points)+1, int ) ;
  //bft_printf("nDistBarCoords %d\n",(*n_dist_points)+1);

  (*nDistBarCoords)[0] = 0;
  tailleDistBarCoords = 4*(*n_dist_points);
  BFT_MALLOC(*distBarCoords, tailleDistBarCoords, double) ;
  //bft_printf("distBarCoords %d\n",tailleDistBarCoords);

  for (int ipoint =  0; ipoint < *n_dist_points; ipoint++ ) {
    //bft_printf("-- Etude du point : %d \n", ipoint);

    /* Initialisation - Copie locale */

    int isOnEdge = 0;
    int isVertex = 0;
    int ielt = dist_locations[ipoint] - 1;
    int nbr_som_fac = nMeshElts[ielt+1] - nMeshElts[ielt];
    coo_point_dist[0] = dist_coords[3*ipoint];
    coo_point_dist[1] = dist_coords[3*ipoint + 1];
    coo_point_dist[2] = dist_coords[3*ipoint + 2];

    if (ipoint == 0) {
      taille_coo_som_fac = 3*nbr_som_fac;
      BFT_MALLOC(coo_som_fac, taille_coo_som_fac, double) ;
      BFT_MALLOC(s, taille_coo_som_fac, double) ;
      BFT_MALLOC(dist, nbr_som_fac, double) ;
      BFT_MALLOC(aire, nbr_som_fac, double) ;
      BFT_MALLOC(proScal, nbr_som_fac, double) ;
    }
    else
      if (taille_coo_som_fac < 3*nbr_som_fac) {
        taille_coo_som_fac = 3*nbr_som_fac;
        BFT_REALLOC(coo_som_fac, taille_coo_som_fac, double);
        BFT_REALLOC(s, taille_coo_som_fac, double) ;
        BFT_REALLOC(dist, nbr_som_fac, double) ;
        BFT_REALLOC(aire, nbr_som_fac, double) ;
        BFT_REALLOC(proScal, nbr_som_fac, double) ;
      }

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      coo_som_fac[3*isom]   = meshCoords[3*(meshElts[nMeshElts[ielt]+isom]-1)];
      coo_som_fac[3*isom+1] = meshCoords[3*(meshElts[nMeshElts[ielt]+isom]-1)+1];
      coo_som_fac[3*isom+2] = meshCoords[3*(meshElts[nMeshElts[ielt]+isom]-1)+2];
    }

    /* Projection sur un plan moyen */

    projection_plan_moyen (nbr_som_fac, coo_som_fac, coo_point_dist);

    /* Calcul des coordonnnees barycentriques */

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      s[3*isom]   = coo_som_fac[3*isom]   - coo_point_dist[0];
      s[3*isom+1] = coo_som_fac[3*isom+1] - coo_point_dist[1];
      s[3*isom+2] = coo_som_fac[3*isom+2] - coo_point_dist[2];
      dist[isom] = sqrt(s[3*isom]*s[3*isom] +
                        s[3*isom+1]*s[3*isom+1] +
                        s[3*isom+2]*s[3*isom+2]);
    }

    int currentVertex;
    for (int isom = 0; isom < nbr_som_fac; isom++) {
      if (isom != (nbr_som_fac - 1)) {
        aire[isom] = s[3*isom]*s[3*(isom+1)+1] - s[3*(isom+1)]*s[3*isom+1];
        proScal[isom] = s[3*isom] * s[3*(isom+1)] +
                        s[3*isom+1] * s[3*(isom+1)+1] +
                        s[3*isom+2] * s[3*(isom+1)+2];
      }
      else {
        aire[isom] = s[3*isom]*s[1] - s[0]*s[3*isom+1];
        proScal[isom] = s[3*isom] * s[0] +
                        s[3*isom+1] * s[1] +
                        s[3*isom+2] * s[2];
      }
      if (dist[isom] <= eps) {
        isVertex = 1;
        currentVertex = isom;
        break;
      }
      /* faire un test avec eps pour proScal */
      else if (aire[isom] <= eps && proScal[isom] < 0.) {
        isOnEdge = 1;
        currentVertex = isom;
        break;
      }
    }

    /* Mise a jour de la taille du tableau de stockage des coordonnees barycentriques */

    (*nDistBarCoords)[ipoint+1] = (*nDistBarCoords)[ipoint] + nbr_som_fac;

    if (tailleDistBarCoords <= (*nDistBarCoords)[ipoint+1]) {
      tailleDistBarCoords *= 2;
      BFT_REALLOC(*distBarCoords, tailleDistBarCoords, double) ;
    }

    /* Le point distant est un sommet */

    if (isVertex) {

      for (int isom = 0; isom < nbr_som_fac; isom++)
        (*distBarCoords)[(*nDistBarCoords)[ipoint]+isom] = 0.;
      (*distBarCoords)[(*nDistBarCoords)[ipoint]+currentVertex] = 1.;
    }

    /* Le point distant est sur arete */

    else if (isOnEdge) {
      for (int isom = 0; isom < nbr_som_fac; isom++)
        (*distBarCoords)[(*nDistBarCoords)[ipoint]+isom] = 0.;

      int nextPoint;
      if (currentVertex == (nbr_som_fac - 1))
        nextPoint = 0;
      else
        nextPoint = currentVertex + 1;

      (*distBarCoords)[(*nDistBarCoords)[ipoint]+currentVertex] = dist[nextPoint]     / (dist[nextPoint]+dist[currentVertex]);
      (*distBarCoords)[(*nDistBarCoords)[ipoint]+nextPoint]     = dist[currentVertex] / (dist[nextPoint]+dist[currentVertex]);

    }

    /* Cas general */

    else {
      double sigma = 0;
      for (int isom = 0; isom < nbr_som_fac; isom++) {
        double coef = 0.;
        int previousVertex;
        int nextVertex;
        if (isom != 0)
          previousVertex = isom - 1;
        else
          previousVertex = nbr_som_fac - 1;
        if (isom < nbr_som_fac - 1)
          nextVertex = isom + 1;
        else
          nextVertex = 0;
        if (ABS(aire[previousVertex]) > eps)
          coef += (dist[previousVertex] - proScal[previousVertex]/dist[isom]) / aire[previousVertex];
        if (ABS(aire[previousVertex]) > eps)
          coef += (dist[nextVertex] - proScal[isom]/dist[isom]) / aire[isom];
        sigma += coef;
        (*distBarCoords)[(*nDistBarCoords)[ipoint]+isom] = coef;
      }
      assert(ABS(sigma) >= eps );
      for (int isom = 0; isom < nbr_som_fac; isom++)
        (*distBarCoords)[(*nDistBarCoords)[ipoint]+isom] /= sigma;
    }
  }

  if (coo_som_fac != NULL) {
    BFT_FREE(coo_som_fac);
    BFT_FREE(s) ;
    BFT_FREE(aire) ;
    BFT_FREE(dist) ;
    BFT_FREE(proScal) ;
  }
  BFT_REALLOC(*distBarCoords, (*nDistBarCoords)[*n_dist_points], double) ;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
