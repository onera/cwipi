/*============================================================================
 * Sorties Cedre
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_writer.h"
#include "pdm_writer_priv.h"
#include "pdm_writer_ensight.h"
#include "pdm_binary_search.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_fortran_to_c_string.h"

// validation gnum
#include "pdm_gnum.h"
//fin validation gnum


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des macros locales
 *============================================================================*/

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*----------------------------------------------------------------------------
 * NOMBRE DE BLOCS MAX
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_writer_DEB_ID_BLOC_STD    = 0,    
  PDM_writer_DEB_ID_BLOC_POLY2D = 1000000,
  PDM_writer_DEB_ID_BLOC_POLY3D = 2000000

} PDM_writer_deb_id_bloc_t;

/*============================================================================
 * Variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Stockage des objets cs
 *----------------------------------------------------------------------------*/

static PDM_writer_t **cs_tab = NULL; 

/*----------------------------------------------------------------------------
 * Taille du tableau cs
 *----------------------------------------------------------------------------*/

static int l_cs_tab = 0;

/*----------------------------------------------------------------------------
 * Nombre d'objets cs stockes dans cs_tab
 *----------------------------------------------------------------------------*/

static int n_cs_tab = 0;

/*----------------------------------------------------------------------------
 * Stockage des objets cs
 *----------------------------------------------------------------------------*/

static PDM_writer_fmt_t **fmt_tab = NULL; 

/*----------------------------------------------------------------------------
 * Taille du tableau cs
 *----------------------------------------------------------------------------*/

static int l_fmt_tab = 0;

/*----------------------------------------------------------------------------
 * Nombre d'objets cs stockes dans cs_tab
 *----------------------------------------------------------------------------*/

static int n_fmt_tab = 0;

/*----------------------------------------------------------------------------
 * Nombre d'objets cs stockes dans cs_tab
 *----------------------------------------------------------------------------*/

static const int n_intern_fmt = 1;


/*============================================================================
 * Definition des fonctions privees
 *============================================================================*/


/*----------------------------------------------------------------------------
 *
 * Som free
 *
 * parameters:
 *   a           <-- Premiere valeur
 *   b           <-- Seconde valeur                    
 *
 * return:
 *   max(a, b)
 *
 *----------------------------------------------------------------------------*/

static inline PDM_writer_som_t *
_som_free
(
 PDM_writer_som_t *som
)
{
  if (som->parent != NULL) {
    som->parent =_som_free (som->parent);
  }

  if (som->coords != NULL) {
    free (som->coords);
    som->coords = NULL;
  }

  free (som);
  
  return NULL;
}

/*----------------------------------------------------------------------------
 *
 * Fonction max
 *
 * parameters:
 *   a           <-- Premiere valeur
 *   b           <-- Seconde valeur                    
 *
 * return:
 *   max(a, b)
 *
 *----------------------------------------------------------------------------*/

static inline PDM_g_num_t
_lmax
(
 PDM_g_num_t a,
 PDM_g_num_t b
)
{
  return ((a) > (b) ? (a) : (b));
}

/*----------------------------------------------------------------------------
 *
 * Fonction max
 *
 * parameters:
 *   a           <-- Premiere valeur
 *   b           <-- Seconde valeur                    
 *
 * return:
 *   max(a, b)
 *
 *----------------------------------------------------------------------------*/

static inline double
_dmax
(
 double a,
 double b
)
{
  return ((a) > (b) ? (a) : (b));
}


/*----------------------------------------------------------------------------
 *
 * Produit vectoriel
 *
 * parameters:
 *   a           <-- Premier vecteur
 *   b           <-- Second  vecteur                   
 *   c           <-- Resultat a X b                    
 *
 *
 *----------------------------------------------------------------------------*/

static inline void
_p_vect
(
 const double a[3],
 const double b[3],
       double c[3]
)
{
  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = b[0] * a[2] - a[0] * b[2];
  c[2] = a[0] * b[1] - b[0] * a[1];
}


/*----------------------------------------------------------------------------
 *
 * Norme d'un vecteur
 *
 * parameters:
 *   a           <-- Vecteur
 *
 *----------------------------------------------------------------------------*/

static inline double
_p_scal 
(
 const double a[3], 
 const double b[3] 
)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/*----------------------------------------------------------------------------
 *
 * Norme d'un vecteur
 *
 * parameters:
 *   a           <-- Vecteur
 *
 *----------------------------------------------------------------------------*/

static inline double
_norme 
(
 const double a[3]
)
{
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet CS a partir de son identificateur
 *
 * parameters :
 *   id_cs           <-- Identificateur  
 *
 * return :
 *                   --> Objet CS  
 *
 *----------------------------------------------------------------------------*/

static PDM_writer_t *
_PDM_writer_get
(
const int  id_cs
)
{
  PDM_writer_t *cs = NULL;

  if (id_cs >= l_cs_tab) {
    fprintf(stderr, "Erreur PDM_writer_get : Identificateur de l'objet CS trop grand\n");
    abort();
  }

  cs = cs_tab[id_cs];

  if (cs == NULL) {
    fprintf(stderr, "Erreur PDM_writer_get : Identificateur de l'objet CS non defini\n");
    abort();
  }

  return cs;
}


/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet geom a partir de son identificateur
 *
 * parameters :
 *   id_gom          <-- Identificateur  
 *
 * return :
 *                   --> Objet CS  
 *
 *----------------------------------------------------------------------------*/

static PDM_writer_geom_t *
_geom_get
(
 const PDM_writer_t *cs,
 const int   id_geom
)
{
  PDM_writer_geom_t *geom = NULL;

  if (id_geom >= cs->l_geom_tab) {
    fprintf(stderr, "Erreur _geom_get : Identificateur de geometrie trop grand\n");
    abort();
  }
 
  geom = cs->geom_tab[id_geom];

  if (geom == NULL) {
    fprintf(stderr, "Erreur _geom_get : Identificateur de geometrie incorrect\n");
    abort();
  }

  return geom;
}


/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet var a partir de son identificateur
 *
 * parameters :
 *   id_gom          <-- Identificateur  
 *
 * return :
 *                   --> Objet CS  
 *
 *----------------------------------------------------------------------------*/

static PDM_writer_var_t *
_var_get
(
 const PDM_writer_t *cs,
 const int   id_var
)
{
  PDM_writer_var_t *var = NULL;

  if (id_var >= cs->l_var_tab) {
    fprintf(stderr, "Erreur _var_get : Identificateur de variable trop grand\n");
    abort();
  }
 
  var = cs->var_tab[id_var];

  if (var == NULL) {
    fprintf(stderr, "Erreur _var_get : Identificateur de variable incorrect\n");
    abort();
  }

  return var;
}


/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet bloc_std a partir de son identificateur
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *   id_bloc         <-- Identificateur  
 *
 * return :
 *                   --> Objet CS  
 *
 *----------------------------------------------------------------------------*/

static PDM_writer_bloc_std_t *
_bloc_std_get
(
 const PDM_writer_geom_t *geom,
 const int        id_bloc
)
{
  PDM_writer_bloc_std_t *bloc = NULL;

  if (id_bloc >= geom->l_blocs_std) {
    fprintf(stderr, "Erreur _bloc_get : Identificateur de bloc standard trop grand\n");
    abort();
  }
 
  bloc = geom->blocs_std[id_bloc];

  if (bloc == NULL) {
    fprintf(stderr, "Erreur _bloc_get : Identificateur de bloc standard incorrect\n");
    abort();
  }

  return bloc;
}


/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet bloc_poly2d a partir de son identificateur
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *   id_bloc         <-- Identificateur  
 *
 * return :
 *                   --> Objet CS  
 *
 *----------------------------------------------------------------------------*/

static PDM_writer_bloc_poly2d_t *
_bloc_poly2d_get
(
 const PDM_writer_geom_t *geom,
 const int        id_bloc
)
{
  PDM_writer_bloc_poly2d_t *bloc = NULL;

  int _id_bloc = id_bloc - PDM_writer_DEB_ID_BLOC_POLY2D;

  if (_id_bloc >= geom->l_blocs_poly2d) {
    fprintf(stderr, "Erreur _bloc_poly2d_get : Identificateur de bloc standard trop grand\n");
    abort();
  }
 
  bloc = geom->blocs_poly2d[_id_bloc];

  if (bloc == NULL) {
    fprintf(stderr, "Erreur _bloc_poly2d_get : Identificateur de bloc standard incorrect\n");
    abort();
  }

  return bloc;
}


/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet bloc_poly3d a partir de son identificateur
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *   id_bloc         <-- Identificateur  
 *
 * return :
 *                   --> Objet CS  
 *
 *----------------------------------------------------------------------------*/

static PDM_writer_bloc_poly3d_t *
_bloc_poly3d_get
(
 const PDM_writer_geom_t *geom,
 const int        id_bloc
)
{
  PDM_writer_bloc_poly3d_t *bloc = NULL;

  int _id_bloc = id_bloc - PDM_writer_DEB_ID_BLOC_POLY3D;

  if (_id_bloc >= geom->l_blocs_poly3d) {
    fprintf(stderr, "Erreur _bloc_poly3d_get : Identificateur de bloc standard trop grand\n");
    abort();
  }
 
  bloc = geom->blocs_poly3d[_id_bloc];

  if (bloc == NULL) {
    fprintf(stderr, "Erreur _bloc_poly3d_get : Identificateur de bloc standard incorrect\n");
    abort();
  }

  return bloc;
}


/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet CS a partir de son identificateur
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_geom_init
(
PDM_writer_geom_t *geom 
)
{
  geom->nom_geom= NULL;

  geom->st_decoup_poly2d = PDM_writer_OFF;
  geom->st_decoup_poly3d = PDM_writer_OFF;

  geom->n_som_abs = 0;          /* Nombre absolu de sommets */
  geom->n_som_abs_total = 0;    /* Nombre absolu de sommets 
                                   (sommets issus du decoupage 
                                   des polyedres/polygones compris) */
  geom->n_elt_abs = 0;          /* Nombre absolu d'elements */
  geom->n_part = 0;             /* Nombre de partitions */

  geom->som            = NULL;  /* Description des sommets de chaque partition */

  geom->n_blocs_std    = 0;     /* Nb de blocs d'elements standards */                     
  geom->l_blocs_std    = 0;     /* Taille du tableau des blocs d'elements standards */
  geom->blocs_std      = NULL;  /* Blocs d'elements standards de chaque partition */
  geom->n_blocs_poly2d = 0;     /* Nb de blocs de polygones */                              
  geom->l_blocs_poly2d = 0;     /* Taille du tableau des blocs de polygones */          
  geom->blocs_poly2d   = NULL;  /* Blocs de polygones de chaque partition */          
  geom->n_blocs_poly3d = 0;     /* Nb de blocs de polyedres */                              
  geom->l_blocs_poly3d = 0;     /* Taille du tableau des blocs de polyedres */
  geom->blocs_poly3d   = NULL;  /* Blocs de polyedres de chaque partition */          

  geom->geom_fmt       = NULL;
  geom->_cs            = NULL;
  geom->pdm_mpi_comm       = PDM_MPI_COMM_NULL;
}


/*----------------------------------------------------------------------------
 * Retourne un pointeur un objet CS a partir de son identificateur
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_var_init
(
PDM_writer_var_t *var
)
{
   var->nom_var    = NULL;          /* Nom de la geometrie */
  var->st_dep_tps = PDM_writer_OFF;        /* Variable en temps */
  var->dim        = PDM_WRITER_VAR_CSTE;   /* Dimension de la variable */
  var->loc        = PDM_WRITER_VAR_SOMMETS;/* Dimension de la variable */
  var->_val       = NULL;          /* Valeurs de la variable */          
  var->var_fmt    = NULL;          /* Description propre au format fmt */
  var->_cs        = NULL;
}


/*----------------------------------------------------------------------------
 * Calcul les centres faces d'ensemble de polygones 
 *
 * parameters :
 *   n_som           <-- Nombre de sommets
 *   coords          <-- Coordonnees des sommets
 *   n_elt           <-- Nombre d'elements
 *   connec_idx      <-- Index de connectivite des polygones
 *   connec          <-- Connectivite des polygones           
 *   centre_face     --> Coordonnees des centres faces
 *
 *----------------------------------------------------------------------------*/

static void   
_centre_face
(
const double *coords,
const int     n_elt,
const int    *connec_idx,
const int    *connec,
      double *centre_face
)
{
  
  const double deplacement_min = 1e-9;   /* Deplacement min */
  const int    n_iter_max      = 100;    /* Nombre d'iterations max */
  const double geom_eps_min    = 1e-30;  /* Valeur min du denominateur*/
  double       vecteur_surface[3];

  /* Boucle sur les polygones */            

  int ipoint = 0;
  for (int i = 0; i < n_elt; i++) {
  
    const int ideb      = connec_idx[i];
    const int ifin      = connec_idx[i+1];
    const int n_som_elt = ifin - ideb;

    if (n_som_elt > 4) {

      /* Initialisation par le barycentre des sommets */

      for (int k = 0; k < 3; k++) 
        centre_face[3*ipoint + k] = 0.;
      
      for (int j = ideb; j < ifin; j++) {
        const int idx_som = connec[j] - 1;  
        for (int k = 0; k < 3; k++) 
          centre_face[3*ipoint + k] += coords[3*idx_som + k];
      }
  
      for (int k = 0; k < 3; k++) 
        centre_face[3*ipoint + k] /= n_som_elt;
  
      /* Calcul iteratif */
  
      int n_iter = 0; /* Nombre d'iterations */
      while (1) {
    
        double deplacement[3] = {0, 0, 0}; /* Deplacement du centre cellule */

        n_iter += 1;
    
        for (int k1 = 0; k1 < 3; k1++)
          vecteur_surface[k1] = 0;
    
        double aire_poly = 0; 
    
        /* Boucle sur les aretes */
        
        for (int isom = 0; isom < n_som_elt; isom++) {
      
          const int som1 = connec[isom] - 1;  
          const int som2 = connec[(isom + 1) % n_som_elt] - 1;
        
          /* Centre des aretes */
      
          const double centre_arete[3] = {0.5 * (coords[3*som1    ] + coords[3*som2    ]),
                                          0.5 * (coords[3*som1 + 1] + coords[3*som2 + 1]),
                                          0.5 * (coords[3*som1 + 2] + coords[3*som2 + 2])};
      
          /* Vecteur arete */
          
          const double vect_arete[3] = {coords[3*som2    ] - coords[3*som1    ],
                                        coords[3*som2 + 1] - coords[3*som1 + 1],
                                        coords[3*som2 + 2] - coords[3*som1 + 2]};
      
          /* Vecteur centre face -> centre arete */
          
          const double vect_centres_face_arete[3] = {centre_arete[0] - centre_face[0],
                                                     centre_arete[1] - centre_face[1],
                                                     centre_arete[2] - centre_face[2]};
      
          /* Calcul des aires des triangles */
      
          double vecteur_surface_tri[3];
          _p_vect(vect_centres_face_arete, vect_arete, vecteur_surface_tri);

          for (int k1 = 0; k1 < 3; k1++) 
            vecteur_surface_tri[k1] *= 0.5;
      
          const double aire_tri = _norme(vecteur_surface_tri);

          aire_poly += aire_tri;
          for (int k1 = 0; k1 < 3; k1++) {
            vecteur_surface[k1] += vecteur_surface_tri[k1];
            deplacement[k1] += aire_tri * vect_centres_face_arete[k1];
          }

          aire_poly += aire_tri;
      
        }

        double denom_aire_poly = 1. / _dmax(fabs(aire_poly), geom_eps_min);
        
        for (int k = 0; k < 3; k++) {
          deplacement[k] = 2./3. * denom_aire_poly * deplacement[k];
          centre_face[3*ipoint +k] += deplacement[k];
        }
        
        /* Verification de la convergence */
        
        const double norme_deplacement = _norme(deplacement);
        
        if (norme_deplacement < deplacement_min) {
          break;
        }
        
        /* Verification du nombre d'iterations */  
        
        else if (n_iter_max < n_iter) {
          break;
        }
      } // while (1)
      ipoint += 1;
    }
  }

  if (1 == 1) {
    printf("centres : ");
    for (int i = 0; i < ipoint; i++) {
      printf("%12.5e %12.5e %12.5e\n", 
             centre_face[3*i  ],
             centre_face[3*i+1],
             centre_face[3*i+2]);
    }
    printf("\n");
  }
}


/*----------------------------------------------------------------------------
 * Calcule la numerotation aboslue des elements  
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_calcul_numabs_bloc
(
 const PDM_MPI_Comm    pdm_mpi_comm,
 const int         n_procs,
 const int         i_proc,
       int        *sendBuffN,
       int        *sendBuffIdx,
       int        *recvBuffN,
       int        *recvBuffIdx,
       PDM_g_num_t  *n_elt_stocke_procs,
 const int         l_numabs_tmp,
       PDM_g_num_t  *numabs_tmp,
       PDM_g_num_t  *d_elt_proc, 
 const PDM_l_num_t    n_part,
       PDM_l_num_t   *n_elt,
       PDM_g_num_t **numabs,
       PDM_g_num_t **numabs_int
)
{
  /* Calcul du nombre total d'elements du bloc */

  PDM_l_num_t n_elt_loc_total = 0;

  for (int j = 0; j < n_part; j++) {
    n_elt_loc_total += n_elt[j];
  }

  /* Comptage du nombre d'elements a envoyer a chaque processus */
  
  for (int j = 0; j < n_procs; j++) {
    sendBuffN[j] = 0;
    sendBuffIdx[j] = 0;
    recvBuffN[j] = 0;
    recvBuffIdx[j] = 0;
  }

  for (int j = 0; j < n_part; j++) {
    for (int k = 0; k < n_elt[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long (numabs[j][k],
                                                         d_elt_proc,
                                                         n_procs + 1);
      sendBuffN[i_elt_proc] += 1;
    }
  }

  sendBuffIdx[0] = 0;
  for (int j = 1; j < n_procs; j++) {
    sendBuffIdx[j] = sendBuffIdx[j-1] + sendBuffN[j-1];
  }
   
  /* Determination du nombre d'elements recu de chaque processus */

  PDM_MPI_Alltoall(sendBuffN, 
               1, 
               PDM_MPI_INT, 
               recvBuffN, 
               1, 
               PDM_MPI_INT, 
               pdm_mpi_comm);

  recvBuffIdx[0] = 0;
  for(int j = 1; j < n_procs; j++) {
    recvBuffIdx[j] = recvBuffIdx[j-1] + recvBuffN[j-1];
  }

  /* Transmission des numeros absolus  */
  
  PDM_g_num_t *sendBuffNumabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                                       n_elt_loc_total);
  PDM_g_num_t *recvBuffNumabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                                      (recvBuffIdx[n_procs - 1] + 
                                                       recvBuffN[n_procs - 1]));
                                                          
  for (int j = 0; j < n_procs; j++) {
    sendBuffN[j] = 0;
  }

  for (int j = 0; j < n_part; j++) {
    for (int k = 0; k < n_elt[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long(numabs[j][k],
                                                        d_elt_proc,
                                                        n_procs+1);
      sendBuffNumabs[sendBuffIdx[i_elt_proc] + sendBuffN[i_elt_proc]] = numabs[j][k];
      sendBuffN[i_elt_proc] += 1;
    }
  }

  PDM_MPI_Alltoallv((void *) sendBuffNumabs, 
                sendBuffN, 
                sendBuffIdx, 
                PDM__PDM_MPI_G_NUM,
                (void *) recvBuffNumabs, 
                recvBuffN, 
                recvBuffIdx,
                PDM__PDM_MPI_G_NUM, 
                pdm_mpi_comm);
  
  /* Echange du nombre d'elements stockes sur chaque processus */

  const PDM_g_num_t n_elt_stocke = 
    (PDM_g_num_t) (recvBuffIdx[n_procs - 1] + recvBuffN[n_procs - 1]);

  PDM_MPI_Allgather((void *) &n_elt_stocke, 
                1,
                PDM__PDM_MPI_G_NUM,
                (void *) (n_elt_stocke_procs + 1), 
                1,
                PDM__PDM_MPI_G_NUM,
                pdm_mpi_comm);

  n_elt_stocke_procs[0] = 1;
  for (int j = 1; j < n_procs + 1; j++) {
    n_elt_stocke_procs[j] += n_elt_stocke_procs[j-1];
  }    

  /* Stockage du resultat et determination de la nouvelle numerotation absolue
     independante du parallelisme */
    
  for (int j = 0; j < l_numabs_tmp; j++) 
    numabs_tmp[j] = 0;

  for (int j = 0; j < n_procs; j++) {
    
    const int ideb = recvBuffIdx[j];
    const int ifin = recvBuffIdx[j] + recvBuffN[j];
    
    for (int k = ideb; k < ifin; k++) {
      
      PDM_g_num_t _idx = recvBuffNumabs[k] - d_elt_proc[i_proc];
      const int idx = (int) _idx;
      assert((idx < l_numabs_tmp) && (idx >= 0));

      numabs_tmp[idx] = 1; /* On marque les elements */
    }
  }

  int cpt_elt_proc = 0;
  for (int j = 0; j < l_numabs_tmp; j++) {
    if (numabs_tmp[j] == 1) {

      /* On fournit une numerotation independante du parallelisme */
      
      numabs_tmp[j] = n_elt_stocke_procs[i_proc] + cpt_elt_proc;
      cpt_elt_proc += 1;
    }
  }

  /* On remplit le buffer de reception qui devient le buffer d'envoi
     Le buffer d'envoi devient lui le buffer de reception */

  cpt_elt_proc = 0;
  for (int j = 0; j < n_procs; j++) {
    
    const int ideb = recvBuffIdx[j];
    const int ifin = recvBuffIdx[j] + recvBuffN[j];
    
    for (int k = ideb; k < ifin; k++) {
      
      PDM_g_num_t _idx = recvBuffNumabs[k] - d_elt_proc[i_proc];
      const int idx = (int) _idx;
      
      recvBuffNumabs[cpt_elt_proc] = numabs_tmp[idx];
      
      cpt_elt_proc += 1;
    }
  }

  PDM_MPI_Alltoallv((void *) recvBuffNumabs, 
                recvBuffN, 
                recvBuffIdx, 
                PDM__PDM_MPI_G_NUM,
                (void *) sendBuffNumabs, 
                sendBuffN, 
                sendBuffIdx,
                PDM__PDM_MPI_G_NUM, 
                pdm_mpi_comm);

  /* On Stocke l'information recue */

  for (int j = 0; j < n_procs; j++)
    sendBuffN[j] = 0;

  for (int j = 0; j < n_part; j++) {
    
    numabs_int[j] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_elt[j]);

    for (int k = 0; k < n_elt[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long(numabs[j][k],
                                                        d_elt_proc,
                                                        n_procs+1);
      numabs_int[j][k] = sendBuffNumabs[sendBuffIdx[i_elt_proc] + sendBuffN[i_elt_proc]];
      sendBuffN[i_elt_proc] += 1;
    }
  }

  /* Liberation memoire */

  free(sendBuffNumabs);
  free(recvBuffNumabs);
  
}

/*----------------------------------------------------------------------------
 * Calcule la numerotation aboslue des elements  
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_calcul_numabs_elt
(
PDM_writer_geom_t *geom 
)
{
  /* Calcul de la taille du bloc du processus courant */
    
  int n_procs = 0;
  PDM_MPI_Comm_size(geom->pdm_mpi_comm,
                &n_procs);
  
  int i_proc = 0;
  PDM_MPI_Comm_rank(geom->pdm_mpi_comm,
                &i_proc);

  PDM_g_num_t *d_elt_proc = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_procs + 1));

  /* Calcul du nombre absolu d'�l�ments */

  /* Tri a l'aide d'un tableau defini par blocs continus
     repartis sur l'ensemble des processus */

  PDM_g_num_t div_entiere = geom->n_elt_abs / n_procs;
  PDM_g_num_t div_reste = geom->n_elt_abs % n_procs;

  d_elt_proc[0] = 1;
  for (int i = 0; i < n_procs; i++) {
    d_elt_proc[i+1] =  div_entiere;
    if (i < div_reste)
      d_elt_proc[i+1] += 1;
  }

  /* Calcul de la repartition des elements sur les processus */

  for (int j = 0; j < n_procs; j++) {
    d_elt_proc[j+1] += d_elt_proc[j] ;
  }

  PDM_g_num_t _l_numabs_tmp = d_elt_proc[i_proc+1] - d_elt_proc[i_proc];
  int l_numabs_tmp = (int) _l_numabs_tmp;
  PDM_g_num_t *numabs_tmp = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * l_numabs_tmp);

  /* Allocation des tableaux pour echanges MPI */

  PDM_g_num_t *n_elt_stocke_procs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_procs+1));

  int *sendBuffN   = (int *) malloc(sizeof(int) * n_procs);
  int *sendBuffIdx = (int *) malloc(sizeof(int) * n_procs);

  int *recvBuffN   = (int *) malloc(sizeof(int) * n_procs);
  int *recvBuffIdx = (int *) malloc(sizeof(int) * n_procs);

  /* Boucle sur les blocs standard */

  for (int i = 0; i < geom->n_blocs_std; i++) {
    PDM_writer_bloc_std_t *_bloc_std = geom->blocs_std[i];
    _calcul_numabs_bloc(geom->pdm_mpi_comm,
                        n_procs,
                        i_proc,
                        sendBuffN,
                        sendBuffIdx,
                        recvBuffN,
                         recvBuffIdx,
                        n_elt_stocke_procs,
                         l_numabs_tmp,
                        numabs_tmp,
                         d_elt_proc,
                        _bloc_std->n_part,
                        _bloc_std->n_elt,
                        _bloc_std->_numabs,
                        _bloc_std->numabs_int);
  }    

  /* Boucle sur les blocs de polygones s'ils ne sont pas decoupes 
     en triangle. Dans le cas contraire la numerotation absolue
     est calculee lors du decoupage */

  for (int i = 0; i < geom->n_blocs_poly2d; i++) {
    PDM_writer_bloc_poly2d_t *_bloc_poly2d = geom->blocs_poly2d[i];
    if (_bloc_poly2d->st_decoup_poly2d == PDM_writer_OFF) {
      _calcul_numabs_bloc(geom->pdm_mpi_comm,
                          n_procs,
                          i_proc,
                          sendBuffN,
                          sendBuffIdx,
                          recvBuffN,
                          recvBuffIdx,
                          n_elt_stocke_procs,
                          l_numabs_tmp,
                          numabs_tmp,
                          d_elt_proc,
                          _bloc_poly2d->n_part,
                          _bloc_poly2d->n_elt,
                          _bloc_poly2d->_numabs,
                          _bloc_poly2d->numabs_int);

    }
  }

  /* Boucle sur les blocs de polyedres  s'ils ne sont pas decoupes 
     en triangle. Dans le cas contraire la numerotation absolue
     est calculee lors du decoupage */

  for (int i = 0; i < geom->n_blocs_poly3d; i++) {
    PDM_writer_bloc_poly3d_t *_bloc_poly3d = geom->blocs_poly3d[i];
    if (_bloc_poly3d->st_decoup_poly3d == PDM_writer_OFF) {
      _calcul_numabs_bloc(geom->pdm_mpi_comm,
                          n_procs,
                          i_proc,
                          sendBuffN,
                          sendBuffIdx,
                          recvBuffN,
                          recvBuffIdx,
                          n_elt_stocke_procs,
                          l_numabs_tmp,
                          numabs_tmp,
                          d_elt_proc,
                          _bloc_poly3d->n_part,
                          _bloc_poly3d->n_elt,
                          _bloc_poly3d->_numabs,
                          _bloc_poly3d->numabs_int);
    }
  }

  /* Liberation memoire */

  free(d_elt_proc);
  free(numabs_tmp);
  free(n_elt_stocke_procs);
  free(sendBuffN);
  free(sendBuffIdx);
  free(recvBuffN);
  free(recvBuffIdx);

} 

/*----------------------------------------------------------------------------
 * Calcule la numerotation absolue des sommets supplementaires issues
 * du decoupage des polygones. Cette numerotation se base sur la numerotation
 * absolue unitiale des polygones pour que le resultat soit independant du
 * parallelisme.  
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_calcul_numabs_split_poly2d
(
 PDM_writer_geom_t *geom 
)
{

  /* Calcul de la taille du bloc du processus courant */
  
  int n_procs = 0;
  PDM_MPI_Comm_size(geom->pdm_mpi_comm,
                &n_procs);
  
  int i_proc = 0;
  PDM_MPI_Comm_rank(geom->pdm_mpi_comm,
                &i_proc);

  PDM_g_num_t *d_elt_proc = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_procs + 1));

  /* Tri a l'aide d'un tableau defini par blocs continus
     repartis sur l'ensemble des processus */
  
  PDM_g_num_t div_entiere = geom->n_elt_abs / n_procs;
  PDM_g_num_t div_reste = geom->n_elt_abs % n_procs;
  
  d_elt_proc[0] = 1;
  for (int i = 0; i < n_procs; i++) {
    d_elt_proc[i+1] =  div_entiere;
    if (i < div_reste)
      d_elt_proc[i+1] += 1;
  }
  
  PDM_g_num_t l_numabs_tmp = d_elt_proc[i_proc+1] - d_elt_proc[i_proc];
  PDM_g_num_t *numabs_poly_tmp = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * l_numabs_tmp);
  PDM_g_num_t *numabs_tri_tmp = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * l_numabs_tmp);

  /* Calcul de la repartition des elements sur les processus */

  for (int j = 0; j < n_procs; j++) {
    d_elt_proc[j+1] += d_elt_proc[j] ;
  }

  /* Nombre de donnees a echanger a travers le buffer */

  const int n_donnee_buff = 2;

  /* Allocation des tableaux pour echanges MPI */

  PDM_g_num_t *n_elt_stocke_procs = (PDM_g_num_t *) malloc(3*sizeof(PDM_g_num_t) * (n_procs+1)); 
  /* Nombre d'elements stock�s suivant leur type */

  int *sendBuffN   = (int *) malloc(sizeof(int) * n_procs);
  int *sendBuffIdx = (int *) malloc(sizeof(int) * n_procs);
  
  int *recvBuffN   = (int *) malloc(sizeof(int) * n_procs);
  int *recvBuffIdx = (int *) malloc(sizeof(int) * n_procs);


  for (int i = 0; i < geom->n_blocs_poly2d; i++) {
    
    PDM_writer_bloc_poly2d_t *_bloc_poly2d = geom->blocs_poly2d[i];

    /* Calcul du nombre total de vrais polygones du bloc */

    PDM_l_num_t n_elt_loc_total = 0;

    for (int j = 0; j < geom->n_part; j++) {
      for (int k = 0; k < _bloc_poly2d->n_elt[j]; k++) {
        n_elt_loc_total += 1;
      }
    }

    /* Comptage du nombre de sommets a envoyer a chaque processus
       (1 pour chaque vrai polygone) */
  
    for (int j = 0; j < n_procs; j++) {
      sendBuffN[j] = 0;
      sendBuffIdx[j] = 0;
      recvBuffN[j] = 0;
      recvBuffIdx[j] = 0;
    }

    sendBuffIdx[0] = 0;
    for (int j = 1; j < n_procs; j++) {
      sendBuffIdx[j] = sendBuffIdx[j-1] + sendBuffN[j-1];
    }
   
    /* Determination du nombre d'elements recu de chaque processus */

    PDM_MPI_Alltoall(sendBuffN, 
                 1, 
                 PDM_MPI_INT, 
                 recvBuffN, 
                 1, 
                 PDM_MPI_INT, 
                 geom->pdm_mpi_comm);

    recvBuffIdx[0] = 0;
    for(int j = 1; j < n_procs; j++) {
      recvBuffIdx[j] = recvBuffIdx[j-1] + recvBuffN[j-1];
    }
    
    /* Buffer pour les MPI_alltoall 2 donnees a stocker :
           - numero absolu du polygone
           - nombre de sommets du triangle */

    
    PDM_g_num_t *sendBuffNumabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                                     n_elt_loc_total * 
                                                     n_donnee_buff);
    PDM_g_num_t *recvBuffNumabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 
                                                     (recvBuffIdx[n_procs - 1] + recvBuffN[n_procs - 1]) *
                                                     n_donnee_buff);
    
    for (int j = 0; j < n_procs; j++) {
      sendBuffN[j] = 0;
    }
    
    for (int j = 0; j < geom->n_part; j++) {
      for (int k = 0; k < _bloc_poly2d->n_elt[j]; k++) {
        const int n_som_elt = _bloc_poly2d->_connec_idx[j][k+1] - _bloc_poly2d->_connec_idx[j][k];
        const int i_elt_proc = PDM_binary_search_gap_long(_bloc_poly2d->_numabs[j][k],
                                                          d_elt_proc,
                                                          n_procs+1);
        sendBuffNumabs[sendBuffIdx[i_elt_proc] * n_donnee_buff + sendBuffN[i_elt_proc]    ] =
          _bloc_poly2d->_numabs[j][k];
        sendBuffNumabs[sendBuffIdx[i_elt_proc] * n_donnee_buff + sendBuffN[i_elt_proc] + 1] = 
          n_som_elt;
        sendBuffN[i_elt_proc] += n_donnee_buff;
      }
    }

    PDM_MPI_Alltoallv((void *) sendBuffNumabs, 
                  sendBuffN, 
                  sendBuffIdx, 
                  PDM__PDM_MPI_G_NUM,
                  (void *) recvBuffNumabs, 
                  recvBuffN, 
                  recvBuffIdx,
                  PDM__PDM_MPI_G_NUM, 
                  geom->pdm_mpi_comm);

    /* Stockage du resultat et determination de la numerotation 
       absolue des vrais polygones, des triangles et quadrangles */

    for (int j = 0; j < l_numabs_tmp; j++) {
      numabs_poly_tmp[j] = 0;
      numabs_tri_tmp[j] = 0;
    }

    PDM_g_num_t n_elt_stocke[3] = {0, 0, 0};
    
    for (int j = 0; j < n_procs; j++) {
      
      const int ideb = recvBuffIdx[j];
      const int ifin = recvBuffIdx[j] + recvBuffN[j];
      
      for (int k = ideb; k < ifin; k++) {
      
        PDM_g_num_t _idx = recvBuffNumabs[k*n_donnee_buff] - d_elt_proc[i_proc];
        const int idx = (int) _idx;
      
        numabs_poly_tmp[idx] = 1; /* Numerotation absolue du polygone */
        numabs_tri_tmp[idx] = recvBuffNumabs[k*n_donnee_buff + 1];  /* Initialisation au nb de triangles */
        if (numabs_tri_tmp[idx] == 3)
          n_elt_stocke[0] += 1;
        else if (numabs_tri_tmp[idx] == 4)
          n_elt_stocke[1] += 1;
        else if (numabs_tri_tmp[idx] > 4)
          n_elt_stocke[2] += 1;
        
      }
    }
    
    /* On recupere le nombre d'elements stocke sur chaque processus suivant leur type */
    
    PDM_MPI_Allgather((void *) &n_elt_stocke, 
                  3,
                  PDM__PDM_MPI_G_NUM,
                  (void *) (n_elt_stocke_procs + 1), 
                  3,
                  PDM__PDM_MPI_G_NUM,
                  geom->pdm_mpi_comm);
    
    n_elt_stocke_procs[0] = 1;
    n_elt_stocke_procs[1] = 1;
    n_elt_stocke_procs[2] = 1;
    
    for (int j = 1; j < n_procs + 1; j++) {
      for (int k = 0; k < 3; k++) {
        n_elt_stocke_procs[3*j + k] = n_elt_stocke_procs[3*(j-1) + k];
      }
    }
    
    /* On fournit une numerotation independante du parallelisme pour les elements suivant leur type */
    
    int cpt_tri_proc = 0;
    int cpt_quad_proc = 0;
    int cpt_poly_proc = 0;
    for (int j = 0; j < l_numabs_tmp; j++) {
      if (numabs_poly_tmp[j] == 1) {
        if (numabs_tri_tmp[j] == 3) {
          numabs_tri_tmp[j] = n_elt_stocke_procs[3*i_proc] + cpt_tri_proc;
          numabs_poly_tmp[j] = -1; 
          cpt_tri_proc += 1; 
        }
        else if (numabs_tri_tmp[j] == 4) {
          numabs_tri_tmp[j] = n_elt_stocke_procs[3*i_proc + 1] + cpt_quad_proc;
          numabs_poly_tmp[j] = -1; 
          cpt_quad_proc += 1;
        }
        else if (numabs_tri_tmp[j] > 4) {
          const int n_som = (int) numabs_tri_tmp[j]; 
          /* On fournit le num absolu du premier triangle
             uniquement : les autres en sont deduits */
          numabs_tri_tmp[j] =  n_elt_stocke_procs[3*i_proc + 1] + cpt_tri_proc; 
          /* Numerotation absolue des polygones qui sont en
             fait la numerotation des sommets supplementaires */  
          numabs_poly_tmp[j] =  n_elt_stocke_procs[3*i_proc + 2] + cpt_poly_proc; 
          cpt_tri_proc += n_som; /* On incremente du nombre de triangles dans le polygones */
          cpt_poly_proc += 1;
        }
      }
    }

    /* On remplit le buffer de reception qui devient le buffer d'envoi
       Le buffer d'envoi devient lui le buffer de reception */

    int cpt_elt_proc = 0;
    for (int j = 0; j < n_procs; j++) {
      
      const int ideb = recvBuffIdx[j];
      const int ifin = recvBuffIdx[j] + recvBuffN[j];
      
      for (int k = ideb; k < ifin; k++) {
        
        PDM_g_num_t _idx = recvBuffNumabs[k] - d_elt_proc[i_proc];
        const int idx = (int) _idx;
        
        recvBuffNumabs[n_donnee_buff * cpt_elt_proc    ] = numabs_poly_tmp[idx];
        recvBuffNumabs[n_donnee_buff * cpt_elt_proc + 1] = numabs_tri_tmp[idx];
        
        cpt_elt_proc += 1;
      }
    }
    
    PDM_MPI_Alltoallv((void *) recvBuffNumabs, 
                  recvBuffN, 
                  recvBuffIdx, 
                  PDM__PDM_MPI_G_NUM,
                  (void *) sendBuffNumabs, 
                  sendBuffN, 
                  sendBuffIdx,
                  PDM__PDM_MPI_G_NUM, 
                  geom->pdm_mpi_comm);
    
    /* On Stocke l'information recue */
    
    for (int j = 0; j < n_procs; j++) 
      sendBuffN[j] = 0;
    
    for (int j = 0; j < geom->n_part; j++) {
      for (int k = 0; k < _bloc_poly2d->n_elt[j]; k++) {
        const int n_som_elt = _bloc_poly2d->_connec_idx[j][k+1] - _bloc_poly2d->_connec_idx[j][k];
        const int i_elt_proc = PDM_binary_search_gap_long(_bloc_poly2d->_numabs[j][k],
                                                          d_elt_proc,
                                                          n_procs+1);
        if (n_som_elt > 4) {
          
          _bloc_poly2d->som_sup->numabs[j][k] = geom->n_som_abs_total + 
            sendBuffNumabs[n_donnee_buff * (sendBuffIdx[i_elt_proc] + 
                                            sendBuffN[i_elt_proc])];
          
          const int i_tri = _bloc_poly2d->tri_idx[j][k];
          for (int k1 = 0; k1 < n_som_elt; k1++) {
            _bloc_poly2d->bloc_tri->_numabs[j][i_tri + k1] = 
              sendBuffNumabs[n_donnee_buff * (sendBuffIdx[i_elt_proc] + 
                                              sendBuffN[i_elt_proc]) + 1] + k1;
          }
          
        }
        else if (n_som_elt == 4) {
          const int i_quad = _bloc_poly2d->quad_idx[j][k];
          _bloc_poly2d->bloc_quad->_numabs[j][i_quad] = 
            sendBuffNumabs[n_donnee_buff * (sendBuffIdx[i_elt_proc] + 
                                            sendBuffN[i_elt_proc]) + 1];
          
        }
        else if (n_som_elt == 3) {
          const int i_tri = _bloc_poly2d->tri_idx[j][k];
          _bloc_poly2d->bloc_tri->_numabs[j][i_tri] = 
            sendBuffNumabs[n_donnee_buff * (sendBuffIdx[i_elt_proc] + 
                                            sendBuffN[i_elt_proc]) + 1];
        }
      }
    }

    /* Mise � jour du nombre de sommets absolus total */

    geom->n_som_abs_total +=  n_elt_stocke_procs[3*n_procs+2];
    
    /* Liberation memoire */

    free(sendBuffNumabs);
    free(recvBuffNumabs);

  }

  /* Liberation memoire */

  free(d_elt_proc);
  free(numabs_tri_tmp);
  free(numabs_poly_tmp);
  free(n_elt_stocke_procs);
  free(sendBuffN);
  free(sendBuffIdx);
  free(recvBuffN);
  free(recvBuffIdx);

}


/*----------------------------------------------------------------------------
 * Calcule la numerotation absolue des sommets supplementaires issues
 * du decoupage des polyedres. Cette numerotation se base sur la numerotation
 * absolue unitiale des polygones pour que le resultat soit independant du
 * parallelisme.  
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_calcul_numabs_split_poly3d
(
 PDM_writer_geom_t *geom 
)
{
  fprintf(stderr, "_calcul_numabs_split_poly3d : Not yet impelemented\n");
  geom;
  abort();
}

/*----------------------------------------------------------------------------
 * Decoupe les polygones quelconques en triangles
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_split_poly2d
(
PDM_writer_geom_t *geom 
)
{

  /* Boucle sur les blocs de polygones */

  for (int i = 0; i < geom->n_blocs_poly2d; i++) {
    PDM_writer_bloc_poly2d_t *_bloc_poly2d = geom->blocs_poly2d[i];
    
    if (_bloc_poly2d->st_decoup_poly2d == PDM_writer_ON) {
      
      /* Allocation des structures des sommets supplementaires */
      
      if (_bloc_poly2d->som_sup == NULL)
        _bloc_poly2d->som_sup = (PDM_writer_som_sup_t *) malloc(sizeof(PDM_writer_som_sup_t));
      else {
        fprintf(stderr, "Erreur _split_poly2d : Polygones deja decoupes\n");
        abort();
      }

      PDM_writer_som_sup_t *_som_sup = _bloc_poly2d->som_sup;
      
      _som_sup->n_part = _bloc_poly2d->n_part;
      _som_sup->coords = (double **) malloc(sizeof(double *) * _som_sup->n_part);
      _som_sup->n_som  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * _som_sup->n_part);
      _som_sup->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * _som_sup->n_part);
      
      /* Allocation des structures des triangles supplementaires */

      if (_bloc_poly2d->bloc_tri == NULL)
        _bloc_poly2d->bloc_tri = (PDM_writer_bloc_std_t *) malloc(sizeof(PDM_writer_bloc_std_t));
      else {
        fprintf(stderr, "Erreur _split_poly2d : Polygones deja decoupes\n");
        abort();
      }
      PDM_writer_bloc_std_t  *_bloc_tri = _bloc_poly2d->bloc_tri;
      
      _bloc_tri->t_elt   = PDM_writer_TRIA3;
      _bloc_tri->n_part  = _bloc_poly2d->n_part;
      _bloc_tri->n_elt   = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * _bloc_tri->n_part);
      _bloc_tri->_connec = NULL; /* pas de mapping memoire on met a nul */
      _bloc_tri->_numabs = NULL; /* pas de mapping memoire on met a nul */
      _bloc_tri->numabs_int = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * _bloc_tri->n_part);
      
      if (_bloc_poly2d->tri_idx == NULL)
        _bloc_poly2d->tri_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * _bloc_poly2d->n_part);
      else {
        fprintf(stderr, "Erreur _split_poly2d : Polygones deja decoupes\n");
        abort();
      }
      PDM_l_num_t  **_tri_idx = _bloc_poly2d->tri_idx;
      
      /* Allocation des structures des quadrangles */
      
      if (_bloc_poly2d->bloc_quad == NULL)
        _bloc_poly2d->bloc_quad = (PDM_writer_bloc_std_t *) malloc(sizeof(PDM_writer_bloc_std_t));
      else {
        fprintf(stderr, "Erreur _split_poly2d : Polygones deja decoupes\n");
        abort();
      }
      PDM_writer_bloc_std_t  *_bloc_quad = _bloc_poly2d->bloc_quad;
      
      _bloc_quad->t_elt   = PDM_writer_TRIA3;
      _bloc_quad->n_part  = _bloc_poly2d->n_part;
      _bloc_quad->n_elt   = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * _bloc_quad->n_part);
      _bloc_quad->_connec = NULL; /* pas de mapping memoire on met a nul */
      _bloc_quad->_numabs = NULL; /* pas de mapping memoire on met a nul */
      _bloc_quad->numabs_int = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * _bloc_quad->n_part);
      
      if (_bloc_poly2d->quad_idx == NULL)
        _bloc_poly2d->quad_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * _bloc_poly2d->n_part);
      else {
        fprintf(stderr, "Erreur _split_poly2d : Polygones deja decoupes\n");
        abort();
      }
      PDM_l_num_t  **_quad_idx = _bloc_poly2d->quad_idx;
      
      /* Boucle sur les partitions */

      int n_tri      = 0;
      int n_poly     = 0;
      int n_quad     = 0;
      int n_tri_poly = 0;
      
      for (int j = 0; j < _bloc_poly2d->n_part; j++) {

        /* Initialisation des tableaux sur les sommets 
           Le nombre de sommets ajoute est le nombre de polygone du bloc 
           moins les triangles et quadrangles qui ne seront pas redecoupes */

        for (int k = 0; k < _bloc_poly2d->n_elt[j]; k++) {
          const int n_som_elt = _bloc_poly2d->_connec_idx[j][k+1] - _bloc_poly2d->_connec_idx[j][k];
          if (n_som_elt > 4) {
            n_tri_poly += n_som_elt;
            n_poly += 1;
          }
          else if (n_som_elt == 4)
            n_quad += 1;
          else if (n_som_elt == 3)
            n_tri += 1;
        }
        
        _som_sup->n_som[j]  = n_poly;
        _som_sup->coords[j] = (double *)     malloc(3 * n_poly * sizeof(double));
        _som_sup->numabs[j] = (PDM_g_num_t *)  malloc(n_poly * sizeof(PDM_g_num_t));

        /* Calcul des centres faces et ajout dans la liste des sommets 
         * (les triangles et quadrangles sont ignores)
         *   Remarque : il faudrait reprendre la fonction geofac3
         *              comme elle est dans BIBCEDRE il faudrait la sortir et creer 
         *              une bibliotheque geometrique. La fonction est donc reecrite 
         */

        _centre_face(geom->som[j]->_coords,
                     _bloc_poly2d->n_elt[j],
                     _bloc_poly2d->_connec_idx[j],
                     _bloc_poly2d->_connec[j],
                     _som_sup->coords[j]);
        
        /* Initialisation des tableaux des triangles */

        _tri_idx[j]             = (PDM_l_num_t *) malloc((_bloc_poly2d->n_elt[j] + 1) * sizeof(PDM_l_num_t));
        _bloc_tri->n_elt[j]     = n_tri + n_tri_poly;
        _bloc_tri->numabs_int[j]= (PDM_g_num_t *) malloc(_bloc_tri->n_elt[j] * sizeof(PDM_g_num_t));
        
        /* Initialisation des tableaux des quadrangles */

        _quad_idx[j]             = (PDM_l_num_t *) malloc((_bloc_poly2d->n_elt[j] + 1) * sizeof(PDM_l_num_t));
        _bloc_quad->n_elt[j]     = n_quad;
        _bloc_quad->numabs_int[j]= (PDM_g_num_t *) malloc(_bloc_quad->n_elt[j] * sizeof(PDM_g_num_t));
        
        /* Remplissage des tableaux d'indirection des polygones vers les triangles et quadrangles */

        _tri_idx[j][0] = 0;
        _quad_idx[j][0] = 0;
        for (int k = 1; k < _bloc_poly2d->n_elt[j] + 1; k++) {
          const int n_som = _bloc_poly2d->_connec_idx[j][k] - 
            _bloc_poly2d->_connec_idx[j][k-1];
          if (n_som == 4) {
            _tri_idx[j][k] = _tri_idx[j][k-1];
            _quad_idx[j][k] = _quad_idx[j][k-1] + 1;
          }
          else if (n_som == 3) {
            _tri_idx[j][k] = _tri_idx[j][k-1] + 1;
            _quad_idx[j][k] = _quad_idx[j][k-1];
          }
          else if (n_som > 4) {
            _tri_idx[j][k] = _tri_idx[j][k-1] + n_som;
            _quad_idx[j][k] = _quad_idx[j][k-1];
          }
        }

        /* Remplissage des tableaux de connectivite des triangles et quadrangles */

        /* int idx_tri  = 0; */
        /* int idx_quad = 0; */
        /* int ipoly    = 0; */
        /* for (int k = 0; k < _bloc_poly2d->n_elt[j]; k++) { */
        /*   const int ideb = _bloc_poly2d->_connec_idx[j][k]; */
        /*   const int ifin = _bloc_poly2d->_connec_idx[j][k+1]; */
        /*   const int n_som = ifin - ideb; */
        /*   if (n_som == 4) { */
        /*     for (int k1 = ideb; k1 < ifin; k1++)  */
        /*       _bloc_quad->connec[j][idx_quad++] = _bloc_poly2d->_connec[j][k1]; */
        /*   } */
        /*   else if (n_som == 3) { */
        /*     for (int k1 = ideb; k1 < ifin; k1++)  */
        /*       _bloc_tri->connec[j][idx_tri++] = _bloc_poly2d->_connec[j][k1]; */
        /*   } */
        /*   else if (n_som > 4) { */
        /*     for (int k1 = ideb; k1 < ifin; k1++) { */
        /*       _bloc_tri->connec[j][idx_tri++] = _bloc_poly2d->_connec[j][k1]; */
        /*       if (k1 == (ifin - 1)) */
        /*         _bloc_tri->connec[j][idx_tri++] = _bloc_poly2d->_connec[j][ideb]; */
        /*       else */
        /*         _bloc_tri->connec[j][idx_tri++] = _bloc_poly2d->_connec[j][k1+1]; */
        /*       _bloc_tri->connec[j][idx_tri++] = -(ipoly + 1); /\* Le numero supplementaire est definit  */
        /*                                                          dans la numerotation des sommets suppl�mentaires : */
        /*                                                          on les repere par leur valeur negative *\/ */
        /*       ipoly += 1; */
        /*     } */
        /*   } */
        /* } */
      }
    }
  }
}
 

/*----------------------------------------------------------------------------
 * Decoupe les polyedres en un ensemble de tetraedres et de pyramides 
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void   
_split_poly3d
(
PDM_writer_geom_t *geom 
)
{
  /* Boucle sur les blocs de polygones */

  for (int i = 0; i < geom->n_blocs_poly3d; i++) {
    PDM_writer_bloc_poly3d_t *_bloc_poly3d = geom->blocs_poly3d[i];
    
    if (_bloc_poly3d->st_decoup_poly3d == PDM_writer_ON) {
       
      /* Allocation des structures des sommets supplementaires */
      
      if (_bloc_poly3d->som_sup == NULL)
        _bloc_poly3d->som_sup = (PDM_writer_som_sup_t *) malloc(sizeof(PDM_writer_som_sup_t));
      else {
        fprintf(stderr, "Erreur _split_poly3d : Polyedres deja decoupes\n");
        abort();
      }

      PDM_writer_som_sup_t *_som_sup = _bloc_poly3d->som_sup;
      
      _som_sup->n_part = _bloc_poly3d->n_part;
      _som_sup->coords = (double **) malloc(sizeof(double *) * _som_sup->n_part);
      _som_sup->n_som  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * _som_sup->n_part);
      _som_sup->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * _som_sup->n_part);

      /* Allocation des structures de tetraedres supplementaires */

      if (_bloc_poly3d->bloc_tetra == NULL)
        _bloc_poly3d->bloc_tetra = (PDM_writer_bloc_std_t *) malloc(sizeof(PDM_writer_bloc_std_t));
      else {
        fprintf(stderr, "Erreur _split_poly3d : Polyedres deja decoupes\n");
        abort();
      }
      PDM_writer_bloc_std_t  *_bloc_tetra = _bloc_poly3d->bloc_tetra;

      _bloc_tetra->t_elt   = PDM_writer_TETRA4;
      _bloc_tetra->n_part  = _bloc_poly3d->n_part;
      _bloc_tetra->n_elt   = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * _bloc_tetra->n_part);
      _bloc_tetra->_connec = NULL; /* pas de mapping memoire on met a nul */
      _bloc_tetra->_numabs = NULL; /* pas de mapping memoire on met a nul */
      _bloc_tetra->numabs_int = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * _bloc_tetra->n_part);

      if (_bloc_poly3d->tetra_idx == NULL)
        _bloc_poly3d->tetra_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * _bloc_poly3d->n_part);
      else {
        fprintf(stderr, "Erreur _split_poly3d : Polyedres deja decoupes\n");
        abort();
      }

      /* Allocation des structures de pyramides suppl�mentaires */

      if (_bloc_poly3d->bloc_pyra == NULL)
        _bloc_poly3d->bloc_pyra = (PDM_writer_bloc_std_t *) malloc(sizeof(PDM_writer_bloc_std_t));
      else {
        fprintf(stderr, "Erreur _split_poly3d : Polygones deja decoupes\n");
        abort();
      }
      PDM_writer_bloc_std_t  *_bloc_pyra = _bloc_poly3d->bloc_pyra;

      _bloc_pyra->t_elt   = PDM_writer_PYRAMID5;
      _bloc_pyra->n_part  = _bloc_poly3d->n_part;
      _bloc_pyra->n_elt   = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * _bloc_pyra->n_part);
      _bloc_pyra->_connec = NULL; /* pas de mapping memoire on met a nul */
      _bloc_pyra->_numabs = NULL; /* pas de mapping memoire on met a nul */
      _bloc_pyra->numabs_int = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * _bloc_pyra->n_part);

      if (_bloc_poly3d->pyra_idx == NULL)
        _bloc_poly3d->pyra_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * _bloc_poly3d->n_part);
      else {
        fprintf(stderr, "Erreur _split_poly3d : Polyedres deja decoupes\n");
        abort();
      }

      /* Calcul des centres faces */

      //TODO : Supprimer le decoupage en tetra

      /* const int n_face_elem = _bloc_poly3d->_cellfac_idx[i+1] - _bloc_poly3d->_cellfac_idx[i]; */

      /* double *_centre_face = (double *) malloc(_bloc_poly3d->n_face[i] *  */
      /*                                          sizeof(double)); */

      /* _centre_face (geom->som[j].n_som, */
      /*               geom->som[j].coords, */
      /*               _bloc_poly3d->n_face[i], */
      /*               _bloc_poly3d->_facsom_idx[i], */
      /*               _bloc_poly3d->_facsom[i], */
      /*               _centre_face); */

      /* Calcul des centres cellules */

      /* Fusion des sommets crees en frontiere des sous-domaines
         avec une table de hachage sur la somme des sommets */

  /* Calcul des centres cellules (pas de fusion de points)   
   *   Remarque : il faudrait reprendre la fonction geocel3
   *              comme elle est dans BIBCEDRE il faudrait la sortir et creer 
   *              une bibliotheque geometrique. La fonction est donc reecrite
   *              en C
   */
    }
  }
  fprintf(stderr, "_split_poly3d not yet implemented\n");
  abort();
}


/*----------------------------------------------------------------------------
 *
 * Determine la connectivite nodale d'un tetra
 *
 * parameters :
 *   geom                 <-- Description des sommets
 *   cell_som_tria        <-- Connectivite des 4 triangles
 *   connec_tetra_courant --> Connectivite nodale du tetra
 *
 *----------------------------------------------------------------------------*/

static void
_connec_tetra
(
PDM_writer_som_t *som,
PDM_l_num_t *cell_som_tria,
PDM_l_num_t connec_tetra_courant[]
)
{

  /* Initialisation */

  connec_tetra_courant[0] = cell_som_tria[0];
  connec_tetra_courant[1] = cell_som_tria[1];
  connec_tetra_courant[2] = cell_som_tria[2];

  for (int i = 3; i < 11; i++) {
    if ((cell_som_tria[i] != connec_tetra_courant[0]) &&
        (cell_som_tria[i] != connec_tetra_courant[1]) &&
        (cell_som_tria[i] != connec_tetra_courant[2]))
      connec_tetra_courant[3] = cell_som_tria[i];
  } 

  /* Orientation */

  const PDM_real_t *_coords = som->_coords;
  double v1[3];
  double v2[3];
  double n[3];

  for (int i = 0; i < 3; i++) {
    v1[i] = _coords[3*(cell_som_tria[1] - 1) + i] - _coords[3*(cell_som_tria[0] - 1) + i];
    v2[i] = _coords[3*(cell_som_tria[2] - 1) + i] - _coords[3*(cell_som_tria[0] - 1) + i];
  }

  _p_vect(v1, v2, n);
  double orient = _p_scal(v1, n);

  if (orient < 0) {
    connec_tetra_courant[0] = cell_som_tria[2];
    connec_tetra_courant[1] = cell_som_tria[1];
    connec_tetra_courant[2] = cell_som_tria[0];
  }
}


/*----------------------------------------------------------------------------
 *
 * Determine la connectivite nodale d'un prisme
 *
 * parameters :
 *   geom                 <-- Description des sommets
 *   cell_som_tria        <-- Connectivite des 2 triangles
 *   cell_som_quad        <-- Connectivite des 3 quadrangles
 *   connec_prism_courant --> Connectivite nodale du prisme
 *
 *----------------------------------------------------------------------------*/

static void
_connec_prism
(
PDM_writer_som_t *som,
PDM_l_num_t *cell_som_tria,
PDM_l_num_t *cell_som_quad,
PDM_l_num_t connec_prism_courant[]
)
{

  /* Initialisation */

  for (int i = 0; i < 6; i++)
    connec_prism_courant[i] = cell_som_tria[i];

  /* Orientation des faces */

  const PDM_real_t *_coords = som->_coords;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 3; j++) {
      int isom = connec_prism_courant[3*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += _coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 1.0/3.0;
    
    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;
    
    double v1[3];
    double v2[3];
    int isom3 = connec_prism_courant[3*i+2] - 1 ;
    int isom2 = connec_prism_courant[3*i+1] - 1;
    int isom1 = connec_prism_courant[3*i] - 1;

    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom2+k] - _coords[3*isom1+k]; 
      v2[k] = _coords[3*isom3+k] - _coords[3*isom1+k];
    } 
    _p_vect(v1, v2, n + 3*i);
  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  double orientation = _p_scal(cc, n);
  double orientation2 = _p_scal(cc, n+3);

  if (orientation < 0) {
    int tmp = connec_prism_courant[1];
    connec_prism_courant[1] = connec_prism_courant[2];
    connec_prism_courant[2] = tmp;
  } 

  if (orientation2 < 0) {
    int tmp = connec_prism_courant[4];
    connec_prism_courant[4] = connec_prism_courant[5];
    connec_prism_courant[5] = tmp;
  } 

  /* Permutation circulaire */

  int id1 = -1;
  for (int j = 0; j < 12; j++) {
    if (cell_som_quad[j] == connec_prism_courant[0]) {
      id1 = j;
      break;
    }
  }

  int id2 = (id1 / 4) * 4 + (id1 + 1) % 4;
  if ((cell_som_quad[id2] == connec_prism_courant[1]) ||
      (cell_som_quad[id2] == connec_prism_courant[2]))
    id2 =  (id1 / 4) * 4 + (id1 + 3) % 4;

  int id_deb;
  for (int j = 0; j < 3; j++) {
    if (cell_som_quad[id2] == connec_prism_courant[3+j]) {
      id_deb = j;
      break;
    }
  }

  int tmp[3];
  for (int j = 0; j < 3; j++)
    tmp[j] = connec_prism_courant[3+j];

  for (int j = 0; j < 3; j++) {
    int idx = (id_deb + j) % 3;
    connec_prism_courant[3+j] = tmp[idx];
  }

}

/*----------------------------------------------------------------------------
 *
 * Determine la connectivite nodale d'une pyramide
 *
 * parameters :
 *   geom                   <-- Description des sommets
 *   cell_som_tria          <-- Connectivite des 2 triangles
 *   cell_som_quad          <-- Connectivite des 3 quadrangles
 *   connec_pyramid_courant --> Connectivite nodale de la pyramide
 *
 *----------------------------------------------------------------------------*/

static void
_connec_pyramid
(
PDM_writer_som_t *som,
PDM_l_num_t *cell_som_tria,
PDM_l_num_t *cell_som_quad,
PDM_l_num_t connec_pyramid_courant[]
)
{

  /* Initialisation */

  connec_pyramid_courant[0] = cell_som_quad[0];
  connec_pyramid_courant[1] = cell_som_quad[1];
  connec_pyramid_courant[2] = cell_som_quad[2];
  connec_pyramid_courant[3] = cell_som_quad[3];

  for (int i = 0; i < 9; i++) {
    if ((cell_som_tria[i] != connec_pyramid_courant[0]) &&
        (cell_som_tria[i] != connec_pyramid_courant[1]) &&
        (cell_som_tria[i] != connec_pyramid_courant[2]) &&
        (cell_som_tria[i] != connec_pyramid_courant[3])) {
      connec_pyramid_courant[4] = cell_som_tria[i];
      break;
    }
  }

  /* Orientation */

  const PDM_real_t *_coords = som->_coords;

  double c[3];
  double n[3];

  for (int k = 0; k < 3; k++)
    c[k] = 0.;
  for (int j = 0; j < 4; j++) {
    int isom = connec_pyramid_courant[j] - 1;
    for (int k = 0; k < 3; k++)
      c[k] += _coords[3*isom+k];
  }
  for (int k = 0; k < 3; k++)
    c[k] *= 0.25;
    
  for (int k = 0; k < 3; k++)
    n[k] = 0.;
    
  for (int j = 0; j < 4; j++) {
    int isom = connec_pyramid_courant[j] - 1;
    int suiv = (j+1) % 4;
    int isom_suiv = connec_pyramid_courant[suiv] - 1;

    double v1[3];
    double v2[3];
    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom+k] -  c[k]; 
      v2[k] = _coords[3*isom_suiv+k] -  c[k]; 
    } 
      
    _p_vect(v1, v2, n);

  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = _coords[3*(connec_pyramid_courant[3] - 1) + k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_scal(cc, n);
  
  if (orientation < 0) {
    int tmp = connec_pyramid_courant[0];
    connec_pyramid_courant[0] = connec_pyramid_courant[3];
    connec_pyramid_courant[3] = tmp;
    tmp = connec_pyramid_courant[1];
    connec_pyramid_courant[1] = connec_pyramid_courant[2];
    connec_pyramid_courant[2] = tmp;
  } 

}


/*----------------------------------------------------------------------------
 * Determine la connectivite nodale d'un hexa
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void
_connec_hexa
(
PDM_writer_som_t *som,
PDM_l_num_t *cell_som_quad,
PDM_l_num_t connec_hexa_courant[]
)
{

  /* Initialisation */

  connec_hexa_courant[0] = cell_som_quad[0];
  connec_hexa_courant[1] = cell_som_quad[1];
  connec_hexa_courant[2] = cell_som_quad[2];
  connec_hexa_courant[3] = cell_som_quad[3];

  PDM_l_num_t face_contact[4];

  for (int i = 1; i < 6; i++) {
    int cpt = 0;
    for (int j = 0; j < 4; j++) {
      PDM_l_num_t som_courant = cell_som_quad[4*i+j];
      if ((som_courant != connec_hexa_courant[0]) &&
          (som_courant != connec_hexa_courant[1]) &&
          (som_courant != connec_hexa_courant[2]) &&
          (som_courant != connec_hexa_courant[3]))
        cpt += 1;
    }
    if (cpt == 4) {
      connec_hexa_courant[4] = cell_som_quad[4*i];
      connec_hexa_courant[5] = cell_som_quad[4*i+1];
      connec_hexa_courant[6] = cell_som_quad[4*i+2];
      connec_hexa_courant[7] = cell_som_quad[4*i+3];
    }
    if (cpt == 2) {
      face_contact[0] = cell_som_quad[4*i];
      face_contact[1] = cell_som_quad[4*i+1];
      face_contact[2] = cell_som_quad[4*i+2];
      face_contact[3] = cell_som_quad[4*i+3];
    }
  } 

  /* Calcul des centres et normales de la base et de la face opposee */

  const PDM_real_t *_coords = som->_coords;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 4; j++) {
      int isom = connec_hexa_courant[4*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += _coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 0.25;
    
    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;
    
    for (int j = 0; j < 4; j++) {
      int isom = connec_hexa_courant[4*i+j] - 1;
      int suiv = (j+1) % 4;
      int isom_suiv = connec_hexa_courant[4*i+suiv] - 1;

      double v1[3];
      double v2[3];
      for (int k = 0; k < 3; k++) {
        v1[k] = _coords[3*isom+k] -  c[3*i+k]; 
        v2[k] = _coords[3*isom_suiv+k] -  c[3*i+k]; 
      } 
      
      _p_vect(v1, v2, n + 3*i);

    }

  }
  
  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_scal(cc, n);
  double orientation2 = _p_scal(cc, n+3);

  if (orientation < 0) {
    int tmp = connec_hexa_courant[0];
    connec_hexa_courant[0] = connec_hexa_courant[3];
    connec_hexa_courant[3] = tmp;
    tmp = connec_hexa_courant[1];
    connec_hexa_courant[1] = connec_hexa_courant[2];
    connec_hexa_courant[2] = tmp;
  } 

  if (orientation2 < 0) {
    int tmp = connec_hexa_courant[4];
    connec_hexa_courant[4] = connec_hexa_courant[7];
    connec_hexa_courant[7] = tmp;
    tmp = connec_hexa_courant[5];
    connec_hexa_courant[5] = connec_hexa_courant[6];
    connec_hexa_courant[6] = tmp;
  } 

  /* Permutation circulaire eventuelle de la face sup */

  int id1 = -1;
  int k1 = -1;
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      if (face_contact[j] == connec_hexa_courant[k]) {
        id1 = j;
        k1 = k;
        break;
      }
      if (id1 != -1)
        break;
    }
  }

  if (k1 == -1) {
    printf("Error connect_hexa : %d %d %d %d %d %d %d %d\n",
           connec_hexa_courant[0],
           connec_hexa_courant[1],
           connec_hexa_courant[2],
           connec_hexa_courant[3],
           connec_hexa_courant[4],
           connec_hexa_courant[5],
           connec_hexa_courant[6],
           connec_hexa_courant[7]);

    for (int i10 = 0; i10 < 4; i10++) {
      printf("   face %d : %d %d %d %d\n", i10+1, cell_som_quad[4*i10],  
                                                  cell_som_quad[4*i10+1],
                                                  cell_som_quad[4*i10+2],
                                                  cell_som_quad[4*i10+3]);
    }
    abort();
    
  }
    
  int id2 = (id1 + 1) % 4;
  int k2 = (k1 + 1) % 4;
  int k3 = (k1 + 3) % 4;
  
  if ((face_contact[id2] == connec_hexa_courant[k2]) ||
      (face_contact[id2] == connec_hexa_courant[k3]))
    id2 = (id1 + 3) % 4;

  int id_deb;
  for (int j = 0; j < 4; j++) {
    if (face_contact[id2] == connec_hexa_courant[4+j]) {
      id_deb = (j - k1);
      if (id_deb < 0) 
        id_deb += 4;
      id_deb = id_deb % 4;
      break;
    }
  }

  int tmp[4];
  for (int j = 0; j < 4; j++)
    tmp[j] = connec_hexa_courant[4+j];

  for (int j = 0; j < 4; j++) {
    int idx = (id_deb + j) % 4;
    connec_hexa_courant[4+j] = tmp[idx];
  }
}


/*----------------------------------------------------------------------------
 * Lib�re partiellement un bloc standard (conserve les numerotation absolues)
 *
 * parameters :
 *   bloc            <-- Bloc a liberer
 *
 *----------------------------------------------------------------------------*/

static
void
_bloc_std_free_partial
(
PDM_writer_bloc_std_t *_bloc_std
)
{
  if (_bloc_std->_connec != NULL) {
    if (_bloc_std->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_std->n_part; i++) {
        if (_bloc_std->_connec[i] != NULL)
          free(_bloc_std->_connec[i]);
        _bloc_std->_connec[i] = NULL;
      }
    }
    free(_bloc_std->_connec);
    _bloc_std->_connec = NULL;
  }
  
  if (_bloc_std->_numabs != NULL) {
    if (_bloc_std->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_std->n_part; i++) {
        if (_bloc_std->_numabs[i] != NULL)
          free(_bloc_std->_numabs[i]);
        _bloc_std->_numabs[i] = NULL;
      }
    }
    free(_bloc_std->_numabs);
    _bloc_std->_numabs = NULL;
  }
  
  if (_bloc_std->_num_part != NULL) {
    if (_bloc_std->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_std->n_part; i++) {
        if (_bloc_std->_num_part[i] != NULL)
          free(_bloc_std->_num_part[i]);
        _bloc_std->_num_part[i] = NULL;
      }
    }
    free(_bloc_std->_num_part);
    _bloc_std->_num_part = NULL;
  }
}


/*----------------------------------------------------------------------------
 * Lib�re un bloc standard 
 *
 * parameters :
 *   bloc            <-- Bloc a liberer
 *
 *----------------------------------------------------------------------------*/

static
void
_bloc_std_free
(
PDM_writer_bloc_std_t *_bloc_std
)
{
  _bloc_std_free_partial(_bloc_std);

  if (_bloc_std->n_elt != NULL) {
    free(_bloc_std->n_elt);
    _bloc_std->n_elt = NULL;
  }

  if (_bloc_std->numabs_int != NULL) {
    for (int j = 0; j < _bloc_std->n_part; j++) {
      if (_bloc_std->numabs_int[j] != NULL) {
        free(_bloc_std->numabs_int[j]);
      }
    }
    free(_bloc_std->numabs_int);
    _bloc_std->numabs_int = NULL;
  }

  free(_bloc_std);
}


/*----------------------------------------------------------------------------
 * Lib�re partiellement des sommets suppl�mentaires (conserve les numerotaion absolues)
 *
 * parameters :
 *   bloc            <-> Sommets � lib�rer
 *
 *----------------------------------------------------------------------------*/

static
void
_som_sup_free_partial
(
PDM_writer_som_sup_t *_som_sup
)
{
  if (_som_sup != NULL) {
    if (_som_sup->coords != NULL) {
      for (int j = 0; j < _som_sup->n_part; j++) {
        if (_som_sup->coords[j] != NULL) {
          free(_som_sup->coords[j]);
        }
      }
      free(_som_sup->coords);
      _som_sup->coords = NULL;
    }

    if (_som_sup->n_som != NULL) {
      free(_som_sup->n_som);
      _som_sup->n_som = NULL;
    }

  }
}


/*----------------------------------------------------------------------------
 * Lib�re les sommets suppl�mentaires 
 *
 * parameters :
 *   bloc            <-> Sommets � lib�rer
 *
 *----------------------------------------------------------------------------*/

static
void
_som_sup_free
(
PDM_writer_som_sup_t *_som_sup
)
{
  _som_sup_free_partial(_som_sup);

  if (_som_sup != NULL) {
    
    if (_som_sup->numabs != NULL) {
      for (int j = 0; j < _som_sup->n_part; j++) {
        if (_som_sup->numabs[j] != NULL) {
          free(_som_sup->numabs[j]);
        }
      }
      free(_som_sup->numabs);
    }
    
    free(_som_sup);

  }
}


/*----------------------------------------------------------------------------
 * Lib�re partiellement un bloc de polygones (conserve les numerotaion absolues)
 *
 * parameters :
 *   bloc            <-- Bloc a liberer
 *
 *----------------------------------------------------------------------------*/

static
void
_bloc_poly2d_free_partial
(
PDM_writer_bloc_poly2d_t *_bloc_poly2d
)
{

  if (_bloc_poly2d->_connec_idx != NULL) {
    if (_bloc_poly2d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly2d->n_part; i++) {
        if (_bloc_poly2d->_connec_idx[i] != NULL)
          free(_bloc_poly2d->_connec_idx[i]);
        _bloc_poly2d->_connec_idx[i] = NULL;
      }
    }
    free(_bloc_poly2d->_connec_idx);
    _bloc_poly2d->_connec_idx = NULL;
  }

  if (_bloc_poly2d->_connec != NULL) {
    if (_bloc_poly2d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly2d->n_part; i++) {
        if (_bloc_poly2d->_connec[i] != NULL)
          free(_bloc_poly2d->_connec[i]);
        _bloc_poly2d->_connec[i] = NULL;
      }
    }
    free(_bloc_poly2d->_connec);
    _bloc_poly2d->_connec = NULL;
  }
  
  if (_bloc_poly2d->_num_part != NULL) {
    if (_bloc_poly2d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly2d->n_part; i++) {
        if (_bloc_poly2d->_num_part[i] != NULL)
          free(_bloc_poly2d->_num_part[i]);
        _bloc_poly2d->_num_part[i] = NULL;
      }
    }
    free(_bloc_poly2d->_num_part);
    _bloc_poly2d->_num_part = NULL;
  }
  
  if (_bloc_poly2d->_numabs != NULL) {
    if (_bloc_poly2d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly2d->n_part; i++) {
        if (_bloc_poly2d->_numabs[i] != NULL)
          free(_bloc_poly2d->_numabs[i]);
        _bloc_poly2d->_numabs[i] = NULL;
      }
    }
    free(_bloc_poly2d->_numabs);
    _bloc_poly2d->_numabs = NULL;
  }

  if (_bloc_poly2d->bloc_tri != NULL) 
    _bloc_std_free_partial(_bloc_poly2d->bloc_tri);

  if (_bloc_poly2d->bloc_quad != NULL)
    _bloc_std_free_partial(_bloc_poly2d->bloc_quad);

  if (_bloc_poly2d->som_sup != NULL) {
    _som_sup_free_partial(_bloc_poly2d->som_sup);
  }

}


/*-----------------------------------------------------------------------------
 * Lib�re un bloc de polygones 
 *
 * parameters :
 *   bloc            <-- Bloc a liberer
 *
 *----------------------------------------------------------------------------*/

static
void
_bloc_poly2d_free
(
PDM_writer_bloc_poly2d_t *_bloc_poly2d
)
{
  _bloc_poly2d_free_partial(_bloc_poly2d);
  
  if (_bloc_poly2d->n_elt != NULL) {
    free(_bloc_poly2d->n_elt);
    _bloc_poly2d->n_elt = NULL;
  }

  if (_bloc_poly2d->numabs_int != NULL) {
    for (int j = 0; j < _bloc_poly2d->n_part; j++) {
      if (_bloc_poly2d->numabs_int[j] != NULL) {
        free(_bloc_poly2d->numabs_int[j]);
      }
    }
    free(_bloc_poly2d->numabs_int);
    _bloc_poly2d->numabs_int = NULL;
  }
  
  if (_bloc_poly2d->tri_idx != NULL) {
    for (int j = 0; j < _bloc_poly2d->n_part; j++) {
      if (_bloc_poly2d->tri_idx[j] != NULL) {
        free(_bloc_poly2d->tri_idx[j]);
      }
    }
    free(_bloc_poly2d->tri_idx);
    _bloc_poly2d->tri_idx = NULL;
  }

  if (_bloc_poly2d->bloc_tri != NULL) {
    _bloc_std_free(_bloc_poly2d->bloc_tri);
    _bloc_poly2d->bloc_tri = NULL;
  }

  if (_bloc_poly2d->quad_idx != NULL) {
    for (int j = 0; j < _bloc_poly2d->n_part; j++) {
      if (_bloc_poly2d->quad_idx[j] != NULL) {
        free(_bloc_poly2d->quad_idx[j]);
      }
    }
    free(_bloc_poly2d->quad_idx);
    _bloc_poly2d->quad_idx = NULL;
  }

  if (_bloc_poly2d->bloc_quad != NULL) {
    _bloc_std_free(_bloc_poly2d->bloc_quad);
    _bloc_poly2d->bloc_quad = NULL;
  }

  if (_bloc_poly2d->som_sup != NULL) {
    _som_sup_free(_bloc_poly2d->som_sup);
    _bloc_poly2d->som_sup = NULL;
  }

  free(_bloc_poly2d);
}


/*----------------------------------------------------------------------------
 * Lib�re partiellement un bloc de polygones (conserve les numerotaion absolues)
 *
 * parameters :
 *   bloc            <-- Bloc a liberer
 *
 *----------------------------------------------------------------------------*/

static
void
_bloc_poly3d_free_partial
(
PDM_writer_bloc_poly3d_t *_bloc_poly3d
)
{
  
  if (_bloc_poly3d->_facsom_idx != NULL) {
    if (_bloc_poly3d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly3d->n_part; i++) {
        if (_bloc_poly3d->_facsom_idx[i] != NULL)
          free(_bloc_poly3d->_facsom_idx[i]);
        _bloc_poly3d->_facsom_idx[i] = NULL;
      }
    }
    free(_bloc_poly3d->_facsom_idx);
    _bloc_poly3d->_facsom_idx = NULL;
  }
  
  if (_bloc_poly3d->_facsom != NULL) {
    if (_bloc_poly3d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly3d->n_part; i++) {
        if (_bloc_poly3d->_facsom[i] != NULL)
          free(_bloc_poly3d->_facsom[i]);
        _bloc_poly3d->_facsom[i] = NULL;
      }
    }
    free(_bloc_poly3d->_facsom);
    _bloc_poly3d->_facsom = NULL;
  }
  
  if (_bloc_poly3d->_cellfac_idx != NULL) {
    if (_bloc_poly3d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly3d->n_part; i++) {
        if (_bloc_poly3d->_cellfac_idx[i] != NULL)
          free(_bloc_poly3d->_cellfac_idx[i]);
        _bloc_poly3d->_cellfac_idx[i] = NULL;
      }
    }
    free(_bloc_poly3d->_cellfac_idx);
    _bloc_poly3d->_cellfac_idx = NULL;
  }
  
  if (_bloc_poly3d->_cellfac != NULL) {
    if (_bloc_poly3d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly3d->n_part; i++) {
        if (_bloc_poly3d->_cellfac[i] != NULL)
          free(_bloc_poly3d->_cellfac[i]);
        _bloc_poly3d->_cellfac[i] = NULL;
      }
    }
    free(_bloc_poly3d->_cellfac);
    _bloc_poly3d->_cellfac = NULL;
  }
  
  /* if (_bloc_poly3d->_num_part != NULL) { */
  /*   if (_bloc_poly3d->st_free_data == PDM_writer_ON) { */
  /*     for (int i = 0; i < _bloc_poly3d->n_part; i++) { */
  /*       if (_bloc_poly3d->_num_part[i] != NULL) */
  /*         free(_bloc_poly3d->_num_part[i]); */
  /*       _bloc_poly3d->_num_part[i] = NULL; */
  /*     } */
  /*   } */
  /*   free(_bloc_poly3d->_num_part); */
  /*   _bloc_poly3d->_num_part = NULL; */
  /* } */
  
  if (_bloc_poly3d->_numabs != NULL) {
    if (_bloc_poly3d->st_free_data == PDM_writer_ON) {
      for (int i = 0; i < _bloc_poly3d->n_part; i++) {
        if (_bloc_poly3d->_numabs[i] != NULL)
          free(_bloc_poly3d->_numabs[i]);
        _bloc_poly3d->_numabs[i] = NULL;
      }
    }
    free(_bloc_poly3d->_numabs);
    _bloc_poly3d->_numabs = NULL;
  }

  if (_bloc_poly3d->bloc_tetra != NULL) {
    _bloc_std_free_partial(_bloc_poly3d->bloc_tetra);
  }

  if (_bloc_poly3d->bloc_pyra != NULL) {
    _bloc_std_free_partial(_bloc_poly3d->bloc_pyra);
  }

  if (_bloc_poly3d->som_sup != NULL) {
    _som_sup_free_partial(_bloc_poly3d->som_sup);
  }
}


/*----------------------------------------------------------------------------
 * Lib�re un bloc de poly�dres 
 *
 * parameters :
 *   bloc            <-- Bloc a lib�rer
 *
 *----------------------------------------------------------------------------*/

static
void
_bloc_poly3d_free
(
PDM_writer_bloc_poly3d_t *_bloc_poly3d
)
{
  _bloc_poly3d_free_partial(_bloc_poly3d);

  if (_bloc_poly3d->n_elt != NULL) {
    free(_bloc_poly3d->n_elt);
    _bloc_poly3d->n_elt = NULL;
  }

  if (_bloc_poly3d->n_face!= NULL) {
    free(_bloc_poly3d->n_face);
    _bloc_poly3d->n_face= NULL;
  }

  if (_bloc_poly3d->numabs_int != NULL) {
    for (int j = 0; j < _bloc_poly3d->n_part; j++) {
      if (_bloc_poly3d->numabs_int[j] != NULL) {
        free(_bloc_poly3d->numabs_int[j]);
      }
    }
    free(_bloc_poly3d->numabs_int);
    _bloc_poly3d->numabs_int = NULL;
  }

  if (_bloc_poly3d->tetra_idx != NULL) {
    for (int j = 0; j < _bloc_poly3d->n_part; j++) {
      if (_bloc_poly3d->tetra_idx[j] != NULL) {
        free(_bloc_poly3d->tetra_idx[j]);
      }
    }
    free(_bloc_poly3d->tetra_idx);
    _bloc_poly3d->tetra_idx = NULL;
  }

  if (_bloc_poly3d->bloc_tetra != NULL) {
    _bloc_std_free(_bloc_poly3d->bloc_tetra);
    _bloc_poly3d->bloc_tetra = NULL;
  }

  if (_bloc_poly3d->pyra_idx != NULL) {
    for (int j = 0; j < _bloc_poly3d->n_part; j++) {
      if (_bloc_poly3d->pyra_idx[j] != NULL) {
        free(_bloc_poly3d->pyra_idx[j]);
      }
    }
    free(_bloc_poly3d->pyra_idx);
    _bloc_poly3d->pyra_idx = NULL;
  }

  if (_bloc_poly3d->bloc_pyra != NULL) {
    _bloc_std_free(_bloc_poly3d->bloc_pyra);
    _bloc_poly3d->bloc_pyra = NULL;
  }

  if (_bloc_poly3d->som_sup != NULL) {
    _som_sup_free(_bloc_poly3d->som_sup);
    _bloc_poly3d->som_sup = NULL;
  }

  free(_bloc_poly3d);
}


/*----------------------------------------------------------------------------
 * 
 * Type d'une cellule 3D
 *
 * parameters :
 *   bloc            <-- Bloc a lib�rer
 *
 *----------------------------------------------------------------------------*/

inline static
PDM_writer_elt_geom_t
_type_cell_3D
(
 const int          n_face_cell,
 const PDM_l_num_t    *cell_face,
 const PDM_l_num_t    *face_som_idx,
 const PDM_l_num_t    *face_som_nb,
 const PDM_l_num_t    *face_som,
 PDM_l_num_t           cell_som_tria[],
 PDM_l_num_t           cell_som_quad[]
)
{
  PDM_l_num_t  n_trias = 0;
  PDM_l_num_t  n_quads = 0;

  if (n_face_cell > 6) {
    return PDM_writer_POLY_3D;
  }

  for (int i = 0; i < n_face_cell; i++) {

    const int face_id = cell_face[i] - 1;
    const int n_som_face = face_som_nb[face_id];
    PDM_l_num_t idx = face_som_idx[face_id] - 1;
 
    if (n_som_face == 3) {
      PDM_l_num_t *cell_som_tria_courant = cell_som_tria + 3*n_trias;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_tria_courant[j-idx] = face_som[j];
      }
      n_trias += 1;
    }
    else if (n_som_face == 4) {
      PDM_l_num_t *cell_som_quad_courant = cell_som_quad + 4*n_quads;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_quad_courant[j-idx] = face_som[j];
      }
      n_quads += 1;
    }
    else 
      return PDM_writer_POLY_3D;

  }
      
  PDM_writer_elt_geom_t cell_type;

  if ((n_quads == 0) && (n_trias == 4))
    cell_type = PDM_writer_TETRA4;
  else if (n_quads == 6)
    cell_type = PDM_writer_HEXA8;
  else if ((n_quads == 1) && (n_trias == 4))
    cell_type = PDM_writer_PYRAMID5;
  else if ((n_quads == 3) && (n_trias == 2)) {
    int trias[6];
    n_trias = 0;
    for (int i = 0; i < n_face_cell; i++) {

      const int face_id = cell_face[i] - 1;
      const int ideb = face_som_idx[face_id] - 1;
      const int n_som_face = face_som_idx[face_id+1] - ideb - 1;
 
      if (n_som_face == 3) {
        for (int j = 0; j < 3; j++) {
          trias[3*n_trias+j] = face_som[ideb+j];
        }
        n_trias += 1;
      }
      if (n_trias >= 2)
        break;
    }

    cell_type = PDM_writer_PRISM6;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (trias[i] == trias[3+j]) {
          cell_type = PDM_writer_POLY_3D;
          break;
        }
      }
      if (cell_type == PDM_writer_POLY_3D)
        break;
    }
  }

  else {
    cell_type = PDM_writer_POLY_3D;
  }

  return cell_type;

}


/*----------------------------------------------------------------------------
 * 
 * Type d'une cellule 3D
 *
 * parameters :
 *   bloc            <-- Bloc a lib�rer
 *
 *----------------------------------------------------------------------------*/

static void
_load_intern_fmt (void)
{
  if (fmt_tab != NULL) {
    return;
  }

  n_fmt_tab = n_intern_fmt;
  l_fmt_tab = 2*n_fmt_tab;
  
  fmt_tab = malloc (sizeof(PDM_writer_fmt_t *) * l_fmt_tab);

  /* Ensight */

  fmt_tab[0] = malloc (sizeof(PDM_writer_fmt_t));
  fmt_tab[0]->name = "Ensight";
  fmt_tab[0]->create_fct       = PDM_writer_ensight_create;
  fmt_tab[0]->free_fct         = PDM_writer_ensight_free;
  fmt_tab[0]->beg_step_fct     = PDM_writer_ensight_step_beg;
  fmt_tab[0]->end_step_fct     = PDM_writer_ensight_step_end;
  fmt_tab[0]->geom_create_fct  = PDM_writer_ensight_geom_create;
  fmt_tab[0]->geom_write_fct   = PDM_writer_ensight_geom_write;
  fmt_tab[0]->geom_free_fct    = PDM_writer_ensight_geom_free;
  fmt_tab[0]->var_create_fct   = PDM_writer_ensight_var_create;
  fmt_tab[0]->var_write_fct    = PDM_writer_ensight_var_write;
  fmt_tab[0]->var_free_fct     = PDM_writer_ensight_var_free;
  
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Cree un objet CS (Cedre Sortie) et retoure un pointeur sur cet objet 
 *
 * parameters :
 *   fmt             <-- Format de sortie
 *   topologie       <-- Indique le maillage est mobile ou non
 *   st_reprise      <-- Complete les sorties des calcul precedents en reprise
 *   rep_sortie      <-- Repertoire de sortie                  
 *   nom_sortie      <-- Nom de la sortie                       
 *   pdm_mpi_com         <-- Communicateur MSG                      
 *
 * return :
 *                   --> Identificateur de l'objet cree
 *
 *----------------------------------------------------------------------------*/

void 
PROCF (pdm_writer_create_cf, PDM_WRITER_CREATE_CF)
(
const char          *fmt,
const int           *l_fmt,
const int           *fmt_fic,
const int           *topologie,
const int           *st_reprise,
const char          *rep_sortie,
const char          *nom_sortie,
const int           *l_rep_sortie,
const int           *l_nom_sortie,
const PDM_MPI_Fint      *pdm_mpi_comm,   
const int           *acces,
const double        *prop_noeuds_actifs,
int                 *id_cs     
ARGF_SUPP_CHAINE
)
{
  char *rep_sortie_c        = PDM_fortran_to_c_string(rep_sortie, *l_rep_sortie);
  char *nom_sortie_c        = PDM_fortran_to_c_string(nom_sortie, *l_nom_sortie);
  char *fmt_c        = PDM_fortran_to_c_string(fmt, *l_fmt);

  const PDM_MPI_Comm pdm_mpi_comm_c = PDM_MPI_Comm_f2c(*pdm_mpi_comm);

  *id_cs = PDM_writer_create(fmt_c,
                             (PDM_writer_fmt_fic_t)   *fmt_fic,
                             (PDM_writer_topologie_t) *topologie,
                             (PDM_writer_statut_t)    *st_reprise,
                             rep_sortie_c,
                             nom_sortie_c,
                             pdm_mpi_comm_c,  
                             (PDM_io_acces_t) *acces,
                             *prop_noeuds_actifs);

  free(rep_sortie_c);
  free(nom_sortie_c);
  free(fmt_c);
}

int
PDM_writer_create
(
const char       *fmt,
const PDM_writer_fmt_fic_t   fmt_fic,   
const PDM_writer_topologie_t topologie,
const PDM_writer_statut_t    st_reprise,
const char          *rep_sortie,
const char          *nom_sortie,
const PDM_MPI_Comm       pdm_mpi_comm,
const PDM_io_acces_t acces,
const double         prop_noeuds_actifs
)  
{

  if (fmt_tab == NULL) {
    _load_intern_fmt();
  }

  /* Look for fmt */

  int fmt_id = -1;
  for (int i = 0; i < n_fmt_tab; i++) {
    if (!strcmp(fmt, fmt_tab[i]->name)) {
      fmt_id = i;
      break;
    }
  }

  if (fmt_id == -1) {
    fprintf(stderr, "Error PDM_writer_create : unknown format '%s'", fmt);
    abort();
  }
  
  /* Mise a jour du tableau de stockage */

  PDM_io_mkdir(rep_sortie);

  if (cs_tab == NULL) {
    l_cs_tab = 4;
    cs_tab = (PDM_writer_t **) malloc(l_cs_tab * sizeof(PDM_writer_t *));
    for (int i = 0; i < l_cs_tab; i++) 
      cs_tab[i] = NULL;
  } 

  if (l_cs_tab <= n_cs_tab) {
    int p_l_cs_tab = l_cs_tab;
    l_cs_tab = 2 * l_cs_tab;
    cs_tab = (PDM_writer_t**) realloc((void*) cs_tab, l_cs_tab * sizeof(PDM_writer_t *));

    for (int i = p_l_cs_tab; i < l_cs_tab; i++) 
      cs_tab[i] = NULL;
  }

  /* Creation du r�pertoire de sortie si non cr�� */

  mkdir(rep_sortie, 0775); 

  /* Recherche de la premiere place libre pour stocker le bloc */

  int id_cs = 0;
  while (cs_tab[id_cs] != NULL) 
    id_cs++;

  /* Allocation de la structure PDM_writer_t */

  PDM_writer_t *cs = (PDM_writer_t *) malloc(sizeof(PDM_writer_t));

  n_cs_tab += 1;

  /* Initialisation de la structure PDM_writer_t */

  cs->fmt_id      = fmt_id;     /* Format de sortie */
  cs->fmt_fic     = fmt_fic;    /* Format du fichier de sortie */
  cs->topologie   = topologie;  /* Type de toplogie du maillage */
  cs->st_reprise  = st_reprise; /* Reprise d'une sortie existante */

  size_t l_rep_sortie = strlen(rep_sortie);
  cs->rep_sortie = (char *) malloc(sizeof(char) * (l_rep_sortie + 1));
  strcpy(cs->rep_sortie, rep_sortie);  /* Nom du repertoire de sortie */

  size_t l_nom_sortie = strlen(nom_sortie);
  cs->nom_sortie = (char *) malloc(sizeof(char) * (l_nom_sortie + 1));
  strcpy(cs->nom_sortie, nom_sortie);  /* Nom de la sortie */

  cs->pdm_mpi_comm    = pdm_mpi_comm;  /* Communicateur MPI */
  cs->sortie_fmt  = NULL;      /* Pointeur sur l'objet sortie de format fmt */    

  cs->var_tab     = NULL;      /* Tableau des variables */
  cs->l_var_tab   = 0;         /* Taille du tableau des variables */
  cs->n_var_tab   = 0;         /* Nombre de variables dans le tableau des variables */
  cs->geom_tab    = NULL;      /* Tableau des geometries */
  cs->l_geom_tab  = 0;         /* Taille du tableau des geometries */
  cs->n_geom_tab  = 0;         /* Nombre de geometries dans le tableau des geoemtries */
  cs->physical_time = 0;       /* Temps physique de simulation */
  cs->acces       = acces;
  cs->prop_noeuds_actifs = prop_noeuds_actifs;

  cs_tab[id_cs] = cs;

  /* Appel de la fonction complementaire propre au format */
  
  if (fmt_tab[cs->fmt_id]->create_fct != NULL) {
    (fmt_tab[cs->fmt_id]->create_fct) (cs);
  }

  return id_cs;

}

/*----------------------------------------------------------------------------
 * Libere un objet CS (Cedre Sortie) et retourne un pointeur NULL si pas d'erreur
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *
 *
 *----------------------------------------------------------------------------*/

void 
PROCF (pdm_writer_free, PDM_WRITER_FREE)
(
int        *id_cs
)
{
  PDM_writer_free(*id_cs);
}

void  
PDM_writer_free
(
const int   id_cs
)
{

  /* Recherche de l'objet cs courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  /* Appel de la fonction complementaire propre au format */

  if (fmt_tab[cs->fmt_id]->free_fct != NULL) {
    (fmt_tab[cs->fmt_id]->free_fct) (cs);
  }

  /* Liberation des diff�rents �l�m�nts de la structure */

  free(cs->rep_sortie);
  free(cs->nom_sortie);

  /* Liberation des variables */

  for (int i = 0; i < cs->l_var_tab; i++) {
    if (cs->var_tab[i] != NULL) {
      PDM_writer_var_free(id_cs, i);
    }
  }

  if (cs->var_tab != NULL) {
    free(cs->var_tab);
    cs->var_tab = NULL;
  }
  cs->l_var_tab = 0;


  /* Liberation de la g�om�trie */

  for (int i = 0; i < cs->l_geom_tab; i++) {
    if (cs->geom_tab[i] != NULL) {
      PDM_writer_geom_free(id_cs, i);
    }
  }

  if (cs->geom_tab != NULL) {
    free(cs->geom_tab);
    cs->geom_tab = NULL;
  }
  cs->l_geom_tab = 0;

  /* Liberation de la structure */

  free(cs);

  cs_tab[id_cs] = NULL;
  n_cs_tab -= 1;

  if (n_cs_tab == 0) {
    free(cs_tab);
    cs_tab = NULL;
    if (n_intern_fmt == n_fmt_tab) {
      for (int i = 0; i < n_fmt_tab; i++) {
        free (fmt_tab[i]);
      }
      free (fmt_tab);
      fmt_tab = NULL;
    }
  }
  
}


/*----------------------------------------------------------------------------
 * Debut d'increment
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   delta_t         <-- Delta de temps par rapport au dernier increment
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_step_beg, PDM_WRITER_STEP_BEG)
(
int           *id_cs,
double        *physical_time
)
{
  PDM_writer_step_beg(*id_cs,
              *physical_time);
}

void
PDM_writer_step_beg
(
const int      id_cs,
const double   physical_time
)
{
  /* Recherche de l'objet cs courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);
  cs->physical_time = physical_time;

  /* Appel de la fonction complementaire propre au format */
  
  if (fmt_tab[cs->fmt_id]->beg_step_fct != NULL) {
    (fmt_tab[cs->fmt_id]->beg_step_fct) (cs);
  }

}

/*----------------------------------------------------------------------------
 * Fin d'increment
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_step_end, PDM_WRITER_STEP_END)
(
int          *id_cs
)
{
  PDM_writer_step_end(*id_cs);
}

void
PDM_writer_step_end
(
const int     id_cs
)
{
  /* Recherche de l'objet cs courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  /* Appel de la fonction complementaire propre au format */
  
  if (fmt_tab[cs->fmt_id]->end_step_fct != NULL) {
    (fmt_tab[cs->fmt_id]->end_step_fct) (cs);
  }

}

/*----------------------------------------------------------------------------
 * Cree une nouvelle geometrie dans l'objet CS (Cedre Sortie)
 *
 * parameters :
 *   id_cs            <-- Identificateur de l'objet cs
 *   nom_geom         <-- Nom de l'objet geometrique
 *   st_decoup_poly2d <-- Active le decoupage des polygones 
 *   st_decoup_poly3d <-- Active le decoupage des polyedres
 *
 * return :
 *                   --> Identificateur de l'objet geom dans cs 
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_create_cf, PDM_WRITER_GEOM_CREATE_CF)
(
int           *id_cs,
char          *nom_geom,
int           *st_decoup_poly2d,
int           *st_decoup_poly3d,
int           *l_nom_geom,
int           *n_part,
int           *id_geom
ARGF_SUPP_CHAINE
)
{
  char *nom_geom_c = PDM_fortran_to_c_string(nom_geom, *l_nom_geom);
  
  *id_geom = PDM_writer_geom_create(*id_cs,
                          nom_geom_c,
                          (PDM_writer_statut_t) *st_decoup_poly2d,
                          (PDM_writer_statut_t) *st_decoup_poly3d,
                          *n_part);

  free(nom_geom_c);
}

int 
PDM_writer_geom_create
(
const int               id_cs,
const char             *nom_geom,
const PDM_writer_statut_t       st_decoup_poly2d,
const PDM_writer_statut_t       st_decoup_poly3d,
const int               n_part
)
{
  /* Erreur si le d�coupage des polygones ou polyedres est choisi */
  
  if (n_part <= 0) {
    fprintf(stderr, "Erreur cs_geom_create : Le nombre de partition doit etre >\n"
                    "                      Ajuster le communicateur MPI ou\n"
                    "                      Creer un sous-domaine avec 0 element\n");
  }

  if ((st_decoup_poly2d == 1) || (st_decoup_poly3d == 1)) {
    fprintf(stderr, "Erreur cs_geom_create : Les fonctions de decoupage ne sont pas operationnelles\n");
    abort();
  }

  /* Recherche de l'objet cs courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  /* Mise a jour du tableau de stockage */

  if (cs->geom_tab == NULL) {
    cs->l_geom_tab = 4;
    cs->geom_tab = (PDM_writer_geom_t **) malloc(cs->l_geom_tab * sizeof(PDM_writer_geom_t *));
    for (int i = 0; i < cs->l_geom_tab; i++) 
      cs->geom_tab[i] = NULL;
  } 
  
  if (cs->l_geom_tab <= cs->n_geom_tab) {
    int p_l_geom_tab = cs->l_geom_tab;
    cs->l_geom_tab = 2 * cs->l_geom_tab;
    cs->geom_tab = (PDM_writer_geom_t**) realloc((void*) cs->geom_tab, cs->l_geom_tab * sizeof(PDM_writer_geom_t *));
    
    for (int i = p_l_geom_tab; i < cs->l_geom_tab; i++) 
      cs->geom_tab[i] = NULL;
  }

  /* Recherche de la premiere place libre pour stocker le bloc */

  int id_geom = 0;
  while (cs->geom_tab[id_geom] != NULL) 
    id_geom++;

  /* Allocation de la structure PDM_writer_geom_t */

  PDM_writer_geom_t *geom = (PDM_writer_geom_t *) malloc(sizeof(PDM_writer_geom_t));

  cs->n_geom_tab += 1;

  cs->geom_tab[id_geom] = geom;

  /* Initialisation de la structure PDM_writer_geom_t */

  _geom_init(geom);

  geom->_cs = cs;
  geom->pdm_mpi_comm = cs->pdm_mpi_comm;
  size_t l_nom_geom = strlen(nom_geom);
  geom->nom_geom = (char *) malloc(sizeof(char) * (l_nom_geom + 1));
  strcpy(geom->nom_geom, nom_geom);  /* Nom de la geometrie */

  geom->n_part = n_part;

  geom->som = (PDM_writer_som_t  **) malloc(n_part * sizeof(PDM_writer_som_t *));     /* Nombre de sommets */
  geom->n_cell = (PDM_l_num_t  *) malloc(n_part * sizeof(PDM_l_num_t));     /* Nombre de sommets */

  for (int i = 0; i < n_part; i++) {
    geom->n_cell[i] = 0;
    geom->som[i] = (PDM_writer_som_t *) malloc (sizeof(PDM_writer_som_t));
    PDM_writer_som_t *_som = geom->som[i];
    _som->parent        = NULL;
    _som->n_som         = 0;
    _som->coords        = NULL;
    _som->_coords       = NULL;
    _som->_numabs        = NULL;
    _som->_numparent = NULL;
  }

  geom->prepa_blocs = NULL;
  geom->num_cell_parent_to_local = NULL;

  /* Appel de la fonction complementaire propre au format */
    
  if (fmt_tab[cs->fmt_id]->geom_create_fct != NULL) {
    (fmt_tab[cs->fmt_id]->geom_create_fct) (geom);
  }

  return id_geom;
}

/*----------------------------------------------------------------------------
 * Definition des coordonnees de la partition courante          
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part          <-- Indice de partition
 *   n_som           <-- Nombre de sommets de la partition
 *   coords          <-- Coordonnes des sommets            
 *   numabs          <-- Numerotation absolue des sommets     
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_coord_set, PDM_WRITER_GEOM_COORD_SET)
(
int             *id_cs,
int             *id_geom,  
int             *id_part, 
int             *n_som,  
PDM_real_t       *coords,  
PDM_g_num_t       *numabs
)
{
  PDM_writer_geom_coord_set(*id_cs,
                    *id_geom,  
                    *id_part, 
                    *n_som,  
                    coords,  
                    numabs);
}

void
PDM_writer_geom_coord_set
(
const int        id_cs,
const int        id_geom,  
const int        id_part, 
const int        n_som,  
const PDM_real_t *coords,  
const PDM_g_num_t *numabs
)
{
  

  /* Acces aux sommets de la partition */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);

  if (geom->n_part == 0) {
    fprintf(stderr, "Erreur PDM_writer_geom_coord_set : Le nombre de partitions n'a pas ete defini\n");
    abort();
  } 
  
  PDM_writer_som_t *som = geom->som[id_part];

  if ((som->_coords != NULL) ||
      (som->_numabs != NULL)) {
    fprintf(stderr, "Erreur PDM_writer_geom_coord_set : Les sommets de la partition ont deja ete definis\n");
    abort();
  }

  /* Mapping memoire */

  som->n_som   = n_som;
  som->_coords = coords;
  som->_numabs = numabs;
}




/*----------------------------------------------------------------------------
 * Definition des coordonnees des sommets de la partition courante a partir
 *          
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part         <-- Indice de partition
 *   n_som           <-- Nombre de sommets de la partition
 *   n_som_parent    <-- Nombre de sommets parent
 *   numabs          <-- Numerotation absolue des sommets (size = n_som)    
 *   num_parent      <-- Numerotation des sommets dans la numerotation parente (size = n_som)    
 *   coords_parent   <-- Coordonnes des sommets parents (size = 3 * n_som_parent)            
 *   numabs_parent   <-- Numerotation absolue des sommets parents (size = n_som_parent)    
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_coord_from_parent_set, PDM_WRITER_GEOM_COORD_FROM_PARENT_SET)
(
int             *id_cs,
int             *id_geom,  
int             *id_part, 
int             *n_som,  
int             *n_som_parent,  
PDM_g_num_t     *numabs,
int             *num_parent,
PDM_real_t      *coords_parent,  
PDM_g_num_t     *numabs_parent
)
{
  PDM_writer_geom_coord_from_parent_set (*id_cs,
                                 *id_geom,  
                                 *id_part, 
                                 *n_som,  
                                 *n_som_parent,  
                                 numabs,
                                 num_parent,
                                 coords_parent,
                                 numabs_parent);
}

void
PDM_writer_geom_coord_from_parent_set
(
const int        id_cs,
const int        id_geom,  
const int        id_part, 
const int        n_som,  
const int        n_som_parent,  
const PDM_g_num_t *numabs,
const int       *num_parent,
const PDM_real_t *coords_parent,  
const PDM_g_num_t *numabs_parent
)
{

  /* Acces aux sommets de la partition */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);

  if (geom->n_part == 0) {
    fprintf(stderr, "Erreur PDM_writer_geom_coord_from_parent_set : Le nombre de partitions n'a pas ete defini\n");
    abort();
  } 
  
  PDM_writer_som_t *som = geom->som[id_part];
  
  if ((som->_coords != NULL) ||
      (som->_numabs != NULL)) {
    fprintf(stderr, "Erreur PDM_writer_geom_coord_from_parent_set : Les sommets de la partition ont deja ete definis\n");
    abort();
  }

  /* Mapping memoire et allocation */

  som->parent     = (PDM_writer_som_t *) malloc (sizeof (PDM_writer_som_t));
  PDM_writer_som_t *_parent = som->parent;
  _parent->parent = NULL;
  _parent->n_som = n_som_parent;
  _parent->coords = NULL;
  _parent->_coords = coords_parent ;
  _parent->_numabs = numabs_parent;
  _parent->_numparent = NULL;
  
  som->n_som      = n_som;
  som->coords     = malloc (sizeof(double) * 3 * n_som);
  som->_coords    = som->coords;
  som->_numabs    = numabs;
  som->_numparent = num_parent;

  for (int i = 0; i < n_som; i++) {
    int i_parent = num_parent[i] - 1;
    for (int j = 0; j < 3; j++) {
      som->coords[3*i+j] = _parent->_coords[3*i_parent+j];
    }
  }

}

/*----------------------------------------------------------------------------
 * Ajout d'un bloc d'elements d'un type donne
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   t_elt           <-- Type d'element
 *
 * return :
 *                   --> Identificateur du bloc
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_add, PDM_WRITER_GEOM_BLOC_ADD)
(
int   *id_cs,
int   *id_geom,  
PDM_writer_statut_t  *st_free_data,  
int   *t_elt,
int   *id_bloc  
) 
{
  *id_bloc = PDM_writer_geom_bloc_add(*id_cs,
                              *id_geom,
                              *st_free_data,
                              (PDM_writer_elt_geom_t) *t_elt);
}

int 
PDM_writer_geom_bloc_add 
(
const int            id_cs,
const int            id_geom,   
PDM_writer_statut_t          st_free_data,  
const PDM_writer_elt_geom_t  t_elt
) 
{
  /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);

  /* Creation du bloc */

  int id_bloc = 0;

  switch (t_elt) {

  case PDM_writer_POINT    :    
  case PDM_writer_BAR2     :   
  case PDM_writer_TRIA3    :    
  case PDM_writer_QUAD4    :    
  case PDM_writer_TETRA4   :     
  case PDM_writer_PYRAMID5 :     
  case PDM_writer_PRISM6   :     
  case PDM_writer_HEXA8    : 
    {
      /* Mise a jour du tableau de stockage */

      if (geom->blocs_std == NULL) {
        geom->l_blocs_std = 4;
        geom->blocs_std = (PDM_writer_bloc_std_t **) malloc(geom->l_blocs_std * sizeof(PDM_writer_bloc_std_t *));
        for (int i = 0; i < geom->l_blocs_std; i++) 
          geom->blocs_std[i] = NULL;
      } 
      
      if (geom->l_blocs_std <= geom->n_blocs_std) {
        int p_l_blocs_std = geom->l_blocs_std;
        geom->l_blocs_std = 2 * geom->l_blocs_std;
        geom->blocs_std = (PDM_writer_bloc_std_t **) realloc((void*) geom->blocs_std, geom->l_blocs_std *
                                                     sizeof(PDM_writer_bloc_std_t*));
        
        for (int i = p_l_blocs_std; i < geom->l_blocs_std; i++) 
          geom->blocs_std[i] = NULL;
      }
      
      /* Recherche de la premiere place libre pour stocker le bloc */
      
      
      while (geom->blocs_std[id_bloc] != NULL) 
        id_bloc++;
      
      /* Allocation du bloc */
      
      PDM_writer_bloc_std_t *bloc_std = (PDM_writer_bloc_std_t *) malloc(sizeof(PDM_writer_bloc_std_t));
      geom->n_blocs_std += 1;
      
      geom->blocs_std[id_bloc] = bloc_std;

      /* Intialisation du bloc */

      bloc_std->t_elt = t_elt;
      bloc_std->st_free_data = st_free_data;
      bloc_std->n_part = geom->n_part; 

      bloc_std->n_elt      = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t ) * bloc_std->n_part);
      bloc_std->_connec    = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_std->n_part);
      bloc_std->_num_part  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_std->n_part);
      bloc_std->_numabs    = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * bloc_std->n_part);
      bloc_std->numabs_int = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * bloc_std->n_part);

      for (int i = 0; i < bloc_std->n_part; i++) {
        bloc_std->n_elt[i]     = 0;
        bloc_std->_connec[i]   = NULL;
        bloc_std->_num_part[i] = NULL;
        bloc_std->_numabs[i]   = NULL;
        bloc_std->numabs_int[i]= NULL;
      }
      
      id_bloc += PDM_writer_DEB_ID_BLOC_STD;
      if (id_bloc >= PDM_writer_DEB_ID_BLOC_POLY2D) {
        fprintf(stderr, "Erreur PDM_writer_geom_bloc_add : Le nombre de blocs d'elements standard doit etre inferieur a %d\n", 
               PDM_writer_DEB_ID_BLOC_POLY2D);
        abort();
      }
    }
    
    break;

  case PDM_writer_POLY_2D  :    
    {
      /* Mise a jour du tableau de stockage */

      if (geom->blocs_poly2d == NULL) {
        geom->l_blocs_poly2d = 4;
        geom->blocs_poly2d = (PDM_writer_bloc_poly2d_t **) malloc(geom->l_blocs_poly2d * 
                                                          sizeof(PDM_writer_bloc_poly2d_t *));
        for (int i = 0; i < geom->l_blocs_poly2d; i++) 
          geom->blocs_poly2d[i] = NULL;
      } 
      
      if (geom->l_blocs_poly2d <= geom->n_blocs_poly2d) {
        int p_l_blocs_poly2d = geom->l_blocs_poly2d;
        geom->l_blocs_poly2d = 2 * geom->l_blocs_poly2d;
        geom->blocs_poly2d = (PDM_writer_bloc_poly2d_t **) realloc((void*) geom->blocs_poly2d,
                                                           geom->l_blocs_poly2d * 
                                                           sizeof(PDM_writer_bloc_poly2d_t*));
        
        for (int i = p_l_blocs_poly2d; i < geom->l_blocs_poly2d; i++) 
          geom->blocs_poly2d[i] = NULL;
      }
      
      /* Recherche de la premiere place libre pour stocker le bloc */
      
      while (geom->blocs_poly2d[id_bloc] != NULL) 
        id_bloc++;
      
      /* Allocation du bloc */
      
      PDM_writer_bloc_poly2d_t *bloc_poly2d = (PDM_writer_bloc_poly2d_t *) malloc(sizeof(PDM_writer_bloc_poly2d_t));
      geom->n_blocs_poly2d += 1;
      
      geom->blocs_poly2d[id_bloc] = bloc_poly2d;

      /* Intialisation du bloc */

      bloc_poly2d->st_decoup_poly2d  = geom->st_decoup_poly2d;
      bloc_poly2d->st_free_data = st_free_data;
      bloc_poly2d->n_part            = geom->n_part; 

      bloc_poly2d->n_elt       = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t ) * bloc_poly2d->n_part);
      bloc_poly2d->_connec_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_poly2d->n_part);
      bloc_poly2d->_connec     = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_poly2d->n_part);
      bloc_poly2d->_num_part   = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_poly2d->n_part);
      bloc_poly2d->_numabs     = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * bloc_poly2d->n_part);
      bloc_poly2d->numabs_int = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * bloc_poly2d->n_part);

      for (int i = 0; i < bloc_poly2d->n_part; i++) {
        bloc_poly2d->n_elt[i]     = 0;
        bloc_poly2d->_connec_idx[i] = NULL;
        bloc_poly2d->_connec[i]     = NULL;
        bloc_poly2d->_num_part[i]= NULL;
        bloc_poly2d->_numabs[i]     = NULL;
        bloc_poly2d->numabs_int[i] = NULL;
      }

      if (bloc_poly2d->st_decoup_poly2d == PDM_writer_ON) {
        bloc_poly2d->tri_idx  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_poly2d->n_part);
        bloc_poly2d->bloc_tri = NULL;
        bloc_poly2d->quad_idx  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_poly2d->n_part);
        bloc_poly2d->bloc_quad = NULL;
        bloc_poly2d->som_sup  = NULL;
        for (int i = 0; i < bloc_poly2d->n_part; i++) {
          bloc_poly2d->tri_idx[i]   = NULL;
          bloc_poly2d->quad_idx[i]   = NULL;
        }
      }
      else {
        bloc_poly2d->tri_idx   = NULL;
        bloc_poly2d->bloc_tri  = NULL;
        bloc_poly2d->quad_idx  = NULL;
        bloc_poly2d->bloc_quad = NULL;
        bloc_poly2d->som_sup   = NULL;
      }

      id_bloc += PDM_writer_DEB_ID_BLOC_POLY2D;
      if (id_bloc >= PDM_writer_DEB_ID_BLOC_POLY3D) {
        fprintf(stderr, "Erreur PDM_writer_geom_bloc_add :"
               " Le nombre de blocs d'elements poly2d doit etre inferieur a %d\n",
               PDM_writer_DEB_ID_BLOC_POLY3D - PDM_writer_DEB_ID_BLOC_POLY2D);
        abort();
      }
    }

    break;

  case PDM_writer_POLY_3D  :     
    {
      /* Mise a jour du tableau de stockage */

      if (geom->blocs_poly3d == NULL) {
        geom->l_blocs_poly3d = 4;
        geom->blocs_poly3d = (PDM_writer_bloc_poly3d_t **) malloc(geom->l_blocs_poly3d * 
                                                          sizeof(PDM_writer_bloc_poly3d_t *));
        for (int i = 0; i < geom->l_blocs_poly3d; i++) 
          geom->blocs_poly3d[i] = NULL;
      } 
      
      if (geom->l_blocs_poly3d <= geom->n_blocs_poly3d) {
        int p_l_blocs_poly3d = geom->l_blocs_poly3d;
        geom->l_blocs_poly3d = 2 * geom->l_blocs_poly3d;
        geom->blocs_poly3d = (PDM_writer_bloc_poly3d_t **) realloc((void*) geom->blocs_poly3d, 
                                                           geom->l_blocs_poly3d * 
                                                           sizeof(PDM_writer_bloc_poly3d_t*));
        
        for (int i = p_l_blocs_poly3d; i < geom->l_blocs_poly3d; i++) 
          geom->blocs_poly3d[i] = NULL;
      }
      
      /* Recherche de la premiere place libre pour stocker le bloc */
      
      while (geom->blocs_poly3d[id_bloc] != NULL) 
        id_bloc++;
      
      /* Allocation du bloc */
      
      PDM_writer_bloc_poly3d_t *bloc_poly3d = (PDM_writer_bloc_poly3d_t *) malloc(sizeof(PDM_writer_bloc_poly3d_t));
      geom->n_blocs_poly3d += 1;
      
      geom->blocs_poly3d[id_bloc] = bloc_poly3d;

      /* Intialisation du bloc */

      bloc_poly3d->n_part            = geom->n_part; 
      bloc_poly3d->st_decoup_poly3d  = geom->st_decoup_poly3d;
      bloc_poly3d->st_free_data = st_free_data;

      bloc_poly3d->n_elt        = (PDM_l_num_t *)   malloc(sizeof(PDM_l_num_t ) * bloc_poly3d->n_part);
      bloc_poly3d->n_face       = (PDM_l_num_t *)   malloc(sizeof(PDM_l_num_t ) * bloc_poly3d->n_part);
      bloc_poly3d->_facsom_idx  = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * bloc_poly3d->n_part);
      bloc_poly3d->_facsom      = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * bloc_poly3d->n_part);
      bloc_poly3d->_cellfac_idx = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * bloc_poly3d->n_part);
      bloc_poly3d->_cellfac     = (PDM_l_num_t **)  malloc(sizeof(PDM_l_num_t *) * bloc_poly3d->n_part);
      bloc_poly3d->_numabs      = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * bloc_poly3d->n_part);
      bloc_poly3d->numabs_int   = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * bloc_poly3d->n_part);

      for (int i = 0; i < bloc_poly3d->n_part; i++) {
        bloc_poly3d->n_elt[i]        = 0;
        bloc_poly3d->n_face[i]       = 0;
        bloc_poly3d->_facsom_idx[i]  = NULL;
        bloc_poly3d->_facsom[i]      = NULL;
        bloc_poly3d->_cellfac_idx[i] = NULL;
        bloc_poly3d->_cellfac[i]     = NULL;
        bloc_poly3d->_numabs[i]      = NULL;
        bloc_poly3d->numabs_int[i]   = NULL;
      }

      if (bloc_poly3d->st_decoup_poly3d == PDM_writer_ON) {
        bloc_poly3d->tetra_idx  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_poly3d->n_part);
        bloc_poly3d->bloc_tetra = NULL;
        bloc_poly3d->pyra_idx  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * bloc_poly3d->n_part);
        bloc_poly3d->bloc_pyra = NULL;
        bloc_poly3d->som_sup   = NULL;
        for (int i = 0; i < bloc_poly3d->n_part; i++) {
          bloc_poly3d->tetra_idx[i]   = NULL;
          bloc_poly3d->pyra_idx[i]   = NULL;
        }
      }
      else {
        bloc_poly3d->tetra_idx  = NULL;
        bloc_poly3d->bloc_tetra = NULL;
        bloc_poly3d->pyra_idx  = NULL;
        bloc_poly3d->bloc_pyra = NULL;
        bloc_poly3d->som_sup    = NULL;
      }
      id_bloc += PDM_writer_DEB_ID_BLOC_POLY3D;
    }

    break;

  default :
    fprintf(stderr, "Erreur PDM_writer_geom_bloc_add : Type d'element inconnu\n");
    abort();

  }

  return id_bloc ;

}
 

/*----------------------------------------------------------------------------
 * Definition d'un bloc standard d'elements 
 *
 *  - PDM_writer_POINT :
 *
 *   1 x            
 *
 *  - PDM_writer_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_writer_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_writer_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_writer_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_writer_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_writer_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_writer_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_bloc         <-- Identificateur du bloc
 *   id_part         <-- Indice de partition
 *   n_elt           <-- Nombre d'elements dans la partition 
 *   connec          <-- Table de connectivite des elements
 *   numabs          <-- Numerotation absolue des elements
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_std_set, PDM_WRITER_GEOM_BLOC_STD_SET)
(
int           *id_cs,
int           *id_geom,  
int           *id_bloc,  
int           *id_part, 
int           *n_elt,    
PDM_l_num_t      *connec,   
PDM_g_num_t     *numabs
) 
{
  PDM_writer_geom_bloc_std_set (*id_cs,
                        *id_geom,
                        *id_bloc,
                        *id_part,
                        *n_elt,    
                        connec,   
                        numabs); 
}

void
PDM_writer_geom_bloc_std_set 
(
const int            id_cs,
const int            id_geom,  
const int            id_bloc,     
const int            id_part, 
const int            n_elt,    
      PDM_l_num_t      *connec,   
      PDM_g_num_t     *numabs
) 
{
  /* Acces a l'objet de geometrie courant */


  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);
  
  PDM_writer_bloc_std_t *bloc = _bloc_std_get(geom, id_bloc);
  
  if (id_part >= bloc->n_part) {
    fprintf(stderr, "Erreur PDM_writer_geom_bloc_std_set : Numero de partition trop grand\n");
    abort();
  }
 
  /* Mapping */
  
  //FIXME: Pourquoi geom->n_cell[id_part] += -bloc->n_elt[id_part]; et geom->n_cell[id_part] += n_elt
  geom->n_cell[id_part] += -bloc->n_elt[id_part]; 
  geom->n_cell[id_part] += n_elt;
  bloc->n_elt[id_part] = n_elt;
  bloc->_connec[id_part] = connec;
  bloc->_numabs[id_part] = numabs;
  
  for (int i = 0; i < n_elt; i++)
    geom->n_elt_abs = _lmax(geom->n_elt_abs, numabs[i]);
  
}


/*----------------------------------------------------------------------------
 * Ajout d'un bloc de polygones dans la partition courante
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part          <-- Indice de partition
 *   n_elt           <-- Nombre d'elements dans le bloc 
 *   connec_idx      <-- Index dans la table de connectivite (dim = n_elt+1)
 *   connec          <-- Table de connectivite des elements (dim = connec_idx[n_elt])
 *   numabs          <-- Numerotation absolue des elements
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_poly2d_set, PDM_WRITER_GEOM_BLOC_POLY2D_SET)
(
int           *id_cs,
int           *id_geom,  
int           *id_bloc, 
int           *id_part, 
PDM_l_num_t      *n_elt,    
PDM_l_num_t      *connec_idx,   
PDM_l_num_t      *connec,
PDM_g_num_t     *numabs
) 
{
  PDM_writer_geom_bloc_poly2d_set (*id_cs,
                           *id_geom,  
                           *id_bloc, 
                           *id_part, 
                           *n_elt,    
                           connec_idx,   
                           connec,
                           numabs); 
}
 
void
PDM_writer_geom_bloc_poly2d_set 
(
const int            id_cs,
const int            id_geom,  
const int            id_bloc, 
const int            id_part, 
const PDM_l_num_t       n_elt,    
      PDM_l_num_t      *connec_idx,   
      PDM_l_num_t      *connec,
      PDM_g_num_t     *numabs
) 
{
  /* Acces a l'objet de geometrie courant */


  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);
  
  PDM_writer_bloc_poly2d_t *bloc = _bloc_poly2d_get(geom, id_bloc);
  
  if (id_part >= bloc->n_part) {
    fprintf(stderr, "Erreur PDM_writer_geom_bloc_poly2d_set : Numero de partition trop grand\n");
    abort();
  }
 
  /* Mapping */

  geom->n_cell[id_part]      += -bloc->n_elt[id_part]; 
  geom->n_cell[id_part]      += n_elt;
  bloc->n_elt[id_part]       = n_elt;
  bloc->_connec_idx[id_part] = connec_idx;
  bloc->_connec[id_part]     = connec;
  bloc->_numabs[id_part]     = numabs;

  for (int i = 0; i < n_elt; i++)
    geom->n_elt_abs = _lmax(geom->n_elt_abs, numabs[i]);
}
 

/*----------------------------------------------------------------------------
 * Ajout d'un bloc de polyedres dans la partition courante
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_bloc         <-- Identificateur de bloc
 *   id_part         <-- Indice de partition
 *   n_elt           <-- Nombre d'elements dans le bloc 
 *   n_face          <-- Nombre de faces de chaque element (dim = n_elt)
 *   facsom_idx      <-- Index dans la table de connectivite des faces (dim = n_face_total+1)
 *   facsom          <-- Table de connectivite des faces (dim = facsom_idx[n_face_total}
 *   cellfac_idx     <-- Index dans la table de connectivite des cellules (dim = n_elt+1)
 *   cellfac         <-- Table de connectivite des elements (dim = cellfac_idx[n_elt])
 *   numabs          <-- Numerotation absolue des elements
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_bloc_poly3d_set, PDM_WRITER_GEOM_BLOC_POLY3D_SET)
(
int           *id_cs,
int           *id_geom,  
int           *id_bloc, 
int           *id_part, 
PDM_l_num_t      *n_elt,    
PDM_l_num_t      *n_face,   
PDM_l_num_t      *facsom_idx,   
PDM_l_num_t      *facsom,
PDM_l_num_t      *cellfac_idx,   
PDM_l_num_t      *cellfac,
PDM_g_num_t     *numabs
) 
{
  PDM_writer_geom_bloc_poly3d_set (*id_cs,
                           *id_geom,  
                           *id_bloc, 
                           *id_part, 
                           *n_elt,    
                           *n_face,   
                           facsom_idx,   
                           facsom,
                           cellfac_idx,   
                           cellfac,
                           numabs); 
}

void
PDM_writer_geom_bloc_poly3d_set 
(
const int        id_cs,
const int        id_geom,  
const int        id_bloc, 
const int        id_part, 
const PDM_l_num_t   n_elt,    
const PDM_l_num_t   n_face,   
      PDM_l_num_t  *facsom_idx,   
      PDM_l_num_t  *facsom,
      PDM_l_num_t  *cellfac_idx,   
      PDM_l_num_t  *cellfac,
      PDM_g_num_t *numabs
) 
{
  /* Acces a l'objet de geometrie courant */


  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);
  
  PDM_writer_bloc_poly3d_t *bloc = _bloc_poly3d_get(geom, id_bloc);
  
  if (id_part >= bloc->n_part) {
    fprintf(stderr, "Erreur PDM_writer_geom_bloc_poly3d_set : Numero de partition trop grand\n");
    abort();
  }
 
  /* Mapping */

  geom->n_cell[id_part]       += -bloc->n_elt[id_part]; 
  geom->n_cell[id_part]       += n_elt;
  bloc->n_elt[id_part]        = n_elt;
  bloc->n_face[id_part]       = n_face;
  bloc->_facsom_idx[id_part]  = facsom_idx;
  bloc->_facsom[id_part]      = facsom;
  bloc->_cellfac_idx[id_part] = cellfac_idx;
  bloc->_cellfac[id_part]     = cellfac;
  bloc->_numabs[id_part]      = numabs;

  for (int i = 0; i < n_elt; i++)
    geom->n_elt_abs = _lmax(geom->n_elt_abs, numabs[i]);
}


/*----------------------------------------------------------------------------
 *
 * Ajout de cellules 3D decrites en fonctions des faces. Cette fonction
 * d�termine les types des �l�ments et cr�e des blocs regrouppant les �l�ments
 * de m�me type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   n_cell          <-- Nombre de cellules 3D ajout�es         
 *   n_face          <-- Nombre de faces d�crites               
 *   face_som_idx    <-- Index de connectivite faces -> sommets
 *   face_som        <-- Connectivite faces -> sommets                                       
 *   cell_face_idx   <-- Index de connectivite cellules -> faces  
 *   cell_face       <-- Connectivite cellules -> faces
 *   numabs          <-- Numerotatio absolue des cellules 
 *   ind_num         --> Indirection vers la nouvelle numerotation des cellules 
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_cell3d_cellface_add, PDM_WRITER_GEOM_CELL3D_CELLFACE_ADD)
(
int         *id_cs,
int         *id_geom,
int         *id_part,
int         *n_cell,
int         *n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
) 
{
  PDM_writer_geom_cell3d_cellface_add(*id_cs,
                              *id_geom,
                              *id_part, 
                              *n_cell,
                              *n_face,
                              face_som_idx,
                              face_som_nb,
                              face_som,
                              cell_face_idx,
                              cell_face_nb,
                              cell_face,
                              numabs);
} 

void
PDM_writer_geom_cell3d_cellface_add
(
const int    id_cs,
const int    id_geom,
const int    id_part, 
const int    n_cell,
const int    n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
) 
{
  /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);
  
  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);
  int n_part = 0;

  if (geom->num_cell_parent_to_local == NULL) {
    geom->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * geom->n_part); 
    for (int ipart = 0; ipart < geom->n_part; ipart++) { 
      geom->num_cell_parent_to_local[ipart] = NULL;
    }
  }

  geom->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    geom->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (geom->prepa_blocs == NULL) {
    geom->prepa_blocs = (PDM_writer_geom_prepa_blocs_t *) malloc(sizeof(PDM_writer_geom_prepa_blocs_t));
    geom->prepa_blocs->t_add = 1;
    geom->prepa_blocs->n_tria_proc = 0;    /* Nb de triangles par proc */
    geom->prepa_blocs->n_quad_proc = 0;    /* Nb de quads par proc */
    geom->prepa_blocs->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    geom->prepa_blocs->n_tetra_proc = 0;   /* Nb de tetra par proc */
    geom->prepa_blocs->n_hexa_proc = 0;    /* Nb d'hexa par proc */
    geom->prepa_blocs->n_prism_proc = 0;   /* Nb de prisme par proc */
    geom->prepa_blocs->n_pyramid_proc = 0; /* Nb de pyramide par proc */
    geom->prepa_blocs->n_poly3d_proc = 0;  /* Nb de poly3d par proc */
    geom->prepa_blocs->n_cell = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_tetra = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_hexa = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_prism = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_pyramid = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_poly3d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->face_som_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part); 
    geom->prepa_blocs->face_som_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->face_som = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->cell_face_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->cell_face_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->cell_face = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part);
    geom->prepa_blocs->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*geom->n_part);
    for (int i = 0; i < geom->n_part; i++) {
      geom->prepa_blocs->add_etat[i] = 0;
    }
  }

  if (geom->prepa_blocs->t_add != 1) {
    fprintf(stderr, "Erreur Cs_geom_cell3d_cellface_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  /* Determination du type de chaque element */

  PDM_l_num_t cell_som_tria[18]; /* 6 triangles max in _type_cell_3D */
  PDM_l_num_t cell_som_quad[24]; /* 6 quadrangles max in _type_cell_3D */
  PDM_l_num_t n_tetra   = 0;
  PDM_l_num_t n_hexa    = 0;
  PDM_l_num_t n_prism   = 0;
  PDM_l_num_t n_pyramid = 0;
  PDM_l_num_t n_poly3d  = 0;
  
  for (int i = 0; i < n_cell; i++) {
    PDM_writer_elt_geom_t cell_type = _type_cell_3D(cell_face_nb[i],
                                                    cell_face + cell_face_idx[i] - 1,
                                                    face_som_idx,
                                                    face_som_nb,
                                                    face_som,
                                                    cell_som_tria,
                                                    cell_som_quad);
    switch(cell_type) {
    case PDM_writer_TETRA4 :
      n_tetra += 1;
      break;
    case PDM_writer_PYRAMID5 :
      n_pyramid += 1;
      break;
    case PDM_writer_PRISM6 :
      n_prism += 1;
      break;
    case PDM_writer_HEXA8 :
      n_hexa += 1;
      break;
    case PDM_writer_POLY_3D :
      n_poly3d += 1;
      break;
    default :
      break;
    }
  }

  geom->prepa_blocs->n_tetra_proc          += n_tetra;
  geom->prepa_blocs->n_hexa_proc           += n_hexa;
  geom->prepa_blocs->n_prism_proc          += n_prism;
  geom->prepa_blocs->n_pyramid_proc        += n_pyramid;
  geom->prepa_blocs->n_poly3d_proc         += n_poly3d;
  geom->prepa_blocs->n_tetra[id_part]       = n_tetra;
  geom->prepa_blocs->n_hexa[id_part]        = n_hexa;
  geom->prepa_blocs->n_prism[id_part]       = n_prism;
  geom->prepa_blocs->n_pyramid[id_part]     = n_pyramid;
  geom->prepa_blocs->n_poly3d[id_part]      = n_poly3d;
  geom->prepa_blocs->face_som_idx[id_part]  = face_som_idx;
  geom->prepa_blocs->face_som_nb[id_part]   = face_som_nb;
  geom->prepa_blocs->face_som[id_part]      = face_som;
  geom->prepa_blocs->cell_face_idx[id_part] = cell_face_idx;
  geom->prepa_blocs->cell_face_nb[id_part]  = cell_face_nb;
  geom->prepa_blocs->cell_face[id_part]     = cell_face;
  geom->prepa_blocs->numabs[id_part]        = numabs;
  geom->prepa_blocs->add_etat[id_part]      = 1;
  geom->prepa_blocs->n_face[id_part]        = n_face;
  geom->prepa_blocs->n_cell[id_part]          = n_cell;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < geom->n_part; i++) {
    if (geom->prepa_blocs->add_etat[i] == 1)
      n_part += 1;
  }

  if (geom->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[5];
    PDM_l_num_t som_elts[5];

    elts[0] = geom->prepa_blocs->n_tetra_proc > 0;
    elts[1] = geom->prepa_blocs->n_hexa_proc  > 0;
    elts[2] = geom->prepa_blocs->n_prism_proc  > 0;
    elts[3] = geom->prepa_blocs->n_pyramid_proc  > 0;
    elts[4] = geom->prepa_blocs->n_poly3d_proc  > 0;
    
    PDM_MPI_Allreduce(elts, som_elts, 5, PDM_MPI_INT, PDM_MPI_SUM, cs->pdm_mpi_comm);

    int id_bloc_tetra4;
    int id_bloc_hexa8;
    int id_bloc_prism6;
    int id_bloc_pyramid5;
    int id_bloc_poly_3d;

    if (som_elts[0] > 0)
      id_bloc_tetra4 = PDM_writer_geom_bloc_add(id_cs,
                                        id_geom,
                                        PDM_writer_ON,
                                        PDM_writer_TETRA4);

    if (som_elts[1] > 0)
      id_bloc_hexa8 = PDM_writer_geom_bloc_add(id_cs,
                                       id_geom,
                                       PDM_writer_ON,
                                       PDM_writer_HEXA8);
    
    if (som_elts[2] > 0)
      id_bloc_prism6 = PDM_writer_geom_bloc_add(id_cs,
                                        id_geom,
                                        PDM_writer_ON,
                                        PDM_writer_PRISM6);

    if (som_elts[3] > 0)
      id_bloc_pyramid5 = PDM_writer_geom_bloc_add(id_cs,
                                          id_geom,
                                          PDM_writer_ON,
                                          PDM_writer_PYRAMID5);

    if (som_elts[4] > 0)
      id_bloc_poly_3d = PDM_writer_geom_bloc_add(id_cs,
                                         id_geom,
                                         PDM_writer_ON,
                                         PDM_writer_POLY_3D);
                                                   
    /* Determination de la connectivite de chaque element */
    

    for (int ipart = 0; ipart < geom->n_part; ipart++) {
      
      PDM_l_num_t n_cell_courant = geom->prepa_blocs->n_cell[ipart];
      PDM_l_num_t *num_cell_parent_to_local_courant = geom->num_cell_parent_to_local[ipart];
      PDM_l_num_t *face_som_idx_courant = geom->prepa_blocs->face_som_idx[ipart];
      PDM_l_num_t *face_som_nb_courant = geom->prepa_blocs->face_som_nb[ipart];
      PDM_l_num_t *face_som_courant = geom->prepa_blocs->face_som[ipart];
      PDM_l_num_t *cell_face_idx_courant = geom->prepa_blocs->cell_face_idx[ipart];
      PDM_l_num_t *cell_face_nb_courant = geom->prepa_blocs->cell_face_nb[ipart];
      PDM_l_num_t *cell_face_courant = geom->prepa_blocs->cell_face[ipart];
      PDM_g_num_t *numabs_courant = geom->prepa_blocs->numabs[ipart];
  
      PDM_l_num_t n_face_part   = geom->prepa_blocs->n_face[ipart];
 
      PDM_l_num_t n_tetra_part   = geom->prepa_blocs->n_tetra[ipart];
      PDM_l_num_t n_hexa_part    = geom->prepa_blocs->n_hexa[ipart];
      PDM_l_num_t n_prism_part   = geom->prepa_blocs->n_prism[ipart];
      PDM_l_num_t n_pyramid_part = geom->prepa_blocs->n_pyramid[ipart];
      PDM_l_num_t n_poly3d_part  = geom->prepa_blocs->n_poly3d[ipart];
      
      PDM_l_num_t *connec_tetra = NULL;
      PDM_l_num_t *connec_hexa = NULL;
      PDM_l_num_t *connec_prism = NULL;
      PDM_l_num_t *connec_pyramid = NULL;
 
      PDM_g_num_t *numabs_tetra = NULL;
      PDM_g_num_t *numabs_hexa = NULL;
      PDM_g_num_t *numabs_prism = NULL;
      PDM_g_num_t *numabs_pyramid = NULL;
      PDM_g_num_t *numabs_poly3d = NULL;
    
      if (n_tetra_part > 0) {
        connec_tetra = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 *n_tetra_part);
        numabs_tetra = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tetra_part);
      }

      if (n_hexa_part > 0) {
        connec_hexa = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 8 * n_hexa_part);
        numabs_hexa = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_hexa_part);
      }

      if (n_prism_part > 0) {
        connec_prism = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 6 * n_prism_part);
        numabs_prism = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_prism_part);
      }

      if (n_pyramid_part > 0) {
        connec_pyramid = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 5 * n_pyramid_part);
        numabs_pyramid = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_pyramid_part);
      }

      if (n_poly3d_part > 0) {
        numabs_poly3d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly3d_part);
      }

      PDM_l_num_t *connec_tetra_courant = connec_tetra;
      PDM_l_num_t *connec_hexa_courant = connec_hexa;
      PDM_l_num_t *connec_prism_courant = connec_prism;
      PDM_l_num_t *connec_pyramid_courant = connec_pyramid;

      PDM_g_num_t *numabs_tetra_courant = numabs_tetra;
      PDM_g_num_t *numabs_hexa_courant = numabs_hexa;
      PDM_g_num_t *numabs_prism_courant = numabs_prism;
      PDM_g_num_t *numabs_pyramid_courant = numabs_pyramid;
      PDM_g_num_t *numabs_poly3d_courant = numabs_poly3d;

      PDM_l_num_t *tag_face_poly3d = NULL;
      PDM_l_num_t  n_face_poly = 0;
      PDM_l_num_t *facsom_poly_idx = NULL;
      PDM_l_num_t *facsom_poly = NULL;
      PDM_l_num_t *cellfac_poly_idx = NULL;
      PDM_l_num_t *cellfac_poly = NULL;
      PDM_l_num_t l_cellfac_poly = 0;

      if (n_poly3d_part > 0) {
        tag_face_poly3d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_face_part);
        for (int i = 0; i < n_face_part; i++) {
          tag_face_poly3d[i] = -1;
        }
        cellfac_poly_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly3d_part + 1));
        cellfac_poly_idx[0] = 0;
      }

      PDM_l_num_t idx_tetra = 0;
      PDM_l_num_t idx_hexa = n_tetra_part;
      PDM_l_num_t idx_prism = idx_hexa + n_hexa_part;
      PDM_l_num_t idx_pyramid = idx_prism + n_prism_part;
      PDM_l_num_t idx_poly3d = idx_pyramid + n_pyramid_part;

      n_poly3d_part = 0; 
      for (int i = 0; i < n_cell_courant; i++) {
        num_cell_parent_to_local_courant[i] = 0;
        PDM_writer_elt_geom_t cell_type = _type_cell_3D(cell_face_nb_courant[i],
                                                cell_face_courant + cell_face_idx_courant[i] - 1,
                                                face_som_idx_courant,
                                                face_som_nb_courant,
                                                face_som_courant,
                                                cell_som_tria,
                                                cell_som_quad);
        
        switch(cell_type) {
        case PDM_writer_TETRA4 :
          _connec_tetra(geom->som[ipart],
                        cell_som_tria,
                        connec_tetra_courant);
          *numabs_tetra_courant = numabs_courant[i];
          numabs_tetra_courant += 1;
          connec_tetra_courant += 4;
          num_cell_parent_to_local_courant[i] = idx_tetra++;
          break;
        case PDM_writer_HEXA8 :
          _connec_hexa(geom->som[ipart],
                       cell_som_quad,
                       connec_hexa_courant);
          *numabs_hexa_courant = numabs_courant[i];
          numabs_hexa_courant += 1;
          connec_hexa_courant += 8;
          num_cell_parent_to_local_courant[i] = idx_hexa++;
          break;
        case PDM_writer_PRISM6 :
          _connec_prism(geom->som[ipart],
                        cell_som_tria,
                        cell_som_quad,
                        connec_prism_courant);
          *numabs_prism_courant = numabs_courant[i];
          numabs_prism_courant += 1;
          connec_prism_courant += 6;
          num_cell_parent_to_local_courant[i] = idx_prism++;
          break;
        case PDM_writer_PYRAMID5 :
          _connec_pyramid(geom->som[ipart],
                          cell_som_tria,
                          cell_som_quad,
                          connec_pyramid_courant);
          *numabs_pyramid_courant = numabs_courant[i];
          numabs_pyramid_courant += 1;
          connec_pyramid_courant += 5;
          num_cell_parent_to_local_courant[i] = idx_pyramid++;
          break;
        case PDM_writer_POLY_3D : 
          {
            PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - 1;
            for (int j = 0; j < cell_face_nb_courant[i]; j++) {
              tag_face_poly3d[cell_face_cell[j] - 1] = 0;
            }
            *numabs_poly3d_courant = numabs_courant[i];
            numabs_poly3d_courant += 1;
            l_cellfac_poly += cell_face_nb_courant[i];
            cellfac_poly_idx[n_poly3d_part+1] = l_cellfac_poly;
            n_poly3d_part += 1;
            num_cell_parent_to_local_courant[i] = idx_poly3d++;
            break;
          }
        default :
          break;
        }
      }
        
      if (n_poly3d_part > 0) {
        cellfac_poly = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_cellfac_poly);
        
        /* Stockage des faces du bloc */
        
        n_face_poly = 0;
        PDM_l_num_t l_facsom_poly = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] == 0) {
            tag_face_poly3d[i] = n_face_poly++;
            l_facsom_poly += face_som_nb_courant[i];
          }
        }
        
        facsom_poly_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_face_poly + 1));
        facsom_poly = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_facsom_poly);
        
        facsom_poly_idx[0] = 0;
        PDM_l_num_t idx_facsom_poly = 0;
        PDM_l_num_t idx_facsom = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] >= 0) {
            PDM_l_num_t ideb = face_som_idx_courant[i] - 1;
            PDM_l_num_t ifin = ideb + face_som_nb_courant[i];
            facsom_poly_idx[idx_facsom+1] = facsom_poly_idx[idx_facsom] + face_som_nb_courant[i];
            idx_facsom += 1;
            for (int j = ideb; j < ifin; j++) {
              facsom_poly[idx_facsom_poly++] = face_som_courant[j];
            }
          }
        }

        /* Remplissage de la structure cellfac_poly */

        l_cellfac_poly = 0;
        for (int i = 0; i < n_cell_courant; i++) {
          PDM_writer_elt_geom_t cell_type = _type_cell_3D(cell_face_nb_courant[i],
                                                  cell_face_courant + cell_face_idx_courant[i] - 1,
                                                  face_som_idx_courant,
                                                  face_som_nb_courant,
                                                  face_som_courant,
                                                  cell_som_tria,
                                                  cell_som_quad);
        
          switch(cell_type) {
            
          case PDM_writer_POLY_3D : 
            {
              PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - 1;
              for (int j = 0; j < cell_face_nb_courant[i]; j++) {
                cellfac_poly[l_cellfac_poly++] = tag_face_poly3d[cell_face_cell[j] - 1] + 1;
              }
              break;
            }
         default:
            break;
          }
        }
        free(tag_face_poly3d);
      }

      if (som_elts[0] > 0)
        PDM_writer_geom_bloc_std_set(id_cs,
                                     id_geom,
                                     id_bloc_tetra4,
                                     ipart,
                                     n_tetra_part,
                                     connec_tetra,
                                     numabs_tetra);

      if (som_elts[1] > 0)
        PDM_writer_geom_bloc_std_set(id_cs,
                             id_geom,
                             id_bloc_hexa8,
                             ipart,
                             n_hexa_part,
                             connec_hexa,
                             numabs_hexa);
    
      if (som_elts[2] > 0)
        PDM_writer_geom_bloc_std_set(id_cs,
                             id_geom,
                             id_bloc_prism6,
                             ipart,
                             n_prism_part,
                             connec_prism,
                             numabs_prism);

      if (som_elts[3] > 0)
        PDM_writer_geom_bloc_std_set(id_cs,
                             id_geom,
                             id_bloc_pyramid5,
                             ipart,
                             n_pyramid_part,
                             connec_pyramid,
                             numabs_pyramid);

      if (som_elts[4] > 0)
        PDM_writer_geom_bloc_poly3d_set(id_cs,
                                id_geom,
                                id_bloc_poly_3d,
                                ipart,
                                n_poly3d_part,
                                n_face_poly,
                                facsom_poly_idx,
                                facsom_poly,
                                cellfac_poly_idx,
                                cellfac_poly,
                                numabs_poly3d);
    }

    if (geom->prepa_blocs != NULL) {
      free(geom->prepa_blocs->n_cell);
      free(geom->prepa_blocs->n_face);
      free(geom->prepa_blocs->n_tetra);
      free(geom->prepa_blocs->n_hexa);
      free(geom->prepa_blocs->n_prism);
      free(geom->prepa_blocs->n_pyramid);
      free(geom->prepa_blocs->n_poly3d);
      free(geom->prepa_blocs->face_som_idx);
      free(geom->prepa_blocs->face_som_nb);
      free(geom->prepa_blocs->face_som);
      free(geom->prepa_blocs->cell_face_idx);
      free(geom->prepa_blocs->cell_face_nb);
      free(geom->prepa_blocs->cell_face);
      free(geom->prepa_blocs->add_etat);
      free(geom->prepa_blocs->numabs);
      free(geom->prepa_blocs);
      geom->prepa_blocs = NULL;
    }
  }
} 


/*----------------------------------------------------------------------------
 *
 * Ajout de cellules 2D decrites en fonctions des faces. Cette fonction
 * d�termine les types des �l�ments et cr�e des blocs regrouppant les �l�ments
 * de m�me type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   n_elt           <-- Nombre de cellules 3D ajout�es         
 *   n_face          <-- Nombre de faces d�crites               
 *   face_som_idx    <-- Index de connectivite faces -> sommets
 *   face_som        <-- Connectivite faces -> sommets                                       
 *   cell_face_idx   <-- Index de connectivite cellules -> faces  
 *   cell_face       <-- Connectivite cellules -> faces
 *   numabs          <-- Numerotatio absolue des cellules 
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_cell2d_cellface_add, PDM_WRITER_GEOM_CELL2D_CELLFACE_ADD)
(
int         *id_cs,
int         *id_geom,
int         *id_part,
int         *n_cell,
int         *n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
) 
{
  PDM_writer_geom_cell2d_cellface_add(*id_cs,
                              *id_geom,
                              *id_part, 
                              *n_cell,
                              *n_face,
                              face_som_idx,
                              face_som_nb,
                              face_som,
                              cell_face_idx,
                              cell_face_nb,
                              cell_face,
                              numabs);
} 

void
PDM_writer_geom_cell2d_cellface_add
(
const int          id_cs,
const int          id_geom,
const int          id_part, 
const int          n_cell,
const int          n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_l_num_t    *cell_face_idx,
PDM_l_num_t    *cell_face_nb,
PDM_l_num_t    *cell_face,
PDM_g_num_t   *numabs
) 
{
  /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);
  int n_part = 0;

  if (geom->num_cell_parent_to_local == NULL) {
    geom->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * geom->n_part); 
    for (int ipart = 0; ipart < geom->n_part; ipart++) { 
      geom->num_cell_parent_to_local[ipart] = NULL;
    }
  }

  geom->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    geom->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (geom->prepa_blocs == NULL) {
    geom->prepa_blocs = (PDM_writer_geom_prepa_blocs_t *) malloc(sizeof(PDM_writer_geom_prepa_blocs_t));
    geom->prepa_blocs->t_add = 2;
    geom->prepa_blocs->n_tria_proc = 0;    /* Nb de triangles par proc */
    geom->prepa_blocs->n_quad_proc = 0;    /* Nb de quads par proc */
    geom->prepa_blocs->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    geom->prepa_blocs->n_cell = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->l_connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->face_som_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part); 
    geom->prepa_blocs->face_som_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->face_som = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->cell_face_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->cell_face_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->cell_face = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part);
    geom->prepa_blocs->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*geom->n_part);
    for (int i = 0; i < geom->n_part; i++) {
      geom->prepa_blocs->add_etat[i] = 0;
    }
  }

  if (geom->prepa_blocs->t_add != 2) {
    fprintf(stderr, "Erreur Cs_geom_cell2d_cellface_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  PDM_l_num_t n_tria    = 0;
  PDM_l_num_t n_quad    = 0;
  PDM_l_num_t n_poly2d  = 0;
  PDM_l_num_t l_connec_poly2d = 0;

  for (int i = 0; i < n_cell; i++) {

    PDM_l_num_t n_face_cell = cell_face_nb[i];
    if (n_face_cell == 3)
      n_tria += 1;
    else if (n_face_cell == 4)
      n_quad += 1;
    else {
      n_poly2d  += 1;
      l_connec_poly2d += cell_face_nb[i]; 
    }
  }

  geom->prepa_blocs->n_tria_proc           += n_tria;
  geom->prepa_blocs->n_quad_proc           += n_quad;
  geom->prepa_blocs->n_poly2d_proc         += n_poly2d;
  geom->prepa_blocs->add_etat[id_part]      = 1;
  geom->prepa_blocs->n_cell[id_part]        = n_cell;
  geom->prepa_blocs->n_tria[id_part]        = n_tria; 
  geom->prepa_blocs->n_quad[id_part]        = n_quad; 
  geom->prepa_blocs->n_poly2d[id_part]      = n_poly2d;
  geom->prepa_blocs->l_connec_poly2d[id_part] = l_connec_poly2d;
  geom->prepa_blocs->face_som_idx[id_part]  = face_som_idx;
  geom->prepa_blocs->face_som_nb[id_part]   = face_som_nb;
  geom->prepa_blocs->face_som[id_part]      = face_som;
  geom->prepa_blocs->cell_face_idx[id_part] = cell_face_idx;
  geom->prepa_blocs->cell_face_nb[id_part]  = cell_face_nb;
  geom->prepa_blocs->cell_face[id_part]     = cell_face;
  geom->prepa_blocs->numabs[id_part]        = numabs;
  geom->prepa_blocs->add_etat[id_part]      = 1;
  geom->prepa_blocs->n_face[id_part]        = n_face;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < geom->n_part; i++) {
    if (geom->prepa_blocs->add_etat[i] == 1)
      n_part += 1;
  }

  if (geom->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[3];
    PDM_l_num_t som_elts[3];

    elts[0] = geom->prepa_blocs->n_tria_proc;
    elts[1] = geom->prepa_blocs->n_quad_proc;
    elts[2] = geom->prepa_blocs->n_poly2d_proc;

    PDM_MPI_Allreduce(elts, som_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, cs->pdm_mpi_comm);

    int id_bloc_tria3;
    int id_bloc_quad4;
    int id_bloc_poly_2d;

    if (som_elts[0] > 0)
      id_bloc_tria3 = PDM_writer_geom_bloc_add(id_cs,
                                       id_geom,
                                       PDM_writer_ON,
                                       PDM_writer_TRIA3);

    if (som_elts[1] > 0)
      id_bloc_quad4 = PDM_writer_geom_bloc_add(id_cs,
                                       id_geom,
                                       PDM_writer_ON,
                                       PDM_writer_QUAD4);
    
    if (som_elts[2] > 0)
      id_bloc_poly_2d = PDM_writer_geom_bloc_add(id_cs,
                                         id_geom,
                                         PDM_writer_ON,
                                         PDM_writer_POLY_2D);

    /* Determination de la connectivite de chaque element */

    for (int ipart = 0; ipart < geom->n_part; ipart++) {

      PDM_l_num_t n_cell_courant = geom->prepa_blocs->n_cell[ipart];
      PDM_l_num_t *num_cell_parent_to_local_courant = geom->num_cell_parent_to_local[ipart];
      PDM_l_num_t *face_som_courant = geom->prepa_blocs->face_som[ipart];
      PDM_l_num_t *cell_face_idx_courant = geom->prepa_blocs->cell_face_idx[ipart];
      PDM_l_num_t *cell_face_nb_courant = geom->prepa_blocs->cell_face_nb[ipart];
      PDM_l_num_t *cell_face_courant = geom->prepa_blocs->cell_face[ipart];
      PDM_g_num_t *numabs_courant = geom->prepa_blocs->numabs[ipart];
   
      n_tria   = geom->prepa_blocs->n_tria[ipart];
      n_quad    = geom->prepa_blocs->n_quad[ipart];
      n_poly2d  = geom->prepa_blocs->n_poly2d[ipart];
      l_connec_poly2d = geom->prepa_blocs->l_connec_poly2d[ipart];

      PDM_l_num_t *connec_tria = NULL;
      PDM_l_num_t *connec_quad = NULL;
      PDM_l_num_t *connec_poly2d = NULL;
      PDM_l_num_t *connec_poly2d_idx = NULL;
 
      PDM_g_num_t *numabs_tria = NULL;
      PDM_g_num_t *numabs_quad = NULL;
      PDM_g_num_t *numabs_poly2d = NULL;

      if (n_tria > 0) {
        connec_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 3 *n_tria);
        numabs_tria = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tria);
      }

      if (n_quad > 0) {
        connec_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 * n_quad);
        numabs_quad = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_quad);
      }

      if (n_poly2d > 0) {
        connec_poly2d_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly2d + 1));
        connec_poly2d_idx[0] = 0;
        connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_connec_poly2d);
        numabs_poly2d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly2d);
      }

      PDM_l_num_t *connec_tria_courant = connec_tria;
      PDM_l_num_t *connec_quad_courant = connec_quad;
      PDM_l_num_t *connec_poly2d_idx_courant = connec_poly2d_idx + 1;
      PDM_l_num_t *connec_poly2d_courant = connec_poly2d;

      PDM_g_num_t *numabs_tria_courant = numabs_tria;
      PDM_g_num_t *numabs_quad_courant = numabs_quad;
      PDM_g_num_t *numabs_poly2d_courant = numabs_poly2d;

      /* Construction de la connectivit� sommet-> arrete */

      PDM_l_num_t *connec_som_are = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 2 * geom->som[ipart]->n_som);

      PDM_l_num_t idx_tria   = 0;
      PDM_l_num_t idx_quad   = n_tria;
      PDM_l_num_t idx_poly2d = idx_quad + n_quad;

      for (int j = 0; j < 2 * geom->som[ipart]->n_som; j++) {
        connec_som_are[j] = -1;
      }

      for (int i = 0; i < n_cell_courant; i++) {

        PDM_l_num_t ideb = cell_face_idx_courant[i] -1;
        PDM_l_num_t n_face_cell = cell_face_nb_courant[i];
        PDM_l_num_t ifin = ideb + n_face_cell;
 
        for (int j = ideb; j < ifin; j++) {
          PDM_l_num_t ifac = cell_face_courant[j] - 1;
          PDM_l_num_t isom1 = face_som_courant[2*ifac] - 1;
          PDM_l_num_t isom2 = face_som_courant[2*ifac+1] - 1;

          if (connec_som_are[2*isom1] == -1)
            connec_som_are[2*isom1] = ifac;
          else
            connec_som_are[2*isom1+1] = ifac;

          if (connec_som_are[2*isom2] == -1)
            connec_som_are[2*isom2] = ifac;
          else
            connec_som_are[2*isom2+1] = ifac;
        }      

        PDM_l_num_t *connec_courant;
        if (n_face_cell == 3) {
          num_cell_parent_to_local_courant[i] = idx_tria++;
          *numabs_tria_courant = numabs_courant[i];
          numabs_tria_courant += 1;
          connec_courant = connec_tria_courant;
          connec_tria_courant += n_face_cell;
        }
        else if (n_face_cell == 4) {
          num_cell_parent_to_local_courant[i] = idx_quad++;;
          *numabs_quad_courant = numabs_courant[i];
          numabs_quad_courant += 1;
          connec_courant = connec_quad_courant;
          connec_quad_courant += n_face_cell;
        }
        else {
          num_cell_parent_to_local_courant[i] = idx_poly2d++;
          *numabs_poly2d_courant = numabs_courant[i];
          numabs_poly2d_courant += 1;
          connec_courant = connec_poly2d_courant;
          *connec_poly2d_idx_courant = *(connec_poly2d_idx_courant - 1) +  n_face_cell;
          connec_poly2d_idx_courant += 1;
          connec_poly2d_courant += n_face_cell;
        }

        /* Remplissage de la connectivite */
        
        PDM_l_num_t idx_som = 0;
        PDM_l_num_t face_courant = cell_face_courant[ideb] - 1;
        PDM_l_num_t isom1 = face_som_courant[2*face_courant] - 1;
        PDM_l_num_t isom_suiv = face_som_courant[2*face_courant + 1] - 1;
        connec_courant[idx_som++] = isom1 + 1;

        while (isom1 != isom_suiv) {
          assert(idx_som <= n_face_cell);
          connec_courant[idx_som++] = isom_suiv + 1;
          
          /* Face suivante */
          
          PDM_l_num_t face_suiv = connec_som_are[2*isom_suiv];
          if (face_suiv == face_courant)
            face_suiv = connec_som_are[2*isom_suiv + 1];
          face_courant = face_suiv;

          /* Sommet suivant */

          PDM_l_num_t isom_tmp = face_som_courant[2*face_courant] - 1;
          if (isom_tmp == isom_suiv)
            isom_tmp = face_som_courant[2*face_courant + 1] - 1;
          isom_suiv = isom_tmp;
        }
        
        for (int j= 0; j < n_face_cell; j++) {
          connec_som_are[2*(connec_courant[j] -1)] = - 1;
          connec_som_are[2*(connec_courant[j] -1) + 1] = - 1;
        }
      }
      free(connec_som_are);

      if (som_elts[0] > 0)
        PDM_writer_geom_bloc_std_set(id_cs,
                             id_geom,
                             id_bloc_tria3,
                             ipart,
                             n_tria,
                             connec_tria,
                             numabs_tria);

      if (som_elts[1] > 0)
        PDM_writer_geom_bloc_std_set(id_cs,
                             id_geom,
                             id_bloc_quad4,
                             ipart,
                             n_quad,
                             connec_quad,
                             numabs_quad);
    
      if (som_elts[2] > 0)
        PDM_writer_geom_bloc_poly2d_set(id_cs,
                                id_geom,
                                id_bloc_poly_2d,
                                ipart,
                                n_poly2d,
                                connec_poly2d_idx,
                                connec_poly2d,
                                numabs_poly2d);
    }
    if (geom->prepa_blocs != NULL) {
      free(geom->prepa_blocs->n_cell);
      free(geom->prepa_blocs->n_face);
      free(geom->prepa_blocs->n_tria);
      free(geom->prepa_blocs->n_quad);
      free(geom->prepa_blocs->n_poly2d);
      free(geom->prepa_blocs->l_connec_poly2d);
      free(geom->prepa_blocs->face_som_idx);
      free(geom->prepa_blocs->face_som_nb);
      free(geom->prepa_blocs->face_som);
      free(geom->prepa_blocs->cell_face_idx);
      free(geom->prepa_blocs->cell_face_nb);
      free(geom->prepa_blocs->cell_face);
      free(geom->prepa_blocs->add_etat);
      free(geom->prepa_blocs->numabs);
      free(geom->prepa_blocs);
      geom->prepa_blocs = NULL;
    }
  }
}


/*----------------------------------------------------------------------------
 *
 * Ajout de faces decrites en fonctions des sommets. Cette fonction
 * d�termine les types des �l�ments et cr�e des blocs regrouppant les �l�ments
 * de m�me type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   n_elt           <-- Nombre de cellules 3D ajout�es         
 *   n_face          <-- Nombre de faces d�crites               
 *   face_som_idx    <-- Index de connectivite faces -> sommets
 *   face_som        <-- Connectivite faces -> sommets                                       
 *   numabs          <-- Numerotation absolue des faces    
 *   ind_num         --> Indirection vers la nouvelle numerotation des faces    
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_faces_facesom_add, PDM_WRITER_GEOM_FACES_FACESOM_ADD)
(
int         *id_cs,
int         *id_geom,
int         *id_part,
int         *n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_g_num_t   *numabs
) 
{
  PDM_writer_geom_faces_facesom_add(*id_cs,
                            *id_geom,
                            *id_part, 
                            *n_face,
                            face_som_idx,
                            face_som_nb,
                            face_som,
                            numabs);
} 

void
PDM_writer_geom_faces_facesom_add
(
const int          id_cs,
const int          id_geom,
const int          id_part, 
const int          n_face,
PDM_l_num_t    *face_som_idx,
PDM_l_num_t    *face_som_nb,
PDM_l_num_t    *face_som,
PDM_g_num_t   *numabs
)
{
  /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);
  int n_part = 0;

  if (geom->num_cell_parent_to_local == NULL) {
    geom->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * geom->n_part); 
    for (int ipart = 0; ipart < geom->n_part; ipart++) { 
      geom->num_cell_parent_to_local[ipart] = NULL;
    }
  }

  geom->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_face);
  for (int i = 0; i < n_face; i++) {
    geom->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (geom->prepa_blocs == NULL) {
    geom->prepa_blocs = (PDM_writer_geom_prepa_blocs_t *) malloc(sizeof(PDM_writer_geom_prepa_blocs_t));
    geom->prepa_blocs->t_add = 3;
    geom->prepa_blocs->n_tria_proc = 0;    /* Nb de triangles par proc */
    geom->prepa_blocs->n_quad_proc = 0;    /* Nb de quads par proc */
    geom->prepa_blocs->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    geom->prepa_blocs->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->n_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->l_connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part); 
    geom->prepa_blocs->face_som_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part); 
    geom->prepa_blocs->face_som_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->face_som = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*geom->n_part);
    geom->prepa_blocs->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*geom->n_part);
    geom->prepa_blocs->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*geom->n_part);
    for (int i = 0; i < geom->n_part; i++) {
      geom->prepa_blocs->add_etat[i] = 0;
    }
  }

  if (geom->prepa_blocs->t_add != 3) {
    fprintf(stderr, "Erreur Cs_geom_cell2d_cellface_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  PDM_l_num_t n_tria    = 0;
  PDM_l_num_t n_quad    = 0;
  PDM_l_num_t n_poly2d  = 0;
  PDM_l_num_t l_connec_poly2d  = 0;

  for (int i = 0; i < n_face; i++) {

    PDM_l_num_t n_som_face = face_som_nb[i];
    if (n_som_face == 3)
      n_tria += 1;
    else if (n_som_face == 4)
      n_quad += 1;
    else {
      n_poly2d  += 1;
      l_connec_poly2d += n_som_face;
    }
  }

  geom->prepa_blocs->n_tria_proc           += n_tria;
  geom->prepa_blocs->n_quad_proc           += n_quad;
  geom->prepa_blocs->n_poly2d_proc         += n_poly2d;
  geom->prepa_blocs->add_etat[id_part]      = 1;
  geom->prepa_blocs->n_tria[id_part]        = n_tria; 
  geom->prepa_blocs->n_quad[id_part]        = n_quad; 
  geom->prepa_blocs->n_poly2d[id_part]      = n_poly2d;
  geom->prepa_blocs->l_connec_poly2d[id_part] = l_connec_poly2d;
  geom->prepa_blocs->face_som_idx[id_part]  = face_som_idx;
  geom->prepa_blocs->face_som_nb[id_part]   = face_som_nb;
  geom->prepa_blocs->face_som[id_part]      = face_som;
  geom->prepa_blocs->numabs[id_part]        = numabs;
  geom->prepa_blocs->add_etat[id_part]      = 1;
  geom->prepa_blocs->n_face[id_part]        = n_face;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < geom->n_part; i++) {
    if (geom->prepa_blocs->add_etat[i] == 1)
      n_part += 1;
  }

  if (geom->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[3];
    PDM_l_num_t som_elts[3];

    elts[0] = geom->prepa_blocs->n_tria_proc;
    elts[1] = geom->prepa_blocs->n_quad_proc;
    elts[2] = geom->prepa_blocs->n_poly2d_proc;
    
    PDM_MPI_Allreduce(elts, som_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, cs->pdm_mpi_comm);

    int id_bloc_tria3;
    int id_bloc_quad4;
    int id_bloc_poly_2d;

    if (som_elts[0] > 0)
      id_bloc_tria3 = PDM_writer_geom_bloc_add(id_cs,
                                       id_geom,
                                       PDM_writer_ON,
                                       PDM_writer_TRIA3);

    if (som_elts[1] > 0)
      id_bloc_quad4 = PDM_writer_geom_bloc_add(id_cs,
                                       id_geom,
                                       PDM_writer_ON,
                                       PDM_writer_QUAD4);
    
    if (som_elts[2] > 0)
      id_bloc_poly_2d = PDM_writer_geom_bloc_add(id_cs,
                                         id_geom,
                                         PDM_writer_ON,
                                         PDM_writer_POLY_2D);

    /* Determination de la connectivite de chaque element */

    for (int ipart = 0; ipart < geom->n_part; ipart++) {

      PDM_l_num_t *num_cell_parent_to_local_courant = geom->num_cell_parent_to_local[ipart];
      PDM_l_num_t *face_som_idx_courant = geom->prepa_blocs->face_som_idx[ipart];
      PDM_l_num_t *face_som_nb_courant = geom->prepa_blocs->face_som_nb[ipart];
      PDM_l_num_t *face_som_courant = geom->prepa_blocs->face_som[ipart];
      PDM_g_num_t *numabs_courant = geom->prepa_blocs->numabs[ipart];
 
      n_tria   = geom->prepa_blocs->n_tria[ipart];
      n_quad    = geom->prepa_blocs->n_quad[ipart];
      n_poly2d  = geom->prepa_blocs->n_poly2d[ipart];
      l_connec_poly2d  = geom->prepa_blocs->l_connec_poly2d[ipart];

      PDM_l_num_t *connec_tria = NULL;
      PDM_l_num_t *connec_quad = NULL;
      PDM_l_num_t *connec_poly2d = NULL;
      PDM_l_num_t *connec_poly2d_idx = NULL;
 
      PDM_g_num_t *numabs_tria = NULL;
      PDM_g_num_t *numabs_quad = NULL;
      PDM_g_num_t *numabs_poly2d = NULL;

      if (n_tria > 0) {
        connec_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 3 *n_tria);
        numabs_tria = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tria);
      }

      if (n_quad > 0) {
        connec_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 * n_quad);
        numabs_quad = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_quad);
      }

      if (n_poly2d > 0) {
        connec_poly2d_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly2d + 1));
        connec_poly2d_idx[0] = 0;
        connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_connec_poly2d);
        numabs_poly2d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly2d);
      }

      PDM_l_num_t *connec_tria_courant = connec_tria;
      PDM_l_num_t *connec_quad_courant = connec_quad;
      PDM_l_num_t *connec_poly2d_idx_courant = connec_poly2d_idx + 1;
      PDM_l_num_t *connec_poly2d_courant = connec_poly2d;

      PDM_g_num_t *numabs_tria_courant = numabs_tria;
      PDM_g_num_t *numabs_quad_courant = numabs_quad;
      PDM_g_num_t *numabs_poly2d_courant = numabs_poly2d;

      PDM_l_num_t idx_tria   = 0;
      PDM_l_num_t idx_quad   = n_tria;
      PDM_l_num_t idx_poly2d = idx_quad + n_quad;

      for (int i = 0; i < n_face; i++) {
        PDM_l_num_t n_som_face = face_som_nb_courant[i];
        PDM_l_num_t idx_som_face = face_som_idx_courant[i] - 1;
        PDM_l_num_t *connec_courant;

        if (n_som_face == 3) {
          num_cell_parent_to_local_courant[i] = idx_tria++;
          *numabs_tria_courant = numabs_courant[i];
          numabs_tria_courant += 1;
          connec_courant = connec_tria_courant;
          connec_tria_courant += n_som_face;
        }
        else if (n_som_face == 4) {
          num_cell_parent_to_local_courant[i] = idx_quad++;;
          *numabs_quad_courant = numabs_courant[i];
          numabs_quad_courant += 1;
          connec_courant = connec_quad_courant;
          connec_quad_courant += n_som_face;
        }
        else { 
          num_cell_parent_to_local_courant[i] = idx_poly2d++;
          *numabs_poly2d_courant = numabs_courant[i];
          numabs_poly2d_courant += 1;
          *connec_poly2d_idx_courant = *(connec_poly2d_idx_courant - 1) + n_som_face;
          connec_poly2d_idx_courant += 1;
          connec_courant = connec_poly2d_courant;
          connec_poly2d_courant += n_som_face;
        }

        /* Remplissage de la connectivite */

        for (int j = 0; j < n_som_face; j++)
          connec_courant[j] = face_som_courant[idx_som_face++];
      }

      if (som_elts[0] > 0) 
        PDM_writer_geom_bloc_std_set(id_cs,
                             id_geom,
                             id_bloc_tria3,
                             ipart,
                             n_tria,
                             connec_tria,
                             numabs_tria);

      if (som_elts[1] > 0)
        PDM_writer_geom_bloc_std_set(id_cs,
                             id_geom,
                             id_bloc_quad4,
                             ipart,
                             n_quad,
                             connec_quad,
                             numabs_quad);
    
      if (som_elts[2] > 0) 
        PDM_writer_geom_bloc_poly2d_set(id_cs,
                                id_geom,
                                id_bloc_poly_2d,
                                ipart,
                                n_poly2d,
                                connec_poly2d_idx,
                                connec_poly2d,
                                numabs_poly2d);
    }
    if (geom->prepa_blocs != NULL) {
      free(geom->prepa_blocs->n_face);
      free(geom->prepa_blocs->n_tria);
      free(geom->prepa_blocs->n_quad);
      free(geom->prepa_blocs->n_poly2d);
      free(geom->prepa_blocs->l_connec_poly2d);
      free(geom->prepa_blocs->face_som_idx);
      free(geom->prepa_blocs->face_som_nb);
      free(geom->prepa_blocs->face_som);
      free(geom->prepa_blocs->add_etat);
      free(geom->prepa_blocs->numabs);
      free(geom->prepa_blocs);
      geom->prepa_blocs = NULL;
    }
  }
} 

/*----------------------------------------------------------------------------
 * Ecriture du maillage courant                                  
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_write, PDM_WRITER_GEOM_WRITE)
(
int           *id_cs,
int           *id_geom
) 
{
  PDM_writer_geom_write (*id_cs,
               *id_geom); 
}

void
PDM_writer_geom_write
(
const int            id_cs,
const int            id_geom
) 
{

    /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);

  //TODO  faire un retour si g�om�trie n'est pas dependante du temps
  //       et si on n'est pas au premier incr�ment
  /* Mise a jour du nombre total d'�l�ments */

   PDM_g_num_t n_elt_abs = geom->n_elt_abs;
  
  // validation gnum
//  int id = PDM_gnum_create (3, geom->n_part, geom->pdm_mpi_comm);
//  printf("1\n");
//
//  for (int i = 0; i < geom->n_part; i++) {
//    
//    PDM_gnum_set_from_coords (id, i, geom->som[i]->n_som, geom->som[i]->_coords);
//
//  }
//  
//  printf("2\n");
//  PDM_gnum_compute (id);
//    printf("3\n");
//
//  for (int i = 0; i < geom->n_part; i++) {
//    
//    geom->som[i]->_numabs = PDM_gnum_get (id, i);
//    printf ("pp %ld\n", PDM_gnum_get (id, i));
//    
//  }
//    printf("4\n");
//
//  PDM_gnum_free (id, 1);
//    printf("5\n");

  // fin validation gnum

  PDM_MPI_Allreduce(&n_elt_abs, &geom->n_elt_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, geom->pdm_mpi_comm);

  /* D�termination de la num�rotation absolue interne des elements 
     Independante du parallelisme */

  _calcul_numabs_elt(geom);

  /* Decoupage des polygones */

  if (geom->st_decoup_poly2d) {
    _split_poly2d(geom);
    _calcul_numabs_split_poly2d(geom);
  }

  /* Decoupage des polyedres */

  if (geom->st_decoup_poly3d) {
    _split_poly3d(geom);
    _calcul_numabs_split_poly3d(geom);
  }

  /* Ecriture au format */

  /* Appel de la fonction complementaire propre au format */
    
  if (fmt_tab[cs->fmt_id]->geom_write_fct != NULL) {
    (fmt_tab[cs->fmt_id]->geom_write_fct) (geom);
  }

}


/*----------------------------------------------------------------------------
 * Liberation des donnees decrivant le maillage courant
 *  On conserve uniquement les donn�es sur les indirections vers la num�rotation
 *  absolue
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_free, PDM_WRITER_GEOM_FREE)
(
int           *id_cs,
int           *id_geom
) 
{
  PDM_writer_geom_free(*id_cs,
                       *id_geom); 

}

void
PDM_writer_geom_free
(
const int      id_cs,
const int      id_geom
) 
{
  PDM_writer_geom_data_free(id_cs,
                   id_geom);

  PDM_writer_t *cs = _PDM_writer_get(id_cs);
  
  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);

  /* Lib�ration des sommets */
  
  if (geom->som != NULL) {
    for (int i = 0; i < geom->n_part; i++) {
      geom->som[i] = _som_free (geom->som[i]);
    }

    free(geom->som);
    geom->som = NULL;
  }

  /* Boucle sur les blocs standard */

  for (int i = 0; i < geom->n_blocs_std; i++) {
    PDM_writer_bloc_std_t *_bloc_std = geom->blocs_std[i];
    _bloc_std_free(_bloc_std);
  }

  if (geom->blocs_std != NULL) {
    free(geom->blocs_std);
    geom->blocs_std = NULL;
  }
  geom->n_blocs_std = 0;

  /* Boucle sur les blocs de polygones */ 

  for (int i = 0; i < geom->n_blocs_poly2d; i++) {
    PDM_writer_bloc_poly2d_t *_bloc_poly2d = geom->blocs_poly2d[i];
    _bloc_poly2d_free(_bloc_poly2d);
  }

  if (geom->blocs_poly2d != NULL) {
    free(geom->blocs_poly2d);
    geom->blocs_poly2d = NULL;
  }
  geom->n_blocs_poly2d  = 0;

  /* Boucle sur les blocs de polyedres */ 

  for (int i = 0; i < geom->n_blocs_poly3d; i++) {
    PDM_writer_bloc_poly3d_t *_bloc_poly3d = geom->blocs_poly3d[i];
    _bloc_poly3d_free(_bloc_poly3d);
  }

  if (geom->blocs_poly3d != NULL) {
    free(geom->blocs_poly3d);
    geom->blocs_poly3d = NULL;
  }
  geom->n_blocs_poly3d = 0;

  /* Lib�ration de la structure */ 

  if (geom->num_cell_parent_to_local != NULL) {
    for (int ipart = 0; ipart < geom->n_part; ipart++) {
      if (geom->num_cell_parent_to_local[ipart] != NULL)
        free(geom->num_cell_parent_to_local[ipart]);
    }
    free(geom->num_cell_parent_to_local);
    geom->num_cell_parent_to_local = NULL;
  }

  free(geom->n_cell);
  geom->n_cell = NULL;

  free(geom->nom_geom);

  /* Lib�ration sp�cifique au format */

  /* Appel de la fonction complementaire propre au format */
    
  if (fmt_tab[cs->fmt_id]->geom_free_fct != NULL) {
    (fmt_tab[cs->fmt_id]->geom_free_fct) (geom);
  }

  free(geom);

  cs->geom_tab[id_geom] = NULL;
  cs->n_geom_tab -= 1;

  if (cs->n_geom_tab == 0) {
    free(cs->geom_tab);
    cs->geom_tab = NULL;
    cs->l_geom_tab = 0;
  }
}


/*----------------------------------------------------------------------------
 * Liberation partielle des donnees decrivant le maillage courant
 * les indirections sur les num�rotation absolues sont conserv�es
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_geom_data_free, PDM_WRITER_GEOM_DATA_FREE)
(
int           *id_cs,
int           *id_geom
) 
{
  PDM_writer_geom_data_free(*id_cs,
                   *id_geom); 
}

void
PDM_writer_geom_data_free
(
const int      id_cs,
const int      id_geom
) 
{
  /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);
  
  /* Boucle sur les blocs standard */

  for (int i = 0; i < geom->n_blocs_std; i++) {
    PDM_writer_bloc_std_t *_bloc_std = geom->blocs_std[i];
    _bloc_std_free_partial(_bloc_std);
  }

  /* Boucle sur les blocs de polygones */ 

  for (int i = 0; i < geom->n_blocs_poly2d; i++) {
    PDM_writer_bloc_poly2d_t *_bloc_poly2d = geom->blocs_poly2d[i];
    _bloc_poly2d_free_partial(_bloc_poly2d);
  }

  /* Boucle sur les blocs de polyedres */ 

  for (int i = 0; i < geom->n_blocs_poly3d; i++) {
    PDM_writer_bloc_poly3d_t *_bloc_poly3d = geom->blocs_poly3d[i];
    _bloc_poly3d_free_partial(_bloc_poly3d);
  }
}

/*----------------------------------------------------------------------------
 * Creation d'une variable                                                     
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   st_dep_temps    <-- Indique si la variable est dependante du temps
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   t_var           <-- Type de variable
 *   nom_var         <-- Nom de la variable
 *
 * return :
 *                   --> Identificateur de l'objet variable     
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_create_cf, PDM_WRITER_VAR_CREATE_CF)
(
int         *id_cs,
int         *st_dep_tps,
int         *dim,  
int         *loc,   
char        *nom_var,   
int         *l_nom_var,   
int         *id_var 
ARGF_SUPP_CHAINE
)
{
  char *nom_var_c = PDM_fortran_to_c_string(nom_var, *l_nom_var);

  *id_var = PDM_writer_var_create (*id_cs,
                                   (PDM_writer_statut_t) *st_dep_tps,
                                   (PDM_writer_var_dim_t) *dim,   
                                   (PDM_writer_var_loc_t) *loc,   
                                   nom_var_c);

  free(nom_var_c);
}

int
PDM_writer_var_create
(
const int          id_cs,
const PDM_writer_statut_t  st_dep_tps,
const PDM_writer_var_dim_t dim,  
const PDM_writer_var_loc_t loc,  
const char        *nom_var
)
{
  /* Recherche de l'objet cs courant */

  PDM_writer_t *cs = _PDM_writer_get (id_cs);

  /* Mise a jour du tableau de stockage */

  if (cs->var_tab == NULL) {
    cs->l_var_tab = 4;
    cs->var_tab = (PDM_writer_var_t **) malloc(cs->l_var_tab * sizeof(PDM_writer_var_t *));
    for (int i = 0; i < cs->l_var_tab; i++) 
      cs->var_tab[i] = NULL;
  } 
  
  if (cs->l_var_tab <= cs->n_var_tab) {
    int p_l_var_tab = cs->l_var_tab;
    cs->l_var_tab = 2 * cs->l_var_tab;
    cs->var_tab = (PDM_writer_var_t**) realloc((void*) cs->var_tab, cs->l_var_tab * sizeof(PDM_writer_var_t *));
    
    for (int i = p_l_var_tab; i < cs->l_var_tab; i++) 
      cs->var_tab[i] = NULL;
  }

  /* Recherche de la premiere place libre pour stocker le bloc */

  int id_var = 0;
  while (cs->var_tab[id_var] != NULL) 
    id_var++;

  /* Allocation de la structure PDM_writer_var_t */

  PDM_writer_var_t *var = (PDM_writer_var_t *) malloc(sizeof(PDM_writer_var_t));

  cs->n_var_tab += 1;

  cs->var_tab[id_var] = var;

  /* Initialisation de la structure PDM_writer_var_t */

  _var_init (var);

  size_t l_nom_var = strlen(nom_var);
  var->nom_var = (char *) malloc(sizeof(char) * (l_nom_var + 1));
  strcpy(var->nom_var, nom_var);   /* Nom de la variable */

  var->st_dep_tps = st_dep_tps;    /* Variable en temps */
  var->dim        = dim;           /* Dimension de la variable */
  var->loc        = loc;           /* Dimension de la variable */
  var->_cs        = cs;

  /* Appel de la fonction complementaire propre au format */

  if (fmt_tab[cs->fmt_id]->var_create_fct != NULL) {
    (fmt_tab[cs->fmt_id]->var_create_fct) (var);
  }

  return id_var;
}


/*----------------------------------------------------------------------------
 * Ecriture des valeurs de la variable
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part         <-- Identificateur de la partition dans l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise � jour
 *   val             <-- Valeurs
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_write, PDM_WRITER_VAR_WRITE)
(
int         *id_cs,
int         *id_var
)
{
  PDM_writer_var_write(*id_cs,
             *id_var);
}

void
PDM_writer_var_write
(
const int        id_cs,
const int        id_var 
)
{

  /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_var_t *var = _var_get(cs, id_var);

  /* Ecriture au format */

  if (fmt_tab[cs->fmt_id]->var_write_fct != NULL) {
    (fmt_tab[cs->fmt_id]->var_write_fct) (var);
  }
  
}



/*----------------------------------------------------------------------------
 * Mise a jour des valeurs de la variable
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_part         <-- Identificateur de la partition dans l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise � jour
 *   val             <-- Valeurs
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_set, PDM_WRITER_VAR_SET)
(
int         *id_cs,
int         *id_var,
int         *id_geom,
int         *id_part,
PDM_real_t   *val
)
{
  PDM_writer_var_set(*id_cs,
             *id_var,
             *id_geom,
             *id_part,
             val);
}

void
PDM_writer_var_set
(
const int        id_cs,
const int        id_var,
const int        id_geom,
const int        id_part,
const PDM_real_t *val
)
{

  /* Acces a l'objet de geometrie courant */


  PDM_writer_t *cs = _PDM_writer_get(id_cs);

  PDM_writer_var_t *var = _var_get(cs, id_var);

  PDM_writer_geom_t *geom = _geom_get(cs, id_geom);



  if (var->_val == NULL) {
    var->_val = (double ***) malloc(sizeof(double **) * cs->l_geom_tab);
    for (int i = 0; i < cs->l_geom_tab; i++)
      var->_val[i] = NULL;
  }

  if (cs->l_geom_tab <= id_geom) {
    fprintf(stderr, "Erreur cs_var_set    : Indice de geometrie incorrect\n");
    abort();
  }

  if (var->_val[id_geom] == NULL) {
    var->_val[id_geom] = (double **) malloc(sizeof(double *) * geom->n_part);
    for (int i = 0; i < geom->n_part; i++)
      var->_val[id_geom][i] = NULL;
  }

  double **val_geom = var->_val[id_geom];

  if (geom->n_part <= id_part) {
    fprintf(stderr, "Erreur cs_var_set    : Indice de partition incorrect\n");
    abort();
  }

  if (var->loc == PDM_WRITER_VAR_ELEMENTS) {
    val_geom[id_part] = (double *) malloc(sizeof(double) * var->dim * geom->n_cell[id_part]);
    if (geom->num_cell_parent_to_local != NULL) {
      for (int i = 0; i < geom->n_cell[id_part]; i++) {
        for (int j = 0; j < var->dim; j++)
          val_geom[id_part][var->dim * geom->num_cell_parent_to_local[id_part][i]+j] = val[i*var->dim + j];
      }
    }
    else {
      for (int i = 0; i < geom->n_cell[id_part]; i++) {
        for (int j = 0; j < var->dim; j++)
          val_geom[id_part][var->dim * i+j] = val[i*var->dim + j];
      }
    }
  }
  else {
    val_geom[id_part] = (double *) malloc(sizeof(double) * var->dim * geom->som[id_part]->n_som);
    for (int i = 0; i < geom->som[id_part]->n_som * var->dim; i++) {
      val_geom[id_part][i] = val[i];
    }
  }

}


/*----------------------------------------------------------------------------
 * Liberation du tableau de donnees des variables
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise � jour
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_data_free, PDM_WRITER_VAR_DATA_FREE)
(
int         *id_cs,
int         *id_var
)
{
  PDM_writer_var_data_free(*id_cs,
                  *id_var);
}

void
PDM_writer_var_data_free
(
const int    id_cs,
const int    id_var
)
{
  /* Acces a l'objet de geometrie courant */

  PDM_writer_t *cs = _PDM_writer_get (id_cs);

  PDM_writer_var_t *var = _var_get(cs, id_var);

  if (var->_val != NULL) {
    for (int i = 0; i < cs->l_geom_tab; i++) {
      PDM_writer_geom_t    *geom = cs->geom_tab[i];
      if ((geom != NULL) && (var->_val[i] != NULL)) {
        for (int j = 0; j < geom->n_part; j++) {
          if (var->_val[i][j] != NULL)
            free(var->_val[i][j]);
          var->_val[i][j] = NULL;
        }
      }
    }
  }
}


/*----------------------------------------------------------------------------
 * Liberation du tableau de donnees des variables
 *
 * parameters :
 *   id_cs           <-- Identificateur de l'objet cs
 *   id_geom         <-- Identificateur de l'objet geometrique
 *   id_var          <-- Identificateur de la variable mise � jour
 *
 *----------------------------------------------------------------------------*/

void
PROCF (pdm_writer_var_free, PDM_WRITER_VAR_FREE)
(
int         *id_cs,
int         *id_var
)
{
  PDM_writer_var_free (*id_cs,
              *id_var);
}

void
PDM_writer_var_free
(
const int    id_cs,
const int    id_var
)
{
  PDM_writer_t *cs = _PDM_writer_get (id_cs);

  if (cs->var_tab != NULL) {

    PDM_writer_var_data_free(id_cs,
                    id_var);

    /* Acces a l'objet de geometrie courant */

    PDM_writer_var_t *var = _var_get (cs, id_var);

    free(var->nom_var);

    for (int i = 0; i < cs->l_geom_tab; i++) {
      if (var->_val[i] != NULL)
        free(var->_val[i]);
      var->_val[i] = NULL;
    }
    free(var->_val);
    var->_val = NULL;

    /* Lib�ration sp�cifique au format */


    if (fmt_tab[cs->fmt_id]->var_free_fct != NULL) {
      (fmt_tab[cs->fmt_id]->var_free_fct) (var);
    }

    free(cs->var_tab[id_var]);
    cs->var_tab[id_var] = NULL;
    cs->n_var_tab -= 1;

    if (cs->n_var_tab == 0) {
      free(cs->var_tab);
      cs->var_tab = NULL;
      cs->l_var_tab = 0;
    }
  }
}

/**
 * \brief Add a writer format
 *
 * Define a new format writer
 *
 * \param [in] name            Name                                                    
 * \param [in] create_fct      Customize \ref PDM_writer_create function for the new format  (or NULL)
 * \param [in] free_fct        Customize \ref PDM_writer_free function for the new format (or NULL)
 * \param [in] beg_step_fct    Customize \ref PDM_writer_step_beg function for the new format (or NULL)
 * \param [in] end_step_fct    Customize \ref PDM_writer_step_end function for the new format (or NULL)
 * \param [in] geom_create_fct Customize \ref PDM_writer_geom_create function for the new format (or NULL)
 * \param [in] geom_write_fct  Customize \ref PDM_writer_geom_write function for the new format
 * \param [in] geom_free_fct   Customize \ref PDM_writer_geom_free function for the new format (or NULL)
 * \param [in] var_create_fct  Customize \ref PDM_writer_var_create function for the new format (or NULL)
 * \param [in] var_write_fct   Customize \ref PDM_writer_var_write function for the new format
 * \param [in] var_free_fct    Customize \ref PDM_writer_var_free function for the new format (or NULL)
 *
 */

void
PDM_writer_fmt_add
(
 const char                  *name,           /*!< Name                                                     */
 const PDM_writer_fct_t      create_fct,      /*!< Customize \ref PDM_writer_create function for the format */
 const PDM_writer_fct_t      free_fct,        /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t      beg_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t      end_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t geom_create_fct, /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t geom_write_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t geom_free_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t  var_create_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t  var_write_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t  var_free_fct    /*!< Customize \ref PDM_writer_free function for the format   */
)
{
  _load_intern_fmt();
  if (n_fmt_tab >= l_fmt_tab) {
    l_fmt_tab *= 2;
    fmt_tab = realloc (fmt_tab, sizeof(PDM_writer_fmt_t *) * l_fmt_tab);
  }

  if (geom_write_fct == NULL) {
    fprintf (stderr, "Error PDM_writer_fmt_add : Undefined geom write function\n");
    abort ();
  }

  if (var_write_fct == NULL) {
    fprintf (stderr, "Error PDM_writer_fmt_add : Undefined var write function\n");
    abort ();
  }
  
  fmt_tab[n_fmt_tab] = malloc (sizeof(PDM_writer_fmt_t));
  
  fmt_tab[n_fmt_tab]->name            = name;
  fmt_tab[n_fmt_tab]->create_fct      = create_fct;
  fmt_tab[n_fmt_tab]->free_fct        = free_fct;
  fmt_tab[n_fmt_tab]->beg_step_fct    = beg_step_fct;
  fmt_tab[n_fmt_tab]->end_step_fct    = end_step_fct;
  fmt_tab[n_fmt_tab]->geom_create_fct = geom_create_fct;
  fmt_tab[n_fmt_tab]->geom_write_fct  = geom_write_fct;
  fmt_tab[n_fmt_tab]->geom_free_fct   = geom_free_fct;
  fmt_tab[n_fmt_tab]->var_create_fct  = var_create_fct;
  fmt_tab[n_fmt_tab]->var_write_fct   = var_write_fct; 
  fmt_tab[n_fmt_tab]->var_free_fct    = var_free_fct;

  n_fmt_tab += 1;

}


/**
 * \brief Free formats
 *
 */

void
PDM_writer_fmt_free
(
 void
)
{
  if (fmt_tab != NULL) {

    for (int i = 0; i < l_fmt_tab; i++) {
      if (fmt_tab != NULL) {
        free (fmt_tab[i]);
      }
    }

    l_fmt_tab = 0;
    n_fmt_tab = 0;

    free (fmt_tab);
    fmt_tab = NULL;
    
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
