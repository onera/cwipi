#ifndef __PDM_WRITER_PRIV_H__
#define __PDM_WRITER_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_writer.h"
#include "pdm_io.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Definitions des macro
 *============================================================================*/

/*============================================================================
 * Definition des types 
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Description des sommets d'une partition
 *----------------------------------------------------------------------------*/

typedef struct PDM_writer_som_t PDM_writer_som_t;

struct PDM_writer_som_t {

  PDM_writer_som_t         *parent;   /* sommets parents */
  PDM_l_num_t          n_som;  /* Nombre de sommets de la partition courante */
  double           *coords;  /* Coordonnees des sommets de la partition courante allouee dans le cas d'une
                                definition a partir d'une geometrie parente */
  const double    *_coords;  /* Coordonnees des sommets de la partition courante (mapping memoire) */
  const PDM_g_num_t *_numabs;  /* Numerotation absolue des sommets de la partition courante (mapping memoire) */
  const int       *_numparent; /* Numerotation dans la geometrie parente (mapping memoire) */
};

/*----------------------------------------------------------------------------
 * Description d'un bloc geometrique
 *----------------------------------------------------------------------------*/

typedef struct {

  PDM_writer_elt_geom_t     t_elt;             /* Type d'elements dans le blocs */
  PDM_writer_statut_t       st_free_data;       /* Liberation des donnees a la destruction */
  PDM_l_num_t          n_part;            /* Nombre de partitions */
  PDM_l_num_t         *n_elt;             /* Nombre d'elements */
  PDM_l_num_t        **_connec;           /* Connectivite des elements (mapping memoire) */
  PDM_g_num_t       **_numabs;           /* Numerotation initiale absolue 
                                          (mapping memoire) */
  PDM_l_num_t        **_num_part;         /* Numerotation initiale dans la partition
                                          (mapping memoire) */
  PDM_g_num_t       **numabs_int;        /* Numerotation absolue interne au bloc 
                                          de chaque element de chaque partition  */

} PDM_writer_bloc_std_t;

/*----------------------------------------------------------------------------
 * Description des sommets supplementaires 
 * lies aux decoupages des polyedres/polygones
 *----------------------------------------------------------------------------*/

typedef struct PDM_writer_som_sup_t{

  PDM_l_num_t       n_part;             /* Nombre de partitions */
  double       **coords;             /* Coordonnees des sommets de chaque partition */
  PDM_l_num_t      *n_som;              /* Nombre de sommets de chaque partition */
  PDM_g_num_t    **numabs;             /* Numerotation absolue des sommets de chaque partition */

} PDM_writer_som_sup_t;

/*----------------------------------------------------------------------------
 * Description d'un bloc geometrique de type polygone
 *----------------------------------------------------------------------------*/

typedef struct PDM_writer_bloc_poly2d_t{

  PDM_writer_statut_t      st_free_data;       /* Liberation des donnees a la destruction */
  PDM_writer_statut_t st_decoup_poly2d;   /* Statut du decoupage du bloc en triangles */
  PDM_l_num_t    n_part;             /* Nombre de partitions */
  PDM_l_num_t   *n_elt;              /* Nombre d'elements de chaque partition */
  PDM_l_num_t  **_connec_idx;        /* Index des elments dans la connectivite de chaque partition 
                                     (mapping memoire) */
  PDM_l_num_t  **_connec;            /* Connectivite des elements de chaque partition
                                     (mapping memoire) */
  PDM_l_num_t  **_num_part;          /* Numerotation initiale dans la partition
                                     (mapping memoire) */
  PDM_g_num_t **_numabs;            /* Numerotation absolue de chaque partition
                                     (mapping memoire) */

  PDM_g_num_t **numabs_int;         /* Numerotation absolue interne au bloc de chaque element
                                     de chaque partition */
  
  PDM_l_num_t  **tri_idx;            /* Index d'indirection des polygones 
                                     vers les triangles de chaque partition 
                                     si decoupage des polygones */

  PDM_writer_bloc_std_t  *bloc_tri;       /* Bloc de triangles 
                                     si decoupage des polygones */
  PDM_l_num_t  **quad_idx;           /* Index d'indirection des polygones 
                                     vers les quadrangles de chaque partition 
                                     si decoupage des polygones */

  PDM_writer_bloc_std_t  *bloc_quad;      /* Bloc de quadrangles
                                     si decoupage des polygones */

  PDM_writer_som_sup_t  *som_sup;         /* Sommets supplementaires pour chaque partition
                                     si decoupage des polygones */

} PDM_writer_bloc_poly2d_t;

/*----------------------------------------------------------------------------
 * Description d'un bloc geometrique de type polyedre
 *----------------------------------------------------------------------------*/

typedef struct PDM_writer_bloc_poly3d_t{

  PDM_writer_statut_t  st_free_data;       /* Liberation des donnees a la destruction */
  PDM_writer_statut_t st_decoup_poly3d;   /* Statut du decoupage du bloc en tetraedres */
  PDM_l_num_t    n_part;             /* Nombre de partitions */
  PDM_l_num_t   *n_elt;              /* Nombre d'elements de chaque partition */
  PDM_l_num_t   *n_face;             /* Nombre de faces de chaque polyedre de chaque partition 
                                     (mapping memoire) */
  PDM_l_num_t  **_facsom_idx;        /* Index de connectivite des faces de la partition
                                     (mapping memoire) */

  PDM_l_num_t  **_facsom;            /* Connectivite des faces                    
                                     de chaque polyedre de chaque partition 
                                     (mapping memoire) */
  PDM_l_num_t  **_cellfac_idx;       /* Index de la connectivite cell->face 
                                     de tous les polyedres de chaque partition   
                                     (mapping memoire) */
  PDM_l_num_t  **_cellfac;           /* Connectivite cell->face 
                                     de tous les polyedres de chaque partition   
                                     (mapping memoire) */
  PDM_l_num_t  **_num_part;          /* Numerotation initiale dans la partition
                                     (mapping memoire) */
  PDM_g_num_t **_numabs;            /* Numerotation absolue des polyedres   
                                     (mapping memoire) */

  PDM_g_num_t **numabs_int;         /* Numerotation absolue interne au bloc de chaque element
                                     de chaque partition */
  PDM_l_num_t  **tetra_idx;          /* Index d'indirection des polyedres 
                                     vers les tetraedes de chaque partition */
  PDM_writer_bloc_std_t  *bloc_tetra;     /* Bloc de tetraedres */
  PDM_l_num_t  **pyra_idx;           /* Index d'indirection des polyedres 
                                     vers les tetraedes de chaque partition */
  PDM_writer_bloc_std_t  *bloc_pyra;      /* Bloc de pyramides */

  PDM_writer_som_sup_t  *som_sup;         /* Sommets supplementaires pour chaque partition
                                     si decoupage des polyedres */
} PDM_writer_bloc_poly3d_t;

/*----------------------------------------------------------------------------
 * Etat de preparation des blocs 
 *----------------------------------------------------------------------------*/

typedef struct PDM_writer_geom_prepa_blocs_t{

  PDM_l_num_t  n_tria_proc;    /* Nb de triangles par proc */
  PDM_l_num_t  n_quad_proc;    /* Nb de quads par proc */
  PDM_l_num_t  n_poly2d_proc;  /* Nb de poly2d par proc */
  PDM_l_num_t  n_tetra_proc;   /* Nb de tetra par proc */
  PDM_l_num_t  n_hexa_proc;    /* Nb d'hexa par proc */
  PDM_l_num_t  n_prism_proc;   /* Nb de prisme par proc */
  PDM_l_num_t  n_pyramid_proc; /* Nb de pyramide par proc */
  PDM_l_num_t  n_poly3d_proc;  /* Nb de poly3d par proc */
  PDM_l_num_t *add_etat;       /* Permet de verifier que toutes */
  PDM_l_num_t  t_add;          /* Type d'ajout (1 : cell3d_cellface,
                                             2 : cell2d_cellface,
                                             3 : faces_facesom_add) */
  PDM_l_num_t  *n_tetra;
  PDM_l_num_t  *n_hexa;
  PDM_l_num_t  *n_prism;
  PDM_l_num_t  *n_pyramid;
  PDM_l_num_t  *n_poly3d;
  PDM_l_num_t  *n_tria;
  PDM_l_num_t  *n_quad;
  PDM_l_num_t  *n_poly2d;
  PDM_l_num_t  *l_connec_poly2d;
  PDM_l_num_t  **face_som_idx; 
  PDM_l_num_t  **face_som_nb;
  PDM_l_num_t  **face_som;
  PDM_l_num_t  *n_cell;
  PDM_l_num_t  *n_face;
  PDM_l_num_t  **cell_face_idx;
  PDM_l_num_t  **cell_face_nb;
  PDM_l_num_t  **cell_face;
  PDM_g_num_t **numabs;

} PDM_writer_geom_prepa_blocs_t;

/*----------------------------------------------------------------------------
 * Description de la geometrie
 *----------------------------------------------------------------------------*/

struct _PDM_writer_geom_t {

  char               *nom_geom;           /* Nom de la geometrie */
  PDM_writer_geom_connec_t    t_connec;           /* Type de connectivité */
  PDM_writer_statut_t         st_decoup_poly2d;   /* Decoupage des polygones */
  PDM_writer_statut_t         st_decoup_poly3d;   /* Decoupage des polyedres */
  PDM_g_num_t           n_som_abs;          /* Nombre absolu de sommets */
  PDM_g_num_t           n_som_abs_total;    /* Nombre absolu de sommets 
                                             (sommets issus du decoupage 
                                             des polyedres/polygones compris) */
  PDM_g_num_t           n_elt_abs;          /* Nombre absolu d'elements */
  int                 n_part;             /* Nombre de partitions */
  PDM_writer_som_t           **som;                /* Description des sommmets de chaque partition */
  PDM_l_num_t           *n_cell;             /* Nombre de blocs d'elements standard */
  PDM_l_num_t            n_blocs_std;        /* Nombre de blocs d'elements standard */
  PDM_l_num_t            l_blocs_std;        /* Taille du tableau de blocs d'elements standard */
  PDM_writer_bloc_std_t     **blocs_std;          /* Blocs d'elements standard */
  PDM_l_num_t            n_blocs_poly2d;     /* Nombre de blocs de polygones */
  PDM_l_num_t            l_blocs_poly2d;     /* Taille du tableau de blocs de polygones */  
  PDM_writer_bloc_poly2d_t  **blocs_poly2d;       /* Blocs de polygones */
  PDM_l_num_t            n_blocs_poly3d;     /* Nombre de blocs de polyedres */
  PDM_l_num_t            l_blocs_poly3d;     /* Taille du tableau de polyedres */
  PDM_writer_bloc_poly3d_t  **blocs_poly3d;       /* Blocs de polyedres */
  void               *geom_fmt;           /* Description propre au format fmt */
  PDM_writer_t        *_cs;                /* Pointeur sur la structure cs parente */
  PDM_MPI_Comm            pdm_mpi_comm;           /* Communicateur MPI */
  PDM_writer_geom_prepa_blocs_t *prepa_blocs;     /* Preparation des blocs */
  PDM_l_num_t         **num_cell_parent_to_local;/* Indirection de la numerotation des cellules
                                                initiale vers la numerotation locale
                                                imposee par les blocs */
} ;

/*----------------------------------------------------------------------------
 * Description de la variable
 *----------------------------------------------------------------------------*/

struct _PDM_writer_var_t{

  char          *nom_var;            /* Nom de la geometrie */
  PDM_writer_statut_t    st_dep_tps;         /* Variable en temps */
  PDM_writer_var_dim_t   dim;                /* Dimension de la variable */
  PDM_writer_var_loc_t   loc;                /* Localisation de la variable */
  double ***_val;                    /* Valeurs de la variable 
                                        (par partition) mapping mémoire */ 
  PDM_writer_t   *_cs;                /* Pointeur sur la structure cs parente */
  void          *var_fmt;            /* Description propre au format fmt */

} ;


/*----------------------------------------------------------------------------
 * Type Cedre sortie 
 *----------------------------------------------------------------------------*/

struct _PDM_writer_t {

  int                    fmt_id;        /* Format de la sortie */
  PDM_writer_fmt_fic_t   fmt_fic;    /* Format du fichier ascii ou binaire */
  PDM_writer_topologie_t topologie;  /* Type de toplogie du maillage */
  PDM_writer_statut_t    st_reprise; /* Reprise d'une sortie existante */
  char          *rep_sortie; /* Nom du repertoire de sortie */
  char          *nom_sortie; /* Nom de la sortie */
  PDM_MPI_Comm       pdm_mpi_comm;   /* Communicateur MPI */
  void          *sortie_fmt; /* Description propre au format */    
  PDM_writer_var_t     **var_tab;    /* Tableau des variables */
  int            l_var_tab;  /* Taille du tableau des variables */
  int            n_var_tab;  /* Nombre de variables dans le tableau des variables */
  PDM_writer_geom_t    **geom_tab;   /* Tableau des geometries */
  int            l_geom_tab; /* Taille du tableau des geometries */
  int            n_geom_tab; /* Nombre de geometries dans le tableau des geoemtries */
  double         physical_time; /* Temps physique de la simulation */
  PDM_io_acces_t acces;    /* Type d'acces au fichier (MPIIIO,...) */
  double         prop_noeuds_actifs; /* Proportion des noeuds actifs */
};



/**
 * \struct PDM_writer_fmt_t
 * \brief  Writer format
 *
 */

typedef struct PDM_writer_fmt_t {

  const char           *name;           /*!< Name                                                     */
  PDM_writer_fct_t      create_fct;      /*!< Customize \ref PDM_writer_create function for the format */
  PDM_writer_fct_t      free_fct;        /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_fct_t      beg_step_fct;    /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_fct_t      end_step_fct;    /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_geom_fct_t geom_create_fct; /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_geom_fct_t geom_write_fct;  /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_geom_fct_t geom_free_fct;   /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_var_fct_t  var_create_fct;  /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_var_fct_t  var_write_fct;   /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_var_fct_t  var_free_fct;    /*!< Customize \ref PDM_writer_free function for the format   */
  
} PDM_writer_fmt_t;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_WRITER_PRIV_H__ */
