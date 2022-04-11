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
#include "pdm_mesh_nodal.h"

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
 * Description de la geometrie
 *----------------------------------------------------------------------------*/

struct _PDM_writer_geom_t {

  char                     *nom_geom;           /* Nom de la geometrie */
  void                     *geom_fmt;           /* Description propre au format fmt */
  PDM_writer_t             *_cs;                /* Pointeur sur la structure cs parente */
  PDM_MPI_Comm              pdm_mpi_comm;       /* Communicateur MPI */
  PDM_Mesh_nodal_t         *mesh_nodal;         /* Mesh handle */

};

typedef struct _PDM_writer_geom_tab_t {

  int                 n_geom;
  int                 s_geom;
  PDM_writer_geom_t **geom;

} _PDM_writer_geom_tab_t;

/*----------------------------------------------------------------------------
 * Mapping des noms de variable
 *----------------------------------------------------------------------------*/

typedef struct PDM_writer_name_map_t {

  char *public_name;         /* Nom public */
  char *private_name;        /* Nom privé */

} PDM_writer_name_map_t;


typedef struct _PDM_writer_name_map_tab_t {

  int                     n_name_map;
  int                     s_name_map;
  PDM_writer_name_map_t **name_map;

} _PDM_writer_name_map_tab_t;

/*----------------------------------------------------------------------------
 * Description d'une option : couple nom/valeur
 *----------------------------------------------------------------------------*/

typedef struct {

  char *nom;
  char *val;

} PDM_writer_option_t;

/*----------------------------------------------------------------------------
 * Description de la variable
 *----------------------------------------------------------------------------*/

struct _PDM_writer_var_t{

  char                  *nom_var;        /* Nom de la geometrie */
  PDM_writer_status_t    st_dep_tps;     /* Variable en temps */
  PDM_writer_var_dim_t   dim;            /* Dimension de la variable */
  PDM_writer_var_loc_t   loc;            /* Localisation de la variable */
  double              ***_val;           /* Valeurs de la variable
                                            (par partition) mapping mémoire */
  PDM_writer_t          *_cs;            /* Pointeur sur la structure cs parente */
  void                  *var_fmt;        /* Description propre au format fmt */
  char                  *private_name;   /* Nom privé de la variable (si mapping) */

} ;


typedef struct _PDM_writer_var_tab_t {

  int                 n_var;
  int                 s_var;
  PDM_writer_var_t **var;

} _PDM_writer_var_tab_t;


/*----------------------------------------------------------------------------
 * Type Cedre sortie
 *----------------------------------------------------------------------------*/

struct _PDM_writer_t {

  int                         fmt_id;             /* Format de la sortie */
  PDM_writer_fmt_fic_t        fmt_fic;            /* Format du fichier ascii ou binaire */
  PDM_writer_topology_t      topologie;          /* Type de toplogie du maillage */
  PDM_writer_status_t         st_reprise;         /* Reprise d'une sortie existante */
  char                       *rep_sortie;         /* Nom du repertoire de sortie */
  char                       *nom_sortie;         /* Nom de la sortie */
  PDM_MPI_Comm                pdm_mpi_comm;       /* Communicateur MPI */
  void                       *sortie_fmt;         /* Description propre au format */
  _PDM_writer_var_tab_t      *var_tab;            /* Tableau des variables */
  _PDM_writer_geom_tab_t     *geom_tab;           /* Tableau des geometries */
  double                      physical_time;      /* Temps physique de la simulation */
  PDM_io_kind_t              acces;              /* Type d'acces au fichier (MPIIIO,...) */
  double                      prop_noeuds_actifs; /* Proportion des noeuds actifs */
  _PDM_writer_name_map_tab_t *name_map_tab;       /* Stockage du mapping des noms */
  int                         n_options;          /* Nombre d'options */
  PDM_writer_option_t        *options;            /* Options complementaire */

};


/**
 * \struct PDM_writer_fmt_t
 * \brief  Writer format
 *
 */

typedef struct PDM_writer_fmt_t {

  char                       *name;            /*!< Name                                                     */
  PDM_writer_fct_t      create_fct;      /*!< Customize \ref PDM_writer_create function for the format */
  PDM_writer_fct_t      free_fct;        /*!< Customize \ref PDM_writer_free function for the format   */
  PDM_writer_fct_t      beg_step_fct;    /*!< Customize \ref PDM_writer_beg_step function for the format   */
  PDM_writer_fct_t      end_step_fct;    /*!< Customize \ref PDM_writer_end_step function for the format   */
  PDM_writer_geom_fct_t geom_create_fct; /*!< Customize \ref PDM_writer_geom_create function for the format   */
  PDM_writer_geom_fct_t geom_write_fct;  /*!< Customize \ref PDM_writer_geom_write function for the format   */
  PDM_writer_geom_fct_t geom_free_fct;   /*!< Customize \ref PDM_writer_geom_free function for the format   */
  PDM_writer_var_fct_t  var_create_fct;  /*!< Customize \ref PDM_writer_var_create function for the format   */
  PDM_writer_var_fct_t  var_write_fct;   /*!< Customize \ref PDM_writer_var_write function for the format   */
  PDM_writer_var_fct_t  var_free_fct;    /*!< Customize \ref PDM_writer_var_free function for the format   */

} PDM_writer_fmt_t;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_WRITER_PRIV_H__ */
