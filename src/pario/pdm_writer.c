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
#include "pdm_remove_blank.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_mesh_nodal.h"


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

  PDM_WRITER_DEB_ID_BLOC_STD    = 0,
  PDM_WRITER_DEB_ID_BLOC_POLY2D = 1000000,
  PDM_WRITER_DEB_ID_BLOC_POLY3D = 2000000

} PDM_writer_deb_id_bloc_t;

/*============================================================================
 * Variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Stockage des objets cs
 *----------------------------------------------------------------------------*/

// static PDM_Handles_t *cs_tab = NULL;

/*----------------------------------------------------------------------------
 * Stockage des objets cs
 *----------------------------------------------------------------------------*/

static PDM_writer_fmt_t **fmt_tab = NULL;

/*----------------------------------------------------------------------------
 * Nombre d'objets cs stockes dans fmt_tab
 *----------------------------------------------------------------------------*/

static const int n_intern_fmt = 1;
static       int s_fmt_tab    = 0;
static       int n_fmt_tab    = 0;

/*============================================================================
 * Definition des fonctions privees
 *============================================================================*/

/**
 *
 * \brief Create a \ref _PDM_writer_geom_tab_t object
 *
 * \param [in]  size   Initial size of the array of geometries
 *
 * \return    Pointer to a new \ref _PDM_writer_geom_tab_t object
 *
 */

static _PDM_writer_geom_tab_t *
_pdm_writer_geom_tab_create
(
 const int size
 )
{
  _PDM_writer_geom_tab_t *geom_tab = (_PDM_writer_geom_tab_t *) malloc(sizeof(_PDM_writer_geom_tab_t));

  geom_tab->n_geom = 0;
  geom_tab->s_geom = size;

  geom_tab->geom = (PDM_writer_geom_t **) malloc(sizeof(PDM_writer_geom_t *) * geom_tab->s_geom);
  for (int i = 0; i < geom_tab->s_geom; i++) {
    geom_tab->geom[i] = NULL;
  }

  return geom_tab;
}


/**
 *
 * \brief Add a geometry
 *
 * \param [in] geom_tab   Pointer to \ref _PDM_writer_geom_tab_t object
 * \param [in] geom       Pointer to \ref PDM_writer_geom_t object
 *
 */

static int
_pdm_writer_geom_tab_add
(
 _PDM_writer_geom_tab_t *geom_tab,
 PDM_writer_geom_t      *geom
 )
{
  assert (geom_tab != NULL);

  if (geom_tab->n_geom >= geom_tab->s_geom) {
    geom_tab->s_geom = PDM_MAX(2*geom_tab->s_geom, geom_tab->n_geom+1);
    geom_tab->geom = (PDM_writer_geom_t **) realloc(geom_tab->geom, sizeof(PDM_writer_geom_t *) * geom_tab->s_geom);

    for (int i = geom_tab->n_geom+1; i < geom_tab->s_geom; i++) {
      geom_tab->geom[i] = NULL;
    }
  }

  int id_geom = geom_tab->n_geom;
  geom_tab->n_geom++;

  geom_tab->geom[id_geom] = geom;


  return id_geom;
}


/**
 *
 * \brief Create a \ref _PDM_writer_var_tab_t object
 *
 * \param [in]  size   Initial size of the array of variables
 *
 * \return    Pointer to a new \ref _PDM_writer_var_tab_t object
 *
 */

static _PDM_writer_var_tab_t *
_pdm_writer_var_tab_create
(
 const int size
 )
{
  _PDM_writer_var_tab_t *var_tab = (_PDM_writer_var_tab_t *) malloc(sizeof(_PDM_writer_var_tab_t));

  var_tab->n_var = 0;
  var_tab->s_var = size;

  var_tab->var = (PDM_writer_var_t **) malloc(sizeof(PDM_writer_var_t *) * var_tab->s_var);
  for (int i = 0; i < var_tab->s_var; i++) {
    var_tab->var[i] = NULL;
  }

  return var_tab;
}


/**
 *
 * \brief Add a variable
 *
 * \param [in] var_tab   Pointer to \ref _PDM_writer_var_tab_t object
 * \param [in] var       Pointer to \ref PDM_writer_var_t object
 *
 */

static int
_pdm_writer_var_tab_add
(
 _PDM_writer_var_tab_t *var_tab,
 PDM_writer_var_t      *var
 )
{
  assert (var_tab != NULL);

  if (var_tab->n_var >= var_tab->s_var) {
    var_tab->s_var = PDM_MAX(2*var_tab->s_var, var_tab->n_var+1);
    var_tab->var = (PDM_writer_var_t **) realloc(var_tab->var, sizeof(PDM_writer_var_t *) * var_tab->s_var);

    for (int i = var_tab->n_var+1; i < var_tab->s_var; i++) {
      var_tab->var[i] = NULL;
    }
  }

  int id_var = var_tab->n_var;
  var_tab->n_var++;

  var_tab->var[id_var] = var;


  return id_var;
}


/**
 * \brief Initialize a \ref PDM_writer_geom_t object
 *
 * \param [in] geom            Pointer to \ref PDM_writer_geom_t object
 * \param [in] n_part          Number of partition
 * \param [in] comm            MPI communicator
 *
 */

static void
_geom_init
(
PDM_writer_geom_t  *geom,
const int           n_part,
const PDM_MPI_Comm  comm
)
{
  geom->nom_geom = NULL;

  geom->st_decoup_poly2d = PDM_WRITER_OFF;
  geom->st_decoup_poly3d = PDM_WRITER_OFF;

  geom->mesh_nodal = PDM_Mesh_nodal_create (n_part, comm);

  geom->geom_fmt       = NULL;
  geom->_cs            = NULL;
  geom->pdm_mpi_comm   = comm;
}


/**
 * \brief Initialize a \ref PDM_writer_var_t object
 *
 * \param [in] var            Pointer to \ref PDM_writer_var_t object
 *
 */

static void
_var_init
(
PDM_writer_var_t *var
)
{
  var->nom_var    = NULL;                    /* Nom de la geometrie */
  var->st_dep_tps = PDM_WRITER_OFF;          /* Variable en temps */
  var->dim        = PDM_WRITER_VAR_CSTE;     /* Dimension de la variable */
  var->loc        = PDM_WRITER_VAR_SOMMETS;  /* Dimension de la variable */
  var->_val       = NULL;                    /* Valeurs de la variable */
  var->var_fmt    = NULL;                    /* Description propre au format fmt */
  var->_cs        = NULL;
}


/**
 *
 * \brief Parse la chaine options pour construire la structure CS correspondante
 *
 * \param [in]  options_str           Options_str : chaine provenant de cs_cree
 * \param [out] n_options             Nombre d'options trouvees
 * \param [out] options               Liste des options parsees
 *
 */

static void
_parse_options
(
 const char           *options_str,
 int                  *n_options,
 PDM_writer_option_t **options
)
{

  if (options_str == NULL) {
    *n_options = 0;
    *options = NULL;
  }

  char *_options_str  = malloc (sizeof(char) * (strlen(options_str) + 1));
  strcpy(_options_str, options_str);

  *n_options = 0;
  char *str2 = _options_str;
  char *pch;

  do {
    pch = strtok (str2,"=");
    str2 = NULL;
    if (pch != NULL) {
      pch = strtok (str2, ":");
      if (pch == NULL) {
        PDM_error(__FILE__, __LINE__, 0, "CS_cree : Erreur dans le parsing des options specifiques :"
                 "verifier les separateurs dans la chaine 'options'\n");
        exit(1);
      }
      *n_options += 1;
    }
  } while (pch != NULL);

  strcpy(_options_str, options_str);
  str2 = _options_str;
  *options = malloc (sizeof(PDM_writer_option_t) * (*n_options));
  PDM_writer_option_t *_curr = *options;

  do {
    pch = strtok (str2,"=");
    str2 = NULL;
    if (pch != NULL) {
      _curr->nom = PDM_remove_blank (pch);
      pch = strtok (str2, ":");
      if (pch != NULL) {
        _curr->val = PDM_remove_blank (pch);
      }
    }
    _curr += 1;
  } while (pch != NULL);


  free (_options_str);
}

/**
 *
 * \brief Load built-in formats
 *
 */

static void
_load_intern_fmt (void)
{
  if (fmt_tab != NULL) {
    return;
  }

  s_fmt_tab = 2 * n_intern_fmt;
  n_fmt_tab = 0;
  fmt_tab = (PDM_writer_fmt_t **) malloc (sizeof(PDM_writer_fmt_t *) * s_fmt_tab);

  /* Ensight */

  PDM_writer_fmt_t *fmt = malloc (sizeof(PDM_writer_fmt_t));
  fmt->name = malloc (sizeof(int) * 8);
  strcpy (fmt->name, "Ensight");
  fmt->create_fct       = PDM_writer_ensight_create;
  fmt->free_fct         = PDM_writer_ensight_free;
  fmt->beg_step_fct     = PDM_writer_ensight_step_beg;
  fmt->end_step_fct     = PDM_writer_ensight_step_end;
  fmt->geom_create_fct  = PDM_writer_ensight_geom_create;
  fmt->geom_write_fct   = PDM_writer_ensight_geom_write;
  fmt->geom_free_fct    = PDM_writer_ensight_geom_free;
  fmt->var_create_fct   = PDM_writer_ensight_var_create;
  fmt->var_write_fct    = PDM_writer_ensight_var_write;
  fmt->var_free_fct     = PDM_writer_ensight_var_free;

  fmt_tab[n_fmt_tab++] = fmt;
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/**
 *
 * \brief Cree un objet CS (Cedre Sortie) et retoure un pointeur sur cet objet
 *
 * \param [in] fmt                  Format de sortie
 * \param [in] fmt_fic              Binary or ASCII
 * \param [in] topologie            Indique le maillage est mobile ou non
 * \param [in] st_reprise           Complete les sorties des calculs precedents en reprise
 * \param [in] rep_sortie           Repertoire de sortie
 * \param [in] nom_sortie           Nom de la sortie
 * \param [in] pdm_mpi_com          Communicateur MSG
 * \param [in] acces                Type d'acces
 * \param [in] prop_noeuds_actifs   Proportion des noeuds actifs dans les acces au fichier
 *                                    *  -1 : tous les processus actifs
 *                                    *   1 : un processus par noeud
 *                                    * 0 < val < 1 : un processus par noeud actif
 * \param [in] options              Options complementaires propres au format sous
 *                                 la forme ("nom_1 = val_1 : ... : nom_n = val_n")
 *
 * \return   Pointer to \ref PDM_writer object
 *
 */
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
const PDM_MPI_Fint  *pdm_mpi_comm,
const int           *acces,
const double        *prop_noeuds_actifs,
const char          *options,
const int           *l_options,
int                 *id_cs
ARGF_SUPP_CHAINE
)
{
  char *rep_sortie_c        = PDM_fortran_to_c_string(rep_sortie, *l_rep_sortie);
  char *nom_sortie_c        = PDM_fortran_to_c_string(nom_sortie, *l_nom_sortie);
  char *fmt_c        = PDM_fortran_to_c_string(fmt, *l_fmt);
  char *options_c           = NULL;
  if (*l_options > 0)
    options_c = PDM_fortran_to_c_string(options, *l_options);

  const PDM_MPI_Comm pdm_mpi_comm_c = PDM_MPI_Comm_f2c(*pdm_mpi_comm);

  // *id_cs = PDM_writer_create(fmt_c,
  //                            (PDM_writer_fmt_fic_t)   *fmt_fic,
  //                            (PDM_writer_topologie_t) *topologie,
  //                            (PDM_writer_statut_t)    *st_reprise,
  //                            rep_sortie_c,
  //                            nom_sortie_c,
  //                            pdm_mpi_comm_c,
  //                            (PDM_io_acces_t) *acces,
  //                            *prop_noeuds_actifs,
  //                             options_c);

  if (rep_sortie_c != NULL) {
    free(rep_sortie_c);
  }

  if (nom_sortie_c != NULL) {
    free(nom_sortie_c);
  }

  if (options_c != NULL) {
    free(options_c);
  }

  if (fmt_c != NULL) {
    free(fmt_c);
  }
}

PDM_writer_t *
PDM_writer_create
(
const char                   *fmt,
const PDM_writer_fmt_fic_t    fmt_fic,
const PDM_writer_topologie_t  topologie,
const PDM_writer_statut_t     st_reprise,
const char                   *rep_sortie,
const char                   *nom_sortie,
const PDM_MPI_Comm            pdm_mpi_comm,
const PDM_io_acces_t          acces,
const double                  prop_noeuds_actifs,
const char                   *options
)
{

  if (fmt_tab == NULL) {
    _load_intern_fmt();
  }

  /* Look for fmt */

  int fmt_id = -1;

  for (int i = 0; i < n_fmt_tab; i++) {
    PDM_writer_fmt_t *fmt_ptr = fmt_tab[i];
    if (!strcmp(fmt, fmt_ptr->name)) {
      fmt_id = i;
      break;
    }
  }

  if (fmt_id == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_create : unknown format '%s'", fmt);
    abort();
  }

  /* Mise a jour du tableau de stockage */

  PDM_io_mkdir(rep_sortie);

  // if (cs_tab == NULL) {
  //   cs_tab = PDM_Handles_create (4);
  // }

  /* Creation du repertoire de sortie si non cree */

#ifdef _WIN32
  mkdir(rep_sortie);
#else
  mkdir(rep_sortie, 0775);
#endif

  /* Allocation de la structure PDM_writer_t */

  PDM_writer_t *cs = (PDM_writer_t *) malloc(sizeof(PDM_writer_t));

  // int id_cs = PDM_Handles_store (cs_tab, (void *) cs);

  /* Initialisation de la structure PDM_writer_t */

  cs->fmt_id      = fmt_id;     /* Format de sortie */
  cs->fmt_fic     = fmt_fic;    /* Format du fichier de sortie */
  cs->topologie   = topologie;  /* Type de toplogie du maillage */
  cs->st_reprise  = st_reprise; /* Reprise d'une sortie existante */

  size_t l_rep_sortie = strlen(rep_sortie);
  cs->rep_sortie = (char *) malloc(sizeof(char) * (l_rep_sortie + 1));
  strcpy(cs->rep_sortie, rep_sortie);  /* Nom du repertoire de sortie */
  // Gestion des options

  cs->n_options = 0;
  cs->options = NULL;
  if (options != NULL) {
    _parse_options (options, &(cs->n_options), &(cs->options));
  }

  size_t l_nom_sortie = strlen(nom_sortie);
  cs->nom_sortie = (char *) malloc(sizeof(char) * (l_nom_sortie + 1));
  strcpy(cs->nom_sortie, nom_sortie);  /* Nom de la sortie */

  cs->pdm_mpi_comm    = pdm_mpi_comm;  /* Communicateur MPI */
  cs->sortie_fmt  = NULL;      /* Pointeur sur l'objet sortie de format fmt */

  cs->var_tab     = NULL;      /* Tableau des variables */
  cs->geom_tab    = NULL;      /* Tableau des geometries */
  cs->physical_time = 0;       /* Temps physique de simulation */
  cs->acces       = acces;
  cs->prop_noeuds_actifs = prop_noeuds_actifs;
  cs->name_map   = NULL;

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->create_fct != NULL) {
    (fmt_ptr->create_fct) (cs);
  }

  return cs;

}

/**
 * \brief Libere un objet CS (Cedre Sortie) et retourne un pointeur NULL si pas d'erreur
 *
 * \param [in] cs    Pointer to \ref PDM_writer object
 *
 */
void
PROCF (pdm_writer_free, PDM_WRITER_FREE)
(
int        *id_cs
)
{
  // PDM_writer_free(*id_cs);
}

void
PDM_writer_free
(
 PDM_writer_t *cs
)
{

  /* Recherche de l'objet cs courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get(cs_tab, id_cs);

  if (cs == NULL) {
    return;
  }

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->free_fct != NULL) {
    (fmt_ptr->free_fct) (cs);
  }

  /* Liberation des differents elements de la structure */

  free(cs->rep_sortie);
  free(cs->nom_sortie);

  /* Liberation des variables */

  if (cs->var_tab != NULL) {
    if (cs->var_tab->var != NULL) {
      for (int i = 0; i < cs->var_tab->n_var; i++) {
        if (cs->var_tab->var[i] != NULL) {
          PDM_writer_var_free(cs, i);
        }
      }

      free(cs->var_tab->var);
      cs->var_tab->var = NULL;
    }

    free(cs->var_tab);
    cs->var_tab = NULL;
  }

  if (cs->options != NULL) {
    for (int i = 0; i < cs->n_options; i++) {
      if ((cs->options[i]).nom != NULL) {
        free ((cs->options[i]).nom);
      }
      if ((cs->options[i]).val != NULL) {
        free ((cs->options[i]).val);
      }
    }
    free (cs->options);
  }


  /* Liberation de la geometrie */

  if (cs->geom_tab != NULL) {
    if (cs->geom_tab->geom != NULL) {
      for (int i = 0; i < cs->geom_tab->n_geom; i++) {
        if (cs->geom_tab->geom[i] != NULL) {
          PDM_writer_geom_free(cs, i);
        }
      }

      free(cs->geom_tab->geom);
      cs->geom_tab->geom = NULL;
    }

    free(cs->geom_tab);
    cs->geom_tab = NULL;
  }

  if (cs->name_map != NULL) {

    int n_map_tab = PDM_Handles_n_get (cs->name_map);
    const int *map_index = PDM_Handles_idx_get(cs->name_map);

    while (n_map_tab > 0) {
      PDM_writer_name_map_t *map = (PDM_writer_name_map_t *)
              PDM_Handles_get (cs->name_map, map_index[0]);
      if (map != NULL) {
        free (map->public_name);
        free (map->private_name);
        free (map);
      }
      PDM_Handles_handle_free (cs->name_map, map_index[0], PDM_FALSE);
      n_map_tab = PDM_Handles_n_get (cs->name_map);
    }

    cs->name_map = PDM_Handles_free (cs->name_map);

  }

  /* Liberation de la structure */

  free(cs);
  cs = NULL;

  // PDM_Handles_handle_free (cs_tab, id_cs, PDM_FALSE);
  //
  // int n_cs = PDM_Handles_n_get (cs_tab);
  //
  // if (n_cs == 0) {
  //
  //   cs_tab = PDM_Handles_free (cs_tab);
  //   if (n_intern_fmt == n_fmt_tab) {
  //     PDM_writer_fmt_free();
  //   }
  // }

}


/**
 * \brief Debut d'increment
 *
 * \param [in] cs             Pointer to \ref PDM_writer object
 * \param [in] physical_time  Temps
 */
void
PROCF (pdm_writer_step_beg, PDM_WRITER_STEP_BEG)
(
int           *id_cs,
double        *physical_time
)
{
  // PDM_writer_step_beg(*id_cs,
  //             *physical_time);
}

void
PDM_writer_step_beg
(
 PDM_writer_t  *cs,
 const double   physical_time
)
{
  /* Recherche de l'objet cs courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  cs->physical_time = physical_time;

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->beg_step_fct != NULL) {
    (fmt_ptr->beg_step_fct) (cs);
  }

}

/**
 * \brief Fin d'increment
 *
 * \param [in] cs             Pointer to \ref PDM_writer object
 *
 */
void
PROCF (pdm_writer_step_end, PDM_WRITER_STEP_END)
(
int          *id_cs
)
{
  // PDM_writer_step_end(*id_cs);
}

void
PDM_writer_step_end
(
 PDM_writer_t  *cs
)
{
  /* Recherche de l'objet cs courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->end_step_fct != NULL) {
    (fmt_ptr->end_step_fct) (cs);
  }

}

/**
 * \brief Cree une nouvelle geometrie dans l'objet CS (Cedre Sortie)
 *
 * \param [in]  cs                Pointer to \ref PDM_writer object
 * \param [in]  nom_geom          Nom de l'objet geometrique
 * \param [in]  st_decoup_poly2d  Active le decoupage des polygones
 * \param [in]  st_decoup_poly3d  Active le decoupage des polyedres
 *
 * \return   Identificateur de l'objet geom dans cs
 *
 */
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

  // *id_geom = PDM_writer_geom_create(*id_cs,
  //                         nom_geom_c,
  //                         (PDM_writer_statut_t) *st_decoup_poly2d,
  //                         (PDM_writer_statut_t) *st_decoup_poly3d,
  //                         *n_part);

  if (nom_geom_c != NULL) {
    free(nom_geom_c);
  }
}

int
PDM_writer_geom_create
(
 PDM_writer_t               *cs,
 const char                 *nom_geom,
 const PDM_writer_statut_t   st_decoup_poly2d,
 const PDM_writer_statut_t   st_decoup_poly3d,
 const int                   n_part
)
{
  /* Erreur si le decoupage des polygones ou polyedres est choisi */

  if (n_part <= 0) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_geom_create : Le nombre de partition doit etre >\n"
                    "                      Ajuster le communicateur MPI ou\n"
                    "                      Creer un sous-domaine avec 0 element\n");
  }

  if ((st_decoup_poly2d == 1) || (st_decoup_poly3d == 1)) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_geom_create : Les fonctions de decoupage ne sont pas operationnelles\n");
    abort();
  }

  /* Recherche de l'objet cs courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->geom_tab == NULL) {
    cs->geom_tab = _pdm_writer_geom_tab_create(4);
  }

  /* Allocation de la structure PDM_writer_geom_t */

  PDM_writer_geom_t *geom = (PDM_writer_geom_t *) malloc(sizeof(PDM_writer_geom_t));

  int id_geom = _pdm_writer_geom_tab_add(cs->geom_tab, geom);


  /* Initialisation de la structure PDM_writer_geom_t */

  _geom_init(geom, n_part, cs->pdm_mpi_comm);

  geom->_cs = cs;
  geom->pdm_mpi_comm = cs->pdm_mpi_comm;
  size_t l_nom_geom = strlen(nom_geom);
  geom->nom_geom = (char *) malloc(sizeof(char) * (l_nom_geom + 1));
  strcpy(geom->nom_geom, nom_geom);  /* Nom de la geometrie */

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->geom_create_fct != NULL) {
    (fmt_ptr->geom_create_fct) (geom);
  }

  return id_geom;
}

//-->>
int
PDM_writer_geom_create_from_mesh_nodal
(
 PDM_writer_t              *cs,
 const char                *nom_geom,
 const PDM_writer_statut_t  st_decoup_poly2d,
 const PDM_writer_statut_t  st_decoup_poly3d,
 PDM_Mesh_nodal_t          *mesh
)
{
  /* Erreur si le decoupage des polygones ou polyedres est choisi */

  if ((st_decoup_poly2d == 1) || (st_decoup_poly3d == 1)) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_geom_create : Les fonctions de decoupage ne sont pas operationnelles\n");
    abort();
  }

  /* Recherche de l'objet cs courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->geom_tab == NULL) {
    cs->geom_tab = _pdm_writer_geom_tab_create(4);
  }

  /* Allocation de la structure PDM_writer_geom_t */

  PDM_writer_geom_t *geom = (PDM_writer_geom_t *) malloc(sizeof(PDM_writer_geom_t));

  int id_geom = _pdm_writer_geom_tab_add(cs->geom_tab, geom);

  /* Initialisation de la structure PDM_writer_geom_t */

  //_geom_init(geom, n_part, cs->pdm_mpi_comm);
  geom->nom_geom = NULL;
  geom->st_decoup_poly2d = PDM_WRITER_OFF;
  geom->st_decoup_poly3d = PDM_WRITER_OFF;
  geom->mesh_nodal = mesh;
  geom->geom_fmt       = NULL;

  geom->_cs = cs;
  geom->pdm_mpi_comm = cs->pdm_mpi_comm;
  size_t l_nom_geom = strlen(nom_geom);
  geom->nom_geom = (char *) malloc(sizeof(char) * (l_nom_geom + 1));
  strcpy(geom->nom_geom, nom_geom);  /* Nom de la geometrie */

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->geom_create_fct != NULL) {
    (fmt_ptr->geom_create_fct) (geom);
  }

  return id_geom;
}
//<<--


/**
 * \brief Definition des coordonnees de la partition courante
 *
 * \param [in] cs        Pointer to \ref PDM_writer object
 * \param [in] id_geom   Identificateur de l'objet geometrique
 * \param [in] id_part   Indice de partition
 * \param [in] n_som     Nombre de sommets de la partition
 * \param [in] coords    Coordonnes des sommets
 * \param [in] numabs    Numerotation absolue des sommets
 *
 */
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
  // PDM_writer_geom_coord_set(*id_cs,
  //                   *id_geom,
  //                   *id_part,
  //                   *n_som,
  //                   coords,
  //                   numabs);
}

void
PDM_writer_geom_coord_set
(
 PDM_writer_t      *cs,
 const int          id_geom,
 const int          id_part,
 const int          n_som,
 const PDM_real_t  *coords,
 const PDM_g_num_t *numabs
)
{

  /* Acces aux sommets de la partition */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_Mesh_nodal_coord_set (geom->mesh_nodal, id_part, n_som, coords, numabs);

  if (0 == 1) {
    printf("n_vtx : %d\n", n_som);
    for (int i = 0; i < n_som; i++) {
      printf ("%d "PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e\n", i+1, numabs[i],
              coords[3*i], coords[3*i+1], coords[3*i+2]);
    }
  }

}


/**
 * \brief Definition des coordonnees des sommets de la partition courante a partir
 *
 *
 * \param [in] cs               Pointer to \ref PDM_writer object
 * \param [in] id_geom          Identificateur de l'objet geometrique
 * \param [in] id_part          Indice de partition
 * \param [in] n_som            Nombre de sommets de la partition
 * \param [in] n_som_parent     Nombre de sommets parent
 * \param [in] numabs           Numerotation absolue des sommets (size = n_som)
 * \param [in] num_parent       Numerotation des sommets dans la numerotation parente (size = n_som)
 * \param [in] coords_parent    Coordonnes des sommets parents (size = 3 * n_som_parent)
 * \param [in] numabs_parent    Numerotation absolue des sommets parents (size = n_som_parent)
 *
 */
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
  // PDM_writer_geom_coord_from_parent_set (*id_cs,
  //                                *id_geom,
  //                                *id_part,
  //                                *n_som,
  //                                *n_som_parent,
  //                                numabs,
  //                                num_parent,
  //                                coords_parent,
  //                                numabs_parent);
}

void
PDM_writer_geom_coord_from_parent_set
(
 PDM_writer_t      *cs,
 const int          id_geom,
 const int          id_part,
 const int          n_som,
 const int          n_som_parent,
 const PDM_g_num_t *numabs,
 const int         *num_parent,
 const PDM_real_t  *coords_parent,
 const PDM_g_num_t *numabs_parent
)
{

  /* Acces aux sommets de la partition */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
  }

  PDM_Mesh_nodal_coord_from_parent_set (geom->mesh_nodal,
                                        id_part,
                                        n_som,
                                        n_som_parent,
                                        numabs,
                                        num_parent,
                                        coords_parent,
                                        numabs_parent);
}

/**
 * \brief Ajout d'un bloc d'elements d'un type donne
 *
 * \param [in] cs             Pointer to \ref PDM_writer object
 * \param [in] id_geom        Identificateur de l'objet geometrique
 * \param [in] t_elt          Type d'element
 *
 * \return   Identificateur du bloc
 *
 */
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
  // *id_bloc = PDM_writer_geom_bloc_add(*id_cs,
  //                             *id_geom,
  //                             *st_free_data,
  //                             (PDM_writer_elt_geom_t) *t_elt);
}

int
PDM_writer_geom_bloc_add
(
 PDM_writer_t                *cs,
 const int                    id_geom,
 PDM_writer_statut_t          st_free_data,
 const PDM_writer_elt_geom_t  t_elt
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  int id_block = PDM_Mesh_nodal_block_add (geom->mesh_nodal, (PDM_bool_t) st_free_data,
                                           (PDM_Mesh_nodal_elt_t) t_elt);


  return id_block;

}


/**
 * \brief Ajout d'un bloc d'elements d'un type donne dans la partition courante
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
 * \param [in] cs                  Pointer to \ref PDM_writer object
 * \param [in] id_geom             Identificateur de l'objet geometrique
 * \param [in] id_bloc             Identificateur du bloc
 * \param [in] id_part             Indice de partition
 * \param [in] t_elt               Type d'element
 * \param [in] n_elt               Nombre d'elements dans le bloc
 * \param [in] connec              Table de connectivite des elements
 * \param [in] num_part            Numerotation dans la partition
 *
 */
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
  // PDM_writer_geom_bloc_std_set (*id_cs,
  //                       *id_geom,
  //                       *id_bloc,
  //                       *id_part,
  //                       *n_elt,
  //                       connec,
  //                       numabs);
}

void
PDM_writer_geom_bloc_std_set
(
 PDM_writer_t  *cs,
 const int      id_geom,
 const int      id_bloc,
 const int      id_part,
 const int      n_elt,
 PDM_l_num_t   *connec,
 PDM_g_num_t   *numabs
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_Mesh_nodal_block_std_set (geom->mesh_nodal, id_bloc, id_part,
                                n_elt, connec, numabs, NULL);

}


/**
 * \brief Ajout d'un bloc de polygones dans la partition courante
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] id_part         Indice de partition
 * \param [in] n_elt           Nombre d'elements dans le bloc
 * \param [in] connec_idx      Index dans la table de connectivite (dim = n_elt+1)
 * \param [in] connec          Table de connectivite des elements (dim = connec_idx[n_elt])
 * \param [in] numabs          Numerotation absolue des elements
 *
 */
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
  // PDM_writer_geom_bloc_poly2d_set (*id_cs,
  //                          *id_geom,
  //                          *id_bloc,
  //                          *id_part,
  //                          *n_elt,
  //                          connec_idx,
  //                          connec,
  //                          numabs);
}

void
PDM_writer_geom_bloc_poly2d_set
(
PDM_writer_t        *cs,
const int            id_geom,
const int            id_bloc,
const int            id_part,
const PDM_l_num_t    n_elt,
      PDM_l_num_t   *connec_idx,
      PDM_l_num_t   *connec,
      PDM_g_num_t   *numabs
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_Mesh_nodal_block_poly2d_set (geom->mesh_nodal, id_bloc, id_part,
                                n_elt, connec_idx, connec, numabs, NULL);

}


/**
 * \brief Ajout d'un bloc de polyedres dans la partition courante
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] id_part         Indice de partition
 * \param [in] n_elt           Nombre d'elements dans le bloc
 * \param [in] n_face          Nombre de faces de chaque element (dim = n_elt)
 * \param [in] facsom_idx      Index dans la table de connectivite des faces (dim = n_face_total+1)
 * \param [in] facsom          Table de connectivite des faces (dim = facsom_idx[n_face_total}
 * \param [in] cellfac_idx     Index dans la table de connectivite des cellules (dim = n_elt+1)
 * \param [in] cellfac         Table de connectivite des elements (dim = cellfac_idx[n_elt])
 * \param [in] numabs          Numerotation absolue des elements
 *
 */
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
  // PDM_writer_geom_bloc_poly3d_set (*id_cs,
  //                          *id_geom,
  //                          *id_bloc,
  //                          *id_part,
  //                          *n_elt,
  //                          *n_face,
  //                          facsom_idx,
  //                          facsom,
  //                          cellfac_idx,
  //                          cellfac,
  //                          numabs);
}

void
PDM_writer_geom_bloc_poly3d_set
(
PDM_writer_t        *cs,
const int            id_geom,
const int            id_bloc,
const int            id_part,
const PDM_l_num_t    n_elt,
const PDM_l_num_t    n_face,
      PDM_l_num_t   *facsom_idx,
      PDM_l_num_t   *facsom,
      PDM_l_num_t   *cellfac_idx,
      PDM_l_num_t   *cellfac,
      PDM_g_num_t   *numabs
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_Mesh_nodal_block_poly3d_set (geom->mesh_nodal,
                                   id_bloc,
                                   id_part,
                                   n_elt,
                                   n_face,
                                   facsom_idx,
                                   facsom,
                                   cellfac_idx,
                                   cellfac,
                                   numabs,
                                   NULL);
}

/**
 *
 * \brief Ajout de cellules 3D decrites en fonctions des faces.
 *
 * Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * \param [in]  cs              Pointer to \ref PDM_writer object
 * \param [in]  id_geom         Identificateur de l'objet geometrique
 * \param [in]  id_part         Identificateur de partition
 * \param [in]  n_cell          Nombre de cellules 3D ajoutées
 * \param [in]  n_face          Nombre de faces décrites
 * \param [in]  face_som_idx    Index de connectivite faces -> sommets
 * \param [in]  face_som        Connectivite faces -> sommets
 * \param [in]  cell_face_idx   Index de connectivite cellules -> faces
 * \param [in]  cell_face       Connectivite cellules -> faces
 * \param [in]  numabs          Numerotation absolue des cellules
 *
 */
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
  // PDM_writer_geom_cell3d_cellface_add(*id_cs,
  //                                     *id_geom,
  //                                     *id_part,
  //                                     *n_cell,
  //                                     *n_face,
  //                                     face_som_idx,
  //                                     face_som_nb,
  //                                     face_som,
  //                                     cell_face_idx,
  //                                     cell_face_nb,
  //                                     cell_face,
  //                                     numabs);
}

void
PDM_writer_geom_cell3d_cellface_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_cell,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_l_num_t  *cell_face_idx,
 PDM_l_num_t  *cell_face_nb,
 PDM_l_num_t  *cell_face,
 PDM_g_num_t  *numabs
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_Mesh_nodal_cell3d_cellface_add (geom->mesh_nodal,
                                      id_part,
                                      n_cell,
                                      n_face,
                                      face_som_idx,
                                      face_som_nb,
                                      face_som,
                                      cell_face_idx,
                                      cell_face_nb,
                                      cell_face,
                                      numabs);
  if (0 == 1) {
    printf("n_cell : %d\n", n_cell);
    for (int i = 0; i < n_cell; i++) {
      printf ("%d "PDM_FMT_G_NUM" : \n", i+1, numabs[i]);
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        printf (" %d", cell_face[j]);
      }
      printf ("\n");
    }
    printf("n_face : %d\n", n_face);
    for (int i = 0; i < n_face; i++) {
      printf ("%d: \n", i+1);
      for (int j = face_som_idx[i]; j < face_som_idx[i+1]; j++) {
        printf (" %d", face_som[j]);
      }
      printf ("\n");
    }
  }


}


/**
 *
 * \brief Ajout de cellules 2D decrites en fonctions des faces.
 *
 * Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] n_cell          Nombre de cellules 3D ajoutées
 * \param [in] n_face          Nombre de faces décrites
 * \param [in] face_som_idx    Index de connectivite faces -> sommets
 * \param [in] face_som        Connectivite faces -> sommets
 * \param [in] cell_face_idx   Index de connectivite cellules -> faces
 * \param [in] cell_face       Connectivite cellules -> faces
 * \param [in] numabs          Numerotation absolue des cellules
 *
 */
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
  // PDM_writer_geom_cell2d_cellface_add(*id_cs,
  //                                     *id_geom,
  //                                     *id_part,
  //                                     *n_cell,
  //                                     *n_face,
  //                                     face_som_idx,
  //                                     face_som_nb,
  //                                     face_som,
  //                                     cell_face_idx,
  //                                     cell_face_nb,
  //                                     cell_face,
  //                                     numabs);
}

void
PDM_writer_geom_cell2d_cellface_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_cell,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_l_num_t  *cell_face_idx,
 PDM_l_num_t  *cell_face_nb,
 PDM_l_num_t  *cell_face,
 PDM_g_num_t  *numabs
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_Mesh_nodal_cell2d_celledge_add (geom->mesh_nodal,
                                      id_part,
                                      n_cell,
                                      n_face,
                                      face_som_idx,
                                      face_som_nb,
                                      face_som,
                                      cell_face_idx,
                                      cell_face_nb,
                                      cell_face,
                                      numabs);
}


/**
 *
 * \brief Ajout de faces decrites en fonctions des sommets.
 *
 * Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] n_elt           Nombre de cellules 3D ajoutées
 * \param [in] n_face          Nombre de faces décrites
 * \param [in] face_som_idx    Index de connectivite faces -> sommets
 * \param [in] face_som        Connectivite faces -> sommets
 * \param [in] numabs          Numerotation absolue des faces
 *
 */
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
  // PDM_writer_geom_faces_facesom_add(*id_cs,
  //                           *id_geom,
  //                           *id_part,
  //                           *n_face,
  //                           face_som_idx,
  //                           face_som_nb,
  //                           face_som,
  //                           numabs);
}

void
PDM_writer_geom_faces_facesom_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_g_num_t  *numabs
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_Mesh_nodal_faces_facevtx_add (geom->mesh_nodal,
                                    id_part,
                                    n_face,
                                    face_som_idx,
                                    face_som_nb,
                                    face_som,
                                    numabs);
}

/**
 * \brief Ecriture du maillage courant
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */
void
PROCF (pdm_writer_geom_write, PDM_WRITER_GEOM_WRITE)
(
int           *id_cs,
int           *id_geom
)
{
  // PDM_writer_geom_write (*id_cs,
  //              *id_geom);
}

void
PDM_writer_geom_write
(
 PDM_writer_t *cs,
 const int     id_geom
 )
{

    /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  //TODO  faire un retour si geometrie n'est pas dependante du temps
  //       et si on n'est pas au premier increment
  /* Mise a jour du nombre total d'elements */


  /* Determination de la numerotation absolue interne des elements
     Independante du parallelisme */

  const int n_blocks = PDM_Mesh_nodal_n_blocks_get (geom->mesh_nodal);
  const int *blocks_id = PDM_Mesh_nodal_blocks_id_get (geom->mesh_nodal);

  for (int i = 0; i < n_blocks; i++) {
    PDM_Mesh_nodal_g_num_in_block_compute (geom->mesh_nodal, blocks_id[i]);
  }

  /* Ecriture au format */

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->geom_write_fct != NULL) {
    (fmt_ptr->geom_write_fct) (geom);
  }

}


/**
 * \brief Liberation des donnees decrivant le maillage courant
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */
void
PROCF (pdm_writer_geom_free, PDM_WRITER_GEOM_FREE)
(
int           *id_cs,
int           *id_geom
)
{
  // PDM_writer_geom_free(*id_cs,
  //                      *id_geom);

}

void
PDM_writer_geom_free
(
 PDM_writer_t *cs,
 const int     id_geom
)
{
  PDM_writer_geom_data_free(cs,
                            id_geom);

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (cs->geom_tab == NULL) {
    return;
  }

  if (cs->geom_tab->geom == NULL) {
    return;
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom != NULL) {

    PDM_Mesh_nodal_free (geom->mesh_nodal);
    if (geom->nom_geom != NULL) {
      free(geom->nom_geom);
      geom->nom_geom = NULL;
    }

    /* Liberation specifique au format */

    /* Appel de la fonction complementaire propre au format */

    PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

    if (fmt_ptr->geom_free_fct != NULL) {
      (fmt_ptr->geom_free_fct) (geom);
    }

    free(geom);
    geom = NULL;
  }
}


/**
 * \brief Liberation des donnees decrivant le maillage courant
 * les indirections sur les numérotation absolues sont conservées
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */
void
PROCF (pdm_writer_geom_data_free, PDM_WRITER_GEOM_DATA_FREE)
(
int           *id_cs,
int           *id_geom
)
{
  // PDM_writer_geom_data_free(*id_cs,
  //                  *id_geom);
}

void
PDM_writer_geom_data_free
(
 PDM_writer_t *cs,
 const int     id_geom
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (cs->geom_tab == NULL) {
    return;
  }

  if (cs->geom_tab->geom == NULL) {
    return;
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom != NULL) {

    PDM_Mesh_nodal_partial_free (geom->mesh_nodal);

  }
}

/**
 * \brief Mapping des noms de variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] public_name     <-- Nom Public de la variable
 * \param [in] pivate_name     <-- Nom privé de la variable
 *
 * \return   Identificateur de l'objet variable
 *
 */
void
PROCF (pdm_writer_name_map_add_cf, PDM_WRITER_NAME_MAP_ADD_CF)
(
int         *id_cs,
char        *public_name,
int         *l_public_name,
char        *private_name,
int         *l_private_name
ARGF_SUPP_CHAINE
)
{
  char *private_name_c = PDM_fortran_to_c_string(private_name, *l_private_name);
  char *public_name_c = PDM_fortran_to_c_string(public_name, *l_public_name);

  // PDM_writer_name_map_add (*id_cs,
  //                  public_name_c,
  //                  private_name_c);

  if (private_name_c != NULL) {
    free (private_name_c);
  }

  if (public_name_c != NULL) {
    free (public_name_c);
  }
}

void
PDM_writer_name_map_add
(
 PDM_writer_t *cs,
 const char   *public_name,
 const char   *private_name
)
{
  /* Recherche de l'objet cs courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->name_map == NULL) {
    cs->name_map = PDM_Handles_create(3);
  }

  PDM_writer_name_map_t *name_map = (PDM_writer_name_map_t *) malloc (sizeof(PDM_writer_name_map_t));

  PDM_Handles_store (cs->name_map, name_map);

  name_map->public_name = malloc ((strlen(public_name) + 1) * sizeof(char));
  name_map->private_name = malloc ((strlen(private_name) + 1) * sizeof(char));

  strcpy(name_map->public_name, public_name);
  strcpy(name_map->private_name, private_name);

}

/**
 * \brief Creation d'une variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] st_dep_temps    Indique si la variable est dependante du temps
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] dim             Dimension de la variable
 * \param [in] loc             Localisation de la variable
 * \param [in] nom_var         Nom de la variable
 *
 * \return  Identificateur de l'objet variable
 *
 */
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

  // *id_var = PDM_writer_var_create (*id_cs,
  //                                  (PDM_writer_statut_t) *st_dep_tps,
  //                                  (PDM_writer_var_dim_t) *dim,
  //                                  (PDM_writer_var_loc_t) *loc,
  //                                  nom_var_c);

  if (nom_var_c != NULL) {
    free(nom_var_c);
  }
}

int
PDM_writer_var_create
(
 PDM_writer_t               *cs,
 const PDM_writer_statut_t   st_dep_tps,
 const PDM_writer_var_dim_t  dim,
 const PDM_writer_var_loc_t  loc,
 const char                 *nom_var
)
{
  /* Recherche de l'objet cs courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->var_tab == NULL) {
    cs->var_tab = _pdm_writer_var_tab_create(4);
  }

  /* Allocation de la structure PDM_writer_var_t */

  PDM_writer_var_t *var = (PDM_writer_var_t *) malloc(sizeof(PDM_writer_var_t));
  int id_var = _pdm_writer_var_tab_add(cs->var_tab, var);

  /* Initialisation de la structure PDM_writer_var_t */

  _var_init (var);

  size_t l_nom_var = strlen(nom_var);
  var->nom_var = (char *) malloc(sizeof(char) * (l_nom_var + 1));
  strcpy(var->nom_var, nom_var);   /* Nom de la variable */

  var->st_dep_tps = st_dep_tps;    /* Variable en temps */
  var->dim        = dim;           /* Dimension de la variable */
  var->loc        = loc;           /* Dimension de la variable */
  var->_cs        = cs;
  var->private_name = NULL;

  if (cs->name_map != NULL) {
    const int n_map = PDM_Handles_n_get (cs->name_map);
    const int *ind = PDM_Handles_idx_get (cs->name_map);

    for (int i = 0; i < n_map; i++) {
      PDM_writer_name_map_t *map = (PDM_writer_name_map_t *)
              PDM_Handles_get (cs->name_map, ind[i]);
      if (!strcmp(nom_var, map->public_name)) {
        var->private_name = map->private_name;
      }
    }
  }

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->var_create_fct != NULL) {
    (fmt_ptr->var_create_fct) (var);
  }

  return id_var;
}


/**
 * \brief Ecriture des valeurs de la variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_var          Identificateur de la variable a ecrire
 *
 */
void
PROCF (pdm_writer_var_write, PDM_WRITER_VAR_WRITE)
(
int         *id_cs,
int         *id_var
)
{
  // PDM_writer_var_write(*id_cs,
  //            *id_var);
}

void
PDM_writer_var_write
(
 PDM_writer_t *cs,
 const int     id_var
)
{

  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_var >= cs->var_tab->n_var) {
    PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
    abort();
  }

  PDM_writer_var_t *var = cs->var_tab->var[id_var];

  if (var == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad var identifier\n");
  }

  /* Ecriture au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->var_write_fct != NULL) {
    (fmt_ptr->var_write_fct) (var);
  }

}



/**
 * \brief Mise a jour des valeurs de la variable.
 *
 * Attention, les valeurs définies aux elements doivent être définies suivant l'ordre de définition des blocs !
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] id_part         Identificateur de la partition dans l'objet geometrique
 * \param [in] id_var          Identificateur de la variable mise à jour
 * \param [in] val             Valeurs
 *
 */
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
  // PDM_writer_var_set(*id_cs,
  //            *id_var,
  //            *id_geom,
  //            *id_part,
  //            val);
}

void
PDM_writer_var_set
(
 PDM_writer_t     *cs,
 const int         id_var,
 const int         id_geom,
 const int         id_part,
 const PDM_real_t *val
)
{

  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_var >= cs->var_tab->n_var) {
    PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
    abort();
  }

  PDM_writer_var_t *var = cs->var_tab->var[id_var];

  if (var == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad var identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  const int n_ind = cs->geom_tab->n_geom;

  if (var->_val == NULL) {
    var->_val = (double ***) malloc(sizeof(double **) * n_ind);
    for (int i = 0; i < n_ind; i++) {
      var->_val[i] = NULL;
    }
  }

  if (n_ind <= id_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_var_set    : Indice de geometrie incorrect\n");
    abort();
  }

  int n_part = PDM_Mesh_nodal_n_part_get (geom->mesh_nodal);

  if (var->_val[id_geom] == NULL) {
    var->_val[id_geom] = (double **) malloc(sizeof(double *) * n_part);
    for (int i = 0; i < n_part; i++)
      var->_val[id_geom][i] = NULL;
  }

  double **val_geom = var->_val[id_geom];

  if (n_part <= id_part) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_var_set    : Indice de partition incorrect\n");
    abort();
  }

  int n_cell = PDM_Mesh_nodal_n_cell_get (geom->mesh_nodal, id_part);
  int *num_cell_parent_to_local = PDM_Mesh_nodal_num_cell_parent_to_local_get (geom->mesh_nodal, id_part);
  int n_som = PDM_Mesh_nodal_n_vertices_get(geom->mesh_nodal, id_part);

  if (var->loc == PDM_WRITER_VAR_ELEMENTS) {
    val_geom[id_part] = (double *) malloc(sizeof(double) * var->dim * n_cell);
    if (num_cell_parent_to_local != NULL) {
      for (int i = 0; i < n_cell; i++) {
        for (int j = 0; j < (int) var->dim; j++)
          val_geom[id_part][var->dim * num_cell_parent_to_local[i]+j] = val[i*var->dim + j];
      }
    }
    else {
      for (int i = 0; i < n_cell; i++) {
        for (int j = 0; j < (int) var->dim; j++)
          val_geom[id_part][var->dim * i+j] = val[i*var->dim + j];
      }
    }
  }
  else {
    val_geom[id_part] = (double *) malloc(sizeof(double) * var->dim * n_som);
    for (int i = 0; i < n_som * (int) var->dim; i++) {
      val_geom[id_part][i] = val[i];
    }
  }

}


/**
 * \brief Liberation du tableau de donnees des variables
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_var          Identificateur de la variable
 *
 */
void
PROCF (pdm_writer_var_data_free, PDM_WRITER_VAR_DATA_FREE)
(
int         *id_cs,
int         *id_var
)
{
  // PDM_writer_var_data_free(*id_cs,
  //                 *id_var);
}

void
PDM_writer_var_data_free
(
 PDM_writer_t *cs,
 const int     id_var
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_var >= cs->var_tab->n_var) {
    PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
    abort();
  }

  PDM_writer_var_t *var = cs->var_tab->var[id_var];

  if (var != NULL) {

    if (var->_val != NULL) {

      const int n_ind = cs->geom_tab->n_geom;

      for (int i = 0; i < n_ind; i++) {
        int idx = i;
        PDM_writer_geom_t *geom = cs->geom_tab->geom[idx];

        if (geom == NULL) {
          PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
          abort();
        }

        int n_part = PDM_Mesh_nodal_n_part_get (geom->mesh_nodal);

        if ((geom != NULL) && (var->_val[idx] != NULL)) {
          for (int j = 0; j < n_part; j++) {
            if (var->_val[idx][j] != NULL)
              free(var->_val[idx][j]);
            var->_val[idx][j] = NULL;
          }
        }
      }
    }
  }
}


/**
 * \brief Liberation d'une variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_var          Identificateur de la variable
 *
 */
void
PROCF (pdm_writer_var_free, PDM_WRITER_VAR_FREE)
(
int         *id_cs,
int         *id_var
)
{
  // PDM_writer_var_free (*id_cs,
  //             *id_var);
}

void
PDM_writer_var_free
(
 PDM_writer_t *cs,
 const int     id_var
 )
{

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (cs->var_tab != NULL) {

    PDM_writer_var_data_free(cs, id_var);

    /* Acces a l'objet de geometrie courant */

    if (id_var >= cs->var_tab->n_var) {
      PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
      abort();
    }

    PDM_writer_var_t *var = cs->var_tab->var[id_var];

    if (var != NULL) {

      free(var->nom_var);

      const int n_ind = cs->geom_tab->n_geom;

      if (var->_val != NULL) {
        for (int i = 0; i < n_ind; i++) {
          int idx = i;
          if (var->_val[idx] != NULL)
            free(var->_val[idx]);
          var->_val[idx] = NULL;
        }

        free(var->_val);
        var->_val = NULL;
      }

      /* Liberation specifique au format */

      PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

      if (fmt_ptr->var_free_fct != NULL) {
        (fmt_ptr->var_free_fct) (var);
      }

      free (var);
      var = NULL;
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
 const char                  *name,            /*!< Name                                                     */
 const PDM_writer_fct_t       create_fct,      /*!< Customize \ref PDM_writer_create function for the format */
 const PDM_writer_fct_t       free_fct,        /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t       beg_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t       end_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t  geom_create_fct, /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t  geom_write_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t  geom_free_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t   var_create_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t   var_write_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t   var_free_fct     /*!< Customize \ref PDM_writer_free function for the format   */
)
{
  _load_intern_fmt();

  if (geom_write_fct == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_fmt_add : Undefined geom write function\n");
    abort ();
  }

  if (var_write_fct == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_fmt_add : Undefined var write function\n");
    abort ();
  }

  if (n_fmt_tab >= s_fmt_tab) {
    s_fmt_tab = PDM_MAX (2*s_fmt_tab, n_fmt_tab+1);
    fmt_tab = realloc(fmt_tab, sizeof(PDM_writer_fmt_t *) * s_fmt_tab);
  }

  PDM_writer_fmt_t *fmt_ptr = malloc (sizeof(PDM_writer_fmt_t));
  fmt_tab[n_fmt_tab++] = fmt_ptr;

  fmt_ptr->name            = malloc(sizeof(char) * (strlen(name) + 1));
  strcpy (fmt_ptr->name, name);

  fmt_ptr->create_fct      = create_fct;
  fmt_ptr->free_fct        = free_fct;
  fmt_ptr->beg_step_fct    = beg_step_fct;
  fmt_ptr->end_step_fct    = end_step_fct;
  fmt_ptr->geom_create_fct = geom_create_fct;
  fmt_ptr->geom_write_fct  = geom_write_fct;
  fmt_ptr->geom_free_fct   = geom_free_fct;
  fmt_ptr->var_create_fct  = var_create_fct;
  fmt_ptr->var_write_fct   = var_write_fct;
  fmt_ptr->var_free_fct    = var_free_fct;

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

    for (int i = 0; i < n_fmt_tab; i++) {
      if (fmt_tab[i] != NULL) {
        free(fmt_tab[i]->name);
        free(fmt_tab[i]);
        fmt_tab[i] = NULL;
      }
    }

    n_fmt_tab = 0;

    free(fmt_tab);
    fmt_tab = NULL;

  }
}

/**
 * \brief Réinitialisation des donnees decrivant le maillage courant
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */
void
PROCF (pdm_writer_geom_data_reset, PDM_WRITER_GEOM_DATA_RESET)
(
int           *id_cs,
int           *id_geom
)
{
  // PDM_writer_geom_data_reset(*id_cs,
  //                  *id_geom);
}

void
PDM_writer_geom_data_reset
(
 PDM_writer_t *cs,
 const int     id_geom
)
{
  /* Acces a l'objet de geometrie courant */

  // PDM_writer_t *cs = (PDM_writer_t *) PDM_Handles_get (cs_tab, id_cs);
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom != NULL) {

    PDM_Mesh_nodal_reset (geom->mesh_nodal);

  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
