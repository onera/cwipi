
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal_elmts.h"
#include "pdm_dmesh_nodal_elmts_priv.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/
/**
 *
 * \brief Free a standard section
 *
 * \param [inout]  _bloc_std    Standard section
 *
 * \return         Null
 *
 */

static
void
_section_std_free
(
PDM_DMesh_nodal_section_std_t *_section_std
)
{
  if (_section_std == NULL) {
    return;
  }

  if (_section_std->distrib != NULL) {
    free (_section_std->distrib);
    _section_std->distrib = NULL;
  }

  if(_section_std->owner == PDM_OWNERSHIP_KEEP) {
    if(_section_std->_connec != NULL) {
      free(_section_std->_connec);
      _section_std->_connec = NULL;
    }
  }

  free(_section_std);
}


/**
 *
 * \brief Free a polygon section
 *
 * \param [inout]  _bloc_poly2d    Polygon section
 *
 * \return         Null
 *
 */

static
void
_section_poly2d_free
(
PDM_DMesh_nodal_section_poly2d_t *_section_poly2d
)
{

  if (_section_poly2d == NULL) {
    return;
  }

  if (_section_poly2d->distrib != NULL) {
    free (_section_poly2d->distrib);
    _section_poly2d->distrib = NULL;
  }

  if(_section_poly2d->owner == PDM_OWNERSHIP_KEEP) {
    if(_section_poly2d->_connec != NULL) {
      free(_section_poly2d->_connec);
      _section_poly2d->_connec = NULL;
    }
    if(_section_poly2d->_connec_idx != NULL) {
      free(_section_poly2d->_connec_idx);
      _section_poly2d->_connec_idx = NULL;
    }
  }

  free(_section_poly2d);
}

/**
 *
 * \brief Free a polyhedron section
 *
 * \param [inout]  _section_poly3d    Polyhedron section
 *
 * \return         Null
 *
 */
static
void
_section_poly3d_free
(
PDM_DMesh_nodal_section_poly3d_t *_section_poly3d
)
{
  if (_section_poly3d == NULL) {
    return;
  }

  if (_section_poly3d->distrib != NULL) {
    free (_section_poly3d->distrib);
    _section_poly3d->distrib = NULL;
  }

  if(_section_poly3d->owner == PDM_OWNERSHIP_KEEP) {
    if(_section_poly3d->_face_vtx_idx != NULL) {
      free(_section_poly3d->_face_vtx_idx);
      _section_poly3d->_face_vtx_idx = NULL;
    }
    if(_section_poly3d->_face_vtx != NULL) {
      free(_section_poly3d->_face_vtx);
      _section_poly3d->_face_vtx = NULL;
    }
    if(_section_poly3d->_cell_face_idx != NULL) {
      free(_section_poly3d->_cell_face_idx);
      _section_poly3d->_cell_face_idx = NULL;
    }
    if(_section_poly3d->_cell_face != NULL) {
      free(_section_poly3d->_cell_face);
      _section_poly3d->_cell_face = NULL;
    }
  }
  free(_section_poly3d);
}


/**
 *
 * \brief Free a list of standard section
 *
 * \param [inout]  sections    standard sections
 * \param [inout]  n_sections  Number of standard sections
 */
static
void
_sections_std_free
(
 PDM_DMesh_nodal_section_std_t **sections,
 int                             n_sections
)
{
  for(int i_section = 0; i_section < n_sections; ++i_section) {
    _section_std_free(sections[i_section]);
  }
  if(sections != NULL){
    free(sections);
    sections = NULL;
  }
}

/**
 *
 * \brief Free a list of polygon section
 *
 * \param [inout]  sections     Polygon sections
 * \param [inout]  n_sections   Number of polygon sections
 */
static
void
_sections_poly2d_free
(
 PDM_DMesh_nodal_section_poly2d_t **sections,
 int                                n_sections
)
{
  for(int i_section = 0; i_section < n_sections; ++i_section) {
    _section_poly2d_free(sections[i_section]);
  }
  if(sections != NULL){
    free(sections);
    sections = NULL;
  }
}


/**
 *
 * \brief Free a list of polyhedron section
 *
 * \param [inout]  sections    Polyhedron sections
 * \param [inout]  n_sections  Number of polyhedron sections
 */
static
void
_sections_poly3d_free
(
 PDM_DMesh_nodal_section_poly3d_t **sections,
 int                                n_sections
)
{
  for(int i_section = 0; i_section < n_sections; ++i_section) {
    _section_poly3d_free(sections[i_section]);
  }
  if(sections != NULL){
    free(sections);
    sections = NULL;
  }
}

static void
_update_elmt_sections_id
(
 PDM_dmesh_nodal_elmts_t *dmn_elts
)
{
  int n_section = 0;


  if (dmn_elts->sections_std != NULL) {
    n_section += dmn_elts->n_section_std;
  }

  if (dmn_elts->sections_poly2d != NULL) {
    n_section += dmn_elts->n_section_poly2d;
  }

  if (dmn_elts->sections_poly3d != NULL) {
    n_section += dmn_elts->n_section_poly3d;
  }

  if (dmn_elts->n_section < n_section) {
    dmn_elts->sections_id = (int *) realloc(dmn_elts->sections_id, sizeof(int) * n_section);
  }

  int k = 0;
  if (dmn_elts->sections_std != NULL) {
    for (int i = 0; i < dmn_elts->n_section_std; i++) {
      dmn_elts->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_STD;
    }
  }

  if (dmn_elts->sections_poly2d != NULL) {
    for (int i = 0; i < dmn_elts->n_section_poly2d; i++) {
      dmn_elts->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }

  if (dmn_elts->sections_poly3d != NULL) {
    for (int i = 0; i < dmn_elts->n_section_poly3d; i++) {
      dmn_elts->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }

  dmn_elts->n_section = n_section;

  // printf("_update_elmt_sections_id | n_section                   = %i \n", n_section);
  // printf("_update_elmt_sections_id | dmn_elts->n_section_std     = %i \n", dmn_elts->n_section_std);
  // printf("_update_elmt_sections_id | dmn_elts->n_section_poly2d  = %i \n", dmn_elts->n_section_poly2d);
  // printf("_update_elmt_sections_id | dmn_elts->n_section_poly3d  = %i \n", dmn_elts->n_section_poly3d);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_nodal_elmts_t*
PDM_DMesh_nodal_elmts_create
(
const PDM_MPI_Comm comm,
      int          mesh_dimension,
      PDM_g_num_t  n_g_elmts
)
{
  PDM_dmesh_nodal_elmts_t *dmn_elts = (PDM_dmesh_nodal_elmts_t *) malloc (sizeof(PDM_dmesh_nodal_elmts_t));

  dmn_elts->comm           = comm;
  dmn_elts->mesh_dimension = mesh_dimension;
  dmn_elts->n_g_elmts      = n_g_elmts;

  PDM_MPI_Comm_size(dmn_elts->comm, &dmn_elts->n_rank);
  PDM_MPI_Comm_rank(dmn_elts->comm, &dmn_elts->i_rank);

  dmn_elts->n_section            = 0;
  dmn_elts->n_section_std        = 0;
  dmn_elts->n_section_poly2d     = 0;
  dmn_elts->n_section_poly3d     = 0;

  dmn_elts->sections_id          = NULL;
  dmn_elts->sections_std         = NULL;
  dmn_elts->sections_poly3d      = NULL;
  dmn_elts->sections_poly2d      = NULL;
  dmn_elts->section_distribution = NULL;

  dmn_elts->n_group_elmt         = 0;
  dmn_elts->dgroup_elmt_idx      = NULL;
  dmn_elts->dgroup_elmt          = NULL;
  dmn_elts->dgroup_elmt_owner    = PDM_OWNERSHIP_KEEP;

  dmn_elts->dparent_idx          = NULL;
  dmn_elts->dparent_gnum         = NULL;
  dmn_elts->dparent_sign         = NULL;
  dmn_elts->delmt_child_distrib  = NULL;

  return dmn_elts;
}



void
PDM_DMesh_nodal_elmts_group_set
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                n_group_elmt,
      int               *dgroup_elmt_idx,
      PDM_g_num_t       *dgroup_elmt,
      PDM_ownership_t    owner
)
{
  dmn_elts->n_group_elmt      = n_group_elmt;
  dmn_elts->dgroup_elmt_idx   = dgroup_elmt_idx;
  dmn_elts->dgroup_elmt       = dgroup_elmt;
  dmn_elts->dgroup_elmt_owner = owner;
}


void
PDM_DMesh_nodal_elmts_free
(
PDM_dmesh_nodal_elmts_t* dmn_elts
)
{
  if (dmn_elts->sections_id != NULL) {
    free (dmn_elts->sections_id);
  }
  dmn_elts->sections_id = NULL;

  if(dmn_elts->section_distribution != NULL) {
    free (dmn_elts->section_distribution);
  }
  dmn_elts->section_distribution = NULL;

  _sections_std_free   (dmn_elts->sections_std   , dmn_elts->n_section_std   );
  _sections_poly3d_free(dmn_elts->sections_poly3d, dmn_elts->n_section_poly3d);
  _sections_poly2d_free(dmn_elts->sections_poly2d, dmn_elts->n_section_poly2d);

  if(dmn_elts->dgroup_elmt_owner == PDM_OWNERSHIP_KEEP) {
    if (dmn_elts->dgroup_elmt != NULL) {
      free (dmn_elts->dgroup_elmt);
    }
    if (dmn_elts->dgroup_elmt_idx != NULL) {
      free (dmn_elts->dgroup_elmt_idx);
    }
  }

  if(dmn_elts->dparent_idx != NULL) {
    free(dmn_elts->dparent_idx);
  }
  if(dmn_elts->dparent_gnum != NULL) {
    free(dmn_elts->dparent_gnum);
  }
  if(dmn_elts->dparent_sign != NULL) {
    free(dmn_elts->dparent_sign);
  }

  if(dmn_elts->delmt_child_distrib != NULL) {
    free(dmn_elts->delmt_child_distrib);
  }

  free(dmn_elts);
}

int
PDM_DMesh_nodal_elmts_section_add
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const PDM_Mesh_nodal_elt_t     t_elt
)
{
  if (dmn_elts == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int id_block = -1;
  dmn_elts->n_section++;

  dmn_elts->sections_id  = realloc(dmn_elts->sections_id , sizeof(int) * dmn_elts->n_section);

  if(t_elt == PDM_MESH_NODAL_POINT) {
    if(dmn_elts->mesh_dimension != 0){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", dmn_elts->mesh_dimension, 0);
    }
  } else if(t_elt == PDM_MESH_NODAL_BAR2) {
    if(dmn_elts->mesh_dimension != 1){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", dmn_elts->mesh_dimension, 1);
    }
  } else if(t_elt == PDM_MESH_NODAL_TRIA3 || t_elt == PDM_MESH_NODAL_QUAD4 || t_elt == PDM_MESH_NODAL_POLY_2D) {
    if(dmn_elts->mesh_dimension != 2){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", dmn_elts->mesh_dimension, 2);
    }
  } else {
    if(dmn_elts->mesh_dimension != 3){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", dmn_elts->mesh_dimension, 3);
    }
  }

  switch (t_elt) {

  case PDM_MESH_NODAL_POINT    :
  case PDM_MESH_NODAL_BAR2     :
  case PDM_MESH_NODAL_TRIA3    :
  case PDM_MESH_NODAL_QUAD4    :
  case PDM_MESH_NODAL_TETRA4   :
  case PDM_MESH_NODAL_PYRAMID5 :
  case PDM_MESH_NODAL_PRISM6   :
  case PDM_MESH_NODAL_HEXA8    :
    {
      id_block = dmn_elts->n_section_std++;

      dmn_elts->sections_std = realloc(dmn_elts->sections_std, dmn_elts->n_section_std * sizeof(PDM_DMesh_nodal_section_std_t * ));
      dmn_elts->sections_std[id_block]          = malloc( sizeof(PDM_DMesh_nodal_section_std_t) );
      dmn_elts->sections_std[id_block]->t_elt   = t_elt;
      dmn_elts->sections_std[id_block]->n_elt   = -1;
      dmn_elts->sections_std[id_block]->_connec = NULL;
      dmn_elts->sections_std[id_block]->distrib = NULL;
      dmn_elts->sections_std[id_block]->owner   = PDM_OWNERSHIP_KEEP;

      id_block += PDM_BLOCK_ID_BLOCK_STD;
    }
    break;

  case PDM_MESH_NODAL_POLY_2D  :
    {
      assert(dmn_elts->mesh_dimension == 2);
      id_block = dmn_elts->n_section_poly2d++;

      dmn_elts->sections_poly2d = realloc(dmn_elts->sections_poly2d, dmn_elts->n_section_poly2d * sizeof(PDM_DMesh_nodal_section_poly2d_t *));
      dmn_elts->sections_poly2d[id_block]              = malloc( sizeof(PDM_DMesh_nodal_section_poly2d_t) );
      dmn_elts->sections_poly2d[id_block]->n_elt       = -1;
      dmn_elts->sections_poly2d[id_block]->_connec     = NULL;
      dmn_elts->sections_poly2d[id_block]->_connec_idx = NULL;
      dmn_elts->sections_poly2d[id_block]->distrib     = NULL;
      dmn_elts->sections_poly2d[id_block]->owner       = PDM_OWNERSHIP_KEEP;

      id_block += PDM_BLOCK_ID_BLOCK_POLY2D;

    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :
    {
      id_block = dmn_elts->n_section_poly3d++;
      assert(dmn_elts->mesh_dimension == 3);

      dmn_elts->sections_poly3d = realloc(dmn_elts->sections_poly3d, dmn_elts->n_section_poly3d * sizeof(PDM_DMesh_nodal_section_poly3d_t *));
      dmn_elts->sections_poly3d[id_block]                 = malloc( sizeof(PDM_DMesh_nodal_section_poly3d_t) );
      dmn_elts->sections_poly3d[id_block]->n_elt          = -1;
      dmn_elts->sections_poly3d[id_block]->n_face         = -1;
      dmn_elts->sections_poly3d[id_block]->_face_vtx_idx  = NULL;
      dmn_elts->sections_poly3d[id_block]->_face_vtx      = NULL;
      dmn_elts->sections_poly3d[id_block]->_cell_face_idx = NULL;
      dmn_elts->sections_poly3d[id_block]->_cell_face     = NULL;
      dmn_elts->sections_poly3d[id_block]->distrib        = NULL;
      dmn_elts->sections_poly3d[id_block]->owner          = PDM_OWNERSHIP_KEEP;

      id_block += PDM_BLOCK_ID_BLOCK_POLY3D;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  _update_elmt_sections_id(dmn_elts);
  return id_block ;
}



/* Update the owner ship of the different composents of dmn_elts */
void
PDM_DMesh_nodal_elmts_update_ownership
(
  PDM_dmesh_nodal_elmts_t *dmn_elts,
  PDM_ownership_t          owner
)
{
  for(int i_section = 0; i_section < dmn_elts->n_section_std; ++i_section) {
    dmn_elts->sections_std[i_section]->owner = owner;
  }
  for(int i_section = 0; i_section < dmn_elts->n_section_poly2d; ++i_section) {
    dmn_elts->sections_poly2d[i_section]->owner = owner;
  }
  for(int i_section = 0; i_section < dmn_elts->n_section_poly3d; ++i_section) {
    dmn_elts->sections_poly3d[i_section]->owner = owner;
  }
  dmn_elts->dgroup_elmt_owner = owner;
}


/**
 * \brief Define a standard section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 *
 */
void
PDM_DMesh_nodal_elmts_section_std_set
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                id_section,
const int                n_elt,
      PDM_g_num_t       *connec,
      PDM_ownership_t    owner
)
{
  if (dmn_elts == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  // PDM_printf("PDM_DMesh_nodal_elmts_section_std_set - _id_section : %i  \n", _id_section);
  // PDM_printf("PDM_DMesh_nodal_elmts_section_std_set - n_elt       : %i  \n", n_elt);

  /* Mapping */
  section->n_elt   = n_elt;
  section->_connec = connec;
  section->owner   = owner;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmn_elts->n_rank + 1));

  /* Creation of distribution */
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmn_elts->comm);

  section->distrib[0] = 0;
  for (int i = 1; i < dmn_elts->n_rank + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }
}


/**
 * \brief Define a polygon section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_elmts_section_poly2d_set
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section,
const PDM_l_num_t              n_elt,
      PDM_l_num_t             *connec_idx,
      PDM_g_num_t             *connec,
      PDM_ownership_t          owner
)
{
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_DMesh_nodal_section_poly2d_t *section = dmn_elts->sections_poly2d[_id_section];

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  /* Mapping */

  section->n_elt       = n_elt;
  section->_connec_idx = connec_idx;
  section->_connec     = connec;
  section->owner       = owner;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmn_elts->n_rank + 1));

  /* Creation of distribution */
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmn_elts->comm);

  section->distrib[0] = 0;
  for (int i = 1; i < dmn_elts->n_rank + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }
}


/**
 * \brief Return a polygon section description
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_elmts_section_poly2d_get
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
      PDM_l_num_t             **connec_idx,
      PDM_g_num_t             **connec
)
{
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  PDM_DMesh_nodal_section_poly2d_t *section = dmn_elts->sections_poly2d[_id_section];

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *connec_idx = section->_connec_idx;
  *connec     = section->_connec;
}



/**
 * \brief  Return section distribution
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  A array of size \ref n_ranks + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_elmts_distrib_section_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
)
{
  if (dmn_elts == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_DMesh_nodal_section_poly3d_t *section = dmn_elts->sections_poly3d[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }

    return section->distrib;
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_DMesh_nodal_section_poly2d_t *section = dmn_elts->sections_poly2d[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }

    return section->distrib;
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }

    return section->distrib;
  }
}

PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_elmts_section_type_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
)
{
  if (dmn_elts == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_POLY_3D;

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
    }

    t_elt = section->t_elt;
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    t_elt = PDM_MESH_NODAL_POLY_2D;

  }

  else {

    t_elt = PDM_MESH_NODAL_POLY_3D;

  }

  return t_elt;
}


void
PDM_DMesh_nodal_elmts_section_poly3d_set
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
const PDM_l_num_t               n_elt,
const PDM_l_num_t               n_face,
      PDM_l_num_t              *facvtx_idx,
      PDM_g_num_t              *facvtx,
      PDM_l_num_t              *cellfac_idx,
      PDM_g_num_t              *cellfac,
      PDM_ownership_t           owner
)
{
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

  PDM_DMesh_nodal_section_poly3d_t *section = dmn_elts->sections_poly3d[_id_section];

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  section->n_elt          = n_elt;
  section->n_face         = n_face;
  section->_face_vtx_idx  = facvtx_idx;
  section->_face_vtx      = facvtx;
  section->_cell_face_idx = cellfac_idx;
  section->_cell_face     = cellfac;
  section->owner          = owner;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmn_elts->n_rank + 1));

  /* Creation of distribution */
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmn_elts->comm);

  section->distrib[0] = 0;
  for (int i = 1; i < dmn_elts->n_rank + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }
}



void
PDM_DMesh_nodal_elmts_section_poly3d_get
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
      PDM_l_num_t              *n_face,
      PDM_l_num_t             **facvtx_idx,
      PDM_g_num_t             **facvtx,
      PDM_l_num_t             **cellfac_idx,
      PDM_g_num_t             **cellfac
)
{
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

  PDM_DMesh_nodal_section_poly3d_t *section = dmn_elts->sections_poly3d[_id_section];

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *n_face      = section->n_face;
  *facvtx_idx  = section->_face_vtx_idx;
  *facvtx      = section->_face_vtx;
  *cellfac_idx = section->_cell_face_idx;
  *cellfac     = section->_cell_face;

}


PDM_g_num_t
PDM_DMesh_nodal_elmts_total_n_elmt_get
(
 PDM_dmesh_nodal_elmts_t  *dmn_elts
)
{
  PDM_g_num_t total_n_elmt = 0;
  for(int i_section = 0; i_section < dmn_elts->n_section_std; ++i_section){
    total_n_elmt += dmn_elts->sections_std[i_section]->distrib[dmn_elts->n_rank];
  }

  for(int i_section = 0; i_section < dmn_elts->n_section_poly2d; ++i_section){
    total_n_elmt += dmn_elts->sections_poly2d[i_section]->distrib[dmn_elts->n_rank];
  }

  for(int i_section = 0; i_section < dmn_elts->n_section_poly3d; ++i_section){
    total_n_elmt += dmn_elts->sections_poly3d[i_section]->distrib[dmn_elts->n_rank];
  }

  return total_n_elmt;
}


/**
*
* \brief PDM_sections_decompose_faces
*
* \param [in]     hdl                Distributed nodal mesh handle
* \param [inout]  n_face_elt_tot     Number of faces
* \param [inout]  n_sum_vtx_face_tot Number of vtx for all faces (cumulative)
*
*/
void
PDM_dmesh_nodal_elmts_decompose_faces_get_size
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
int                     *n_face_elt_tot,
int                     *n_sum_vtx_face_tot
)
{
  /* Get current structure to treat */
  *n_face_elt_tot     = 0;
  *n_sum_vtx_face_tot = 0;

  for (int i_section = 0; i_section < dmn_elts->n_section_std; i_section++) {

    int n_face_elt     = PDM_n_face_elt_per_elmt    (dmn_elts->sections_std[i_section]->t_elt);
    int n_sum_vtx_face = PDM_n_sum_vtx_face_per_elmt(dmn_elts->sections_std[i_section]->t_elt);

    *n_face_elt_tot     += dmn_elts->sections_std[i_section]->n_elt * n_face_elt;
    *n_sum_vtx_face_tot += dmn_elts->sections_std[i_section]->n_elt * n_sum_vtx_face;

  }

  for (int i = 0; i < dmn_elts->n_section_poly3d; i++) {
    int _n_face = dmn_elts->sections_poly3d[i]->n_face;
    *n_face_elt_tot     += _n_face;
    *n_sum_vtx_face_tot += dmn_elts->sections_poly3d[i]->_face_vtx[dmn_elts->sections_poly3d[i]->_face_vtx_idx[_n_face]];
  }

  for (int i = 0; i < dmn_elts->n_section_poly2d; i++) {
    *n_face_elt_tot     += dmn_elts->sections_poly2d[i]->n_elt;
    *n_sum_vtx_face_tot += dmn_elts->sections_poly2d[i]->_connec_idx[dmn_elts->sections_poly2d[i]->n_elt];
  }

  assert(dmn_elts->n_section_poly3d == 0); // Not implemented
  assert(dmn_elts->n_section_poly2d == 0); // Not implemented

  // printf("n_face_elt_tot     ::%i\n", *n_face_elt_tot   );
  // printf("n_sum_vtx_face_tot::%i\n" , *n_sum_vtx_face_tot);
}


/**
*
* \brief PDM_dmesh_nodal_decompose_edges_get_size
*
* \param [in]     hdl                Distributed nodal mesh handle
* \param [inout]  n_edge_elt_tot     Number of edges
* \param [inout]  n_sum_vtx_edge_tot Number of vtx for all edges (cumulative)
*
*/
void
PDM_dmesh_nodal_elmts_decompose_edges_get_size
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
int                     *n_edge_elt_tot,
int                     *n_sum_vtx_edge_tot
)
{
  *n_edge_elt_tot     = 0;
  *n_sum_vtx_edge_tot = 0;

  if(dmn_elts == NULL) {
    return;
  }

  for (int i_section = 0; i_section < dmn_elts->n_section_std; i_section++) {

    int n_edge_elt     = PDM_n_nedge_elt_per_elmt   (dmn_elts->sections_std[i_section]->t_elt);
    int n_sum_vtx_edge = PDM_n_sum_vtx_edge_per_elmt(dmn_elts->sections_std[i_section]->t_elt);

    *n_edge_elt_tot     += dmn_elts->sections_std[i_section]->n_elt*n_edge_elt;
    *n_sum_vtx_edge_tot += dmn_elts->sections_std[i_section]->n_elt*n_sum_vtx_edge;

  }

  for (int i = 0; i < dmn_elts->n_section_poly3d; i++) {
    int _n_face = dmn_elts->sections_poly3d[i]->n_face;
    *n_edge_elt_tot     +=     dmn_elts->sections_poly3d[i]->_face_vtx[dmn_elts->sections_poly3d[i]->_face_vtx_idx[_n_face]];
    *n_sum_vtx_edge_tot += 2 * dmn_elts->sections_poly3d[i]->_face_vtx[dmn_elts->sections_poly3d[i]->_face_vtx_idx[_n_face]];
  }

  for (int i = 0; i < dmn_elts->n_section_poly2d; i++) {
    *n_edge_elt_tot     +=     dmn_elts->sections_poly2d[i]->_connec_idx[dmn_elts->sections_poly2d[i]->n_elt];
    *n_sum_vtx_edge_tot += 2 * dmn_elts->sections_poly2d[i]->_connec_idx[dmn_elts->sections_poly2d[i]->n_elt];
  }

}

void
PDM_dmesh_nodal_elmts_generate_distribution
(
 PDM_dmesh_nodal_elmts_t *dmn_elts
)
{
  /* Creation of element distribution among all sections */
  // printf("dmn_elts->n_section : %i \n", dmn_elts->n_section);
  if(dmn_elts->section_distribution == NULL) {
    dmn_elts->section_distribution = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmn_elts->n_section + 1));
  }
  dmn_elts->section_distribution[0] = 0;

  for(int i_section = 0; i_section < dmn_elts->n_section; ++i_section) {

    int id_section = dmn_elts->sections_id[i_section];
    const PDM_g_num_t* distrib = PDM_DMesh_nodal_elmts_distrib_section_get(dmn_elts, id_section);

    dmn_elts->section_distribution[i_section+1] = dmn_elts->section_distribution[i_section] + distrib[dmn_elts->n_rank];
  }

  /* Verbose */
  if(0 == 1)
  {
    PDM_printf(" ------------------------------ \n");
    for(int i_section=0; i_section < dmn_elts->n_section+1; i_section++){
      PDM_printf("%i ", dmn_elts->section_distribution[i_section]);
    }
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
