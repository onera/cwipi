/*============================================================================
 * Definition of entity groups
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2007  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distibuted in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_config_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_group.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
#ifdef FVM_CPPCALLER
using namespace bft;
namespace fvm {
#else
extern "C" {
#endif
#if 0
}} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Opaque structure describing a group class
 *----------------------------------------------------------------------------*/

struct _fvm_group_class_t {
  int              n_groups;            /* Number of groups in class */
  int              n_attributes;        /* Number of attributes */
  char           **group_name;          /* Array of group names */
  int             *attribute;           /* Array of attributes */
};

/*----------------------------------------------------------------------------
 * Opaque structure defining a set of group classes
 *----------------------------------------------------------------------------*/

struct _fvm_group_class_set_t {

  int size;                             /* Number of group classes */

  fvm_group_class_t   *a_class;           /* Array of group classes */
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Dump printout of a group class
 *
 * parameters:
 *   group_class         <-- pointer to group class to be dumped
 *   id                  <-- index of group class in set
 *----------------------------------------------------------------------------*/

static void
_group_class_dump(const fvm_group_class_t  *this_group_class,
                  int                       id)
{
  int i;

  if (this_group_class == NULL) {
    bft_printf("\n    _group_class[%d]: nil\n", id);
    return;
  }

  bft_printf("\n"
             "    _group_class[%3d]: %p\n"
             "    n_groups:          %d\n"
             "    n_attributes:      %d\n",
             id, this_group_class,
             this_group_class->n_groups,
             this_group_class->n_attributes);

  if (this_group_class->n_groups > 0) {
    bft_printf("    group names:\n");
    for (i = 0; i < this_group_class->n_groups; i++)
      bft_printf("     \" %s\"\n", this_group_class->group_name[i]);
  }

  if (this_group_class->n_attributes > 0) {
    bft_printf("    attributes:");
    for (i = 0; i < this_group_class->n_attributes; i++)
      bft_printf(" %d", this_group_class->attribute[i]);
    bft_printf("\n");
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the number of groups of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   number of groups in group class
 *----------------------------------------------------------------------------*/

int
fvm_group_class_get_n_groups(const fvm_group_class_t  *this_group_class)
{
  int retval = 0;

  if (this_group_class != NULL) {
    retval = this_group_class->n_groups;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the array of group names of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   pointer to array of group names in group class
 *----------------------------------------------------------------------------*/

const char **
fvm_group_class_get_group_names(const fvm_group_class_t  *this_group_class)
{
  const char **retval = NULL;

  if (this_group_class != NULL) {
    retval = (const char **)(this_group_class->group_name);
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the number of attributes of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   number of attributes in group class
 *----------------------------------------------------------------------------*/

int
fvm_group_class_get_n_attributes(const fvm_group_class_t  *this_group_class)
{
  int  retval = 0;

  if (this_group_class != NULL) {
    retval = this_group_class->n_attributes;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the array of attributes of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   pointer to array of attributes in group class
 *----------------------------------------------------------------------------*/

const int *
fvm_group_class_get_attributes(const fvm_group_class_t  *this_group_class)
{
  const int  *retval = NULL;

  if (this_group_class != NULL) {
    retval = this_group_class->attribute;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Creation of a group class set structure.
 *
 * returns:
 *   pointer to created group class set structure
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t *
fvm_group_class_set_create(void)
{
  fvm_group_class_set_t *class_set;

  BFT_MALLOC(class_set, 1, fvm_group_class_set_t);

  class_set->size = 0;

  class_set->a_class = NULL;

  return class_set;
}

/*----------------------------------------------------------------------------
 * Add a group class to a set.
 *
 * parameters:
 *   this_group_class_set <-> pointer to group class set structure
 *   n_groups             <-- number of groups in class
 *   group_names          <-- array of group names
 *   n_attributes         <-- number of attributes in class
 *   attributes           <-- list of attributes
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_add(fvm_group_class_set_t   *this_group_class_set,
                        int                      n_groups,
                        int                      n_attributes,
                        const char             **group_names,
                        const int               *attributes)
{
  int i;
  fvm_group_class_set_t *class_set = this_group_class_set;
  fvm_group_class_t *_class = NULL;

  assert(class_set != NULL);

  /* Resize array of group class descriptors */

  BFT_REALLOC(class_set->a_class, class_set->size + 1, fvm_group_class_t);

  /* Initialize new descriptor */

  _class = class_set->a_class + class_set->size;

  _class->n_groups = n_groups;
  _class->n_attributes = n_attributes;

  BFT_MALLOC(_class->group_name, n_groups, char *);
  BFT_MALLOC(_class->attribute, n_attributes, int);

  /* Copy group names */

  for (i = 0; i < n_groups; i++) {
    BFT_MALLOC(_class->group_name[i], strlen(group_names[i]) + 1, char);
    strcpy(_class->group_name[i], group_names[i]);
  }

  /* Copy attributes list */

  for (i = 0; i < n_attributes; i++)
    _class->attribute[i] = attributes[i];

  /* Update the number of classes in set */

  class_set->size += 1;
}

/*----------------------------------------------------------------------------
 * Destruction of a group class set structure.
 *
 * parameters:
 *   this_class_set <-- pointer to structure which should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t *
fvm_group_class_set_destroy(fvm_group_class_set_t  *this_group_class_set)
{
  int i, j;

  for (i = 0; i < this_group_class_set->size; i++) {

    fvm_group_class_t *_class = this_group_class_set->a_class + i;

    for (j = 0; j < _class->n_groups; j++)
      BFT_FREE(_class->group_name[j]);

    _class->n_groups = 0;
    _class->n_attributes = 0;

    BFT_FREE(_class->group_name);
    BFT_FREE(_class->attribute);

  }

  BFT_FREE(this_group_class_set->a_class);
  BFT_FREE(this_group_class_set);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Return number of classes in a group class set.
 *
 * parameters:
 *   this_group_class_set <-> pointer to group class set structure
 *
 * returns:
 *   number of classes in a group class set
 *----------------------------------------------------------------------------*/

int
fvm_group_class_set_size(const fvm_group_class_set_t  *this_group_class_set)
{
  int retval = 0;

  if (this_group_class_set != NULL)
    retval = this_group_class_set->size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return pointer to a given group class in a group class set.
 *
 * parameters:
 *   this_group_class_set <-- pointer to group class set structure
 *   group_class_id       <-- index of group class in set (0 to n-1)
 *
 * returns:
 *   pointer to group class structure
 *----------------------------------------------------------------------------*/

const fvm_group_class_t *
fvm_group_class_set_get(const fvm_group_class_set_t  *this_group_class_set,
                        int                           group_class_id)
{
  const fvm_group_class_t  *retval = NULL;

  if (this_group_class_set != NULL) {
    if (group_class_id > -1 && group_class_id < this_group_class_set->size)
      retval = this_group_class_set->a_class + group_class_id;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Dump printout of a group class set
 *
 * parameters:
 *   this_class_set      <-- pointer to group class set to be dumped
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_dump(const fvm_group_class_set_t  *this_group_class_set)
{
  int i;
  const fvm_group_class_set_t  *class_set = this_group_class_set;

  if (this_group_class_set == NULL) {
    bft_printf("  group_class_set: nil\n");
    return;
  }

  bft_printf("  _group_class_set: %p\n"
             "  size:             %d\n",
             class_set,
             class_set->size);

  if (class_set->size > 0) {
    bft_printf("\n  group_classes:");
    for (i = 0; i < class_set->size; i++)
      _group_class_dump(class_set->a_class + i, i);
  }

  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
