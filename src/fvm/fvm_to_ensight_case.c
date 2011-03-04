/*============================================================================
 * Manage case files associated with the EnSight Gold writer
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2006  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------  
 *  mpi.h must be include before stdio.h to not define SEE_SET for C++ 
 *  binding of MPI
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
#include "fvm_config.h"
#if defined(FVM_HAVE_MPI)
#include <mpi.h>
#endif
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>  /* toupper() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_file.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_config_defs.h"
#include "fvm_defs.h"
#include "fvm_gather.h"
#include "fvm_parall.h"
#include "fvm_writer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_ensight_case.h"

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
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Time set entry structure
 *----------------------------------------------------------------------------*/

typedef struct {

  int           n_time_values;   /* Number of time step values */
  int           last_time_step;  /* Last (current) time step number */

  double       *time_value;      /* Time step values */

} fvm_to_ensight_case_time_t;

/*----------------------------------------------------------------------------
 * Variable entry structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char         *name;            /* Variable name */
  char         *case_line;       /* Line in case file */
  char         *file_name;       /* Variable file name */

  int           time_set;        /* Associated time set index */
  int           last_time_step;  /* Last (current) time step number */
  int           dim;             /* Associated dimension: 0 (constant), 1
                                    1 (scalar), 3 (vector), 6 (symmetrical
                                    tensor), or 9 (asymmetrical tensor) */
  fvm_writer_var_loc_t  loc;  /* variable at nodes, elements, or particles */

} fvm_to_ensight_case_var_t;

/*----------------------------------------------------------------------------
 * EnSight case file structure
 *----------------------------------------------------------------------------*/

struct _fvm_to_ensight_case_t {

  char          *name;              /* Case name */
  char          *case_file_name;    /* Case file name */

  char          *file_name_prefix;  /* File name prefix */
  int            dir_name_length;   /* Associated directory name length
                                       (may be 0); index in file_name_prefix
                                       corresponding to the base file name */

  char          *geom_file_name;    /* Geometry file name */

  int            n_parts;           /* Number of referenced parts */
  char         **part_name;         /* Part names (used as unique identifier) */

  int                           n_time_sets;  /* Number of time sets */
  fvm_to_ensight_case_time_t  **time_set;     /* Time Set entries */

  int                           n_vars;       /* Number of variables */
  fvm_to_ensight_case_var_t   **var;          /* Variable entries */

  int            geom_time_set;               /* Index of time set
                                                 associated with geometry */

  fvm_writer_time_dep_t   time_dependency;    /* Mesh time dependency */

  _Bool          geom_info_queried;           /* Indicated if current
                                                 geometry file name queried */

  _Bool          modified;                    /* Modified since last output ? */

} _fvm_to_ensight_case_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a new time set structure.
 *
 * returns:
 *   pointer to new time set structure
 *----------------------------------------------------------------------------*/

static fvm_to_ensight_case_time_t *
_time_set_create(void)
{
  fvm_to_ensight_case_time_t   *this_time = NULL;

  /* Create and initialize structure */

  BFT_MALLOC(this_time, 1, fvm_to_ensight_case_time_t);

  this_time->n_time_values = 0;
  this_time->last_time_step = -1;

  this_time->time_value = NULL;

  /* Return new structure */

  return this_time;
}

/*----------------------------------------------------------------------------
 * Create a new time set structure.
 *
 * parameters:
 *   this_case  <-- time set structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

static fvm_to_ensight_case_time_t *
_time_set_destroy(fvm_to_ensight_case_time_t  *this_time)
{
  BFT_FREE(this_time->time_value);

  BFT_FREE(this_time);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Add a new time step number and value to a time set if necessary:
 * if the corresponding physical time is not present in the structure,
 * the corresponding elements are added.
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *   time_step  <-- number of time step to add
 *   time_value <-- associated time value
 *
 * returns:
 *   0 if no time was added, 1 if a new time was added
 *----------------------------------------------------------------------------*/

static int
_add_time(fvm_to_ensight_case_time_t  *const time_set,
          const int                          time_step,
          const double                       time_value)
{
  const char time_value_err_string[] =
    N_("The time value associated with time step <%d> equals <%g>,\n"
       "but time value <%g> has already been associated with this time step.\n");

  int modified = 0;

  /* First, make sure only increasing time steps are defined */

  if (time_step < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("The given time step value should be >= 0, and not %d.\n"),
              time_step);

  else if (time_step < time_set->last_time_step)
    bft_error(__FILE__, __LINE__, 0,
              _("The given time step value should be >= %d, and not %d.\n"),
              time_set->last_time_step, time_step);

  /* Second, check that non conflicting definitions
     for the current time step are given */

  else if (time_step == time_set->last_time_step) {

    double last_time_value = time_set->time_value[time_set->n_time_values - 1];

    if (   time_value < last_time_value - (1.0 + 1.e-16)
        || time_value > last_time_value + (1.0 + 1.e-16))
      bft_error(__FILE__, __LINE__, 0,
                _(time_value_err_string),
                time_step, time_value, last_time_value);

  }

  /* Finally, add a new time step and value if necessary */

  else {

    time_set->last_time_step  = time_step;
    time_set->n_time_values += 1;

    BFT_REALLOC(time_set->time_value, time_set->n_time_values, double);

    time_set->time_value[time_set->n_time_values - 1] = time_value;

    modified = 1;

  }

  return modified;
}

/*----------------------------------------------------------------------------
 * Add a new variable entry
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *   name       <-- variable name
 *   dimension  <-- variable dimension (0: constant, 1: scalar, 3: vector,
 *                  6: symmetrical tensor, 9: asymmetrical tensor)
 *   location   <-- variable definition location (nodes, elements, or particles)
 *   time_set   <-- associated time set index
 *----------------------------------------------------------------------------*/

static void
_add_var(fvm_to_ensight_case_t       *const this_case,
         const char                  *const name,
         const int                          dimension,
         const fvm_writer_var_loc_t         location,
         const int                          time_set)
{
  char line[80], description[20];
  int i, l;
  int remain_len, prefix_len, base_len, postfix_len;
  _Bool rename_file;

  fvm_to_ensight_case_var_t  *var;

  l = strlen(name);

  this_case->n_vars += 1;
  BFT_MALLOC(var, 1, fvm_to_ensight_case_var_t);

  BFT_MALLOC(var->name, l + 1, char);
  strcpy(var->name, name);

  /* Create description (19 chars max) */

  if (l > 19)
      l = 19;

  strncpy(description, name, l);
  description[l] = '\0';

  /* Some characters not allowed in format, replaced by '_' */

  for (i = 0 ; i < l ; i++) {
    switch (description[i]) {
    case '(':
    case ')':
    case ']':
    case '[':
    case '+':
    case '-':
    case '@':
    case ' ':
    case '\t':
    case '!':
    case '#':
    case '*':
    case '^':
    case '$':
    case '/':
      description[i] = '_';
      break;
    default:
      break;
    }
  }

  /* Assign time set to obtain file index, and dimension and location
     before case line creation */

  var->time_set = time_set;
  var->last_time_step = -1;
  var->dim = dimension;
  var->loc = location;

  /* Create associated case file line, up to file name
     (which may depend on the number of remaining characters,
     so as to avoid going beyond 79 characters if possible) */

  switch(var->dim) {
  case 0:
    strcpy(line, "constant per case file: ");
    break;
  case 1:
    strcpy(line, "scalar per ");
    break;
  case 3:
    strcpy(line, "vector per ");
    break;
  case 6:
    strcpy(line, "tensor symm per ");
    break;
  case 9:
    strcpy(line, "tensor asym per ");
    break;
  }

  if (var->dim > 0) {
    switch(var->loc) {
    case FVM_WRITER_PER_NODE:
      strcat(line, "node:    ");
      break;
    case FVM_WRITER_PER_ELEMENT:
      strcat(line, "element: ");
      break;
    case FVM_WRITER_PER_PARTICLE:
      strcat(line, "measured node: ");
      break;
    }
  }

  l = strlen(line); /* At this stage, l = 31 at most */

  if (var->time_set > -1)
    sprintf(line + l, "%d ", var->time_set + 1);
  else
    strcat(line, "  ");

  l = strlen(line);  /* At this stage, l = 35 at most with
                        time set number < 100 (EnSight maximum:
                        16, apparently only for measured data) */

  sprintf(line + l, "%19s ", description); /* Description max 19 chars,
                                              so 79 char line limit not
                                              exceeded here */

  for (l = strlen(line) ; l < 48 ; l++)
    line[l] = ' ';
  line[l] = '\0'; /* Line length = 35 + 19 + 1 = 55 max at this stage,
                     usually 48 */

  /* Create (current) file name.
     To avoid lines longer than 79 characters in the case file, we limit the
     local file name to 79 - l characters. If this limit should be exceeded,
     we try to assign a file name based on the variable index :
     '<case_prefix>.var_***'. If this is not possible due to the length
     of the file name prefix (based on the case name), we choose to
     allow generation of an "incorrect" case file and which can be
     corrected by the user (who will need to rename files) rather than risk
     the possibility of confusing EnSight cases in a same directory due to
     name truncation. */

  remain_len = 79 - strlen(line);
  prefix_len =   strlen(this_case->file_name_prefix)
               - this_case->dir_name_length + 1;
  base_len = strlen(name);
  if (var->time_set > -1)
    postfix_len = 5;
  else
    postfix_len = 0;

  if (   (prefix_len + base_len + postfix_len > remain_len)
      && (prefix_len + 7 + postfix_len <= remain_len)
      && (this_case->n_vars < 1000)) {
    base_len = 7;
    rename_file = true;
  }
  else
    rename_file = false;

  BFT_MALLOC(var->file_name,
             (  this_case->dir_name_length + prefix_len
              + base_len + postfix_len + 1),
             char);
  sprintf(var->file_name, "%s.", this_case->file_name_prefix);

  if (rename_file == true)
    sprintf(var->file_name + strlen(var->file_name),
            "var_%03d", this_case->n_vars);
  else {
    strcat(var->file_name, name);
    for (i = this_case->dir_name_length + prefix_len ;
         i < this_case->dir_name_length + prefix_len + base_len ;
         i++) {
      switch (var->file_name[i]) {
      case '@':
      case ' ':
      case '\t':
      case '!':
      case '#':
      case '*':
      case '^':
      case '$':
      case '/':
        var->file_name[i] = '_';
      default:
        var->file_name[i] = tolower(var->file_name[i]);
      }
    }
  }

  if (var->time_set > -1) {
    sprintf(var->file_name + strlen(var->file_name),
            ".%04d", (this_case->time_set[var->time_set])->n_time_values);
  }

  /* Now we may finish associated case file line */

  BFT_MALLOC(var->case_line,
             (  strlen(line) + strlen(var->file_name)
              - this_case->dir_name_length + 1),
             char);

  strcpy(var->case_line, line);
  strcat(var->case_line, var->file_name + this_case->dir_name_length);

  /* Replace current time index by wildcards */

  if (var->time_set > -1)
    strcpy(var->case_line + strlen(var->case_line) -4, "****");

  /* Finally, associate variable entry in case file */

  BFT_REALLOC(this_case->var, this_case->n_vars, fvm_to_ensight_case_var_t *);

  this_case->var[this_case->n_vars - 1] = var;
  this_case->modified = true;
}

/*----------------------------------------------------------------------------
 * Remove variable entries
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *----------------------------------------------------------------------------*/

static void
_del_vars(fvm_to_ensight_case_t  *const this_case)
{
  int i;

  for (i = 0 ; i < this_case->n_vars ; i++) {

    fvm_to_ensight_case_var_t  *var = this_case->var[i];

    BFT_FREE(var->name);
    BFT_FREE(var->case_line);
    BFT_FREE(var->file_name);

    BFT_FREE(var);

  }

  BFT_FREE(this_case->var);
}

/*----------------------------------------------------------------------------
 * Initialize or update FVM to EnSight Gold geometry file name.
 * The corresponding entry in the case description should be set to NULL
 * prior to first use: in this case, it it created. If a geometry file
 * name already exists, its last characters (representing its time
 * step number) may be replaced, but the string's length should not
 * be modified, so it is modified "in place".
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *----------------------------------------------------------------------------*/

static void
_update_geom_file_name(fvm_to_ensight_case_t  *const this_case)
{
  int geom_index = 0;

  /* Set first time */

  if (this_case->geom_file_name == NULL) {

    char extension[16] = ".geo";

    if (this_case->time_dependency != FVM_WRITER_FIXED_MESH) {
      if (this_case->geom_time_set > -1)
        geom_index =
          (this_case->time_set[this_case->geom_time_set])->n_time_values;
      sprintf(extension, ".geo.%04d", geom_index);
    }

    BFT_MALLOC(this_case->geom_file_name,
               strlen(this_case->file_name_prefix) + strlen(extension) + 1,
               char);
    strcpy(this_case->geom_file_name, this_case->file_name_prefix);
    strcat(this_case->geom_file_name, extension);

  }

  /* Change geometry index */

  else if (this_case->time_dependency != FVM_WRITER_FIXED_MESH) {
    if (this_case->geom_time_set > -1) {
      geom_index =
        (this_case->time_set[this_case->geom_time_set])->n_time_values;
      sprintf(this_case->geom_file_name + strlen(this_case->geom_file_name) - 4,
              "%04d", geom_index);
    }
  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a new case file structure.
 *
 * parameters:
 *   name            <-- case name
 *   dir_prefix      <-- associated local or absolute directory name
 *   time_dependency <-- indicates if and how meshes will change with time
 *
 * returns:
 *   pointer to new case file structure
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_t *
fvm_to_ensight_case_create(const char                   *const name,
                           const char                   *const dir_prefix,
                           const fvm_writer_time_dep_t         time_dependency)
{
  int  i, name_len, prefix_len;

  fvm_to_ensight_case_t   *this_case = NULL;

  /* Create and initialize structure */

  BFT_MALLOC(this_case, 1, fvm_to_ensight_case_t);

  /* Initialize base name and partial file names */

  BFT_MALLOC(this_case->name, strlen(name) + 1, char);
  strcpy(this_case->name, name);
  name_len = strlen(name);

  for (i = 0 ; i < name_len ; i++) {
    if (   (this_case->name[i] == ' ')
        || (this_case->name[i] == '\t'))
      this_case->name[i] = '_';
  }

  if (dir_prefix != NULL)
    prefix_len = strlen(dir_prefix);
  else
    prefix_len = 0;

  this_case->dir_name_length = prefix_len;

  BFT_MALLOC(this_case->case_file_name, prefix_len + name_len + 6, char);
  if (dir_prefix != NULL)
    strcpy(this_case->case_file_name, dir_prefix);
      else
    this_case->case_file_name[0] = '\0';

  for (i = 0 ; i < name_len ; i++)
    this_case->case_file_name[prefix_len + i] = toupper(name[i]);
  this_case->case_file_name[prefix_len + name_len] = '\0';

  BFT_MALLOC(this_case->file_name_prefix,
             strlen(this_case->case_file_name) + 1,
             char);
  strcpy(this_case->file_name_prefix, this_case->case_file_name);
  for (i = 0 ; i < name_len ; i++)
    this_case->file_name_prefix[prefix_len + i]
      = tolower(this_case->case_file_name[prefix_len + i]);

  strcat(this_case->case_file_name, ".case");

  /* Initialize other members */

  this_case->n_parts = 0;
  this_case->part_name = NULL;

  this_case->n_time_sets = 0;
  this_case->time_set = 0;

  this_case->n_vars = 0;
  this_case->var = NULL;

  this_case->geom_time_set = -1; /* no time set yet */

  this_case->time_dependency = time_dependency;

  /* Geometry file name (after time dependency) */

  this_case->geom_file_name = NULL;
  _update_geom_file_name(this_case);

  /* Status information */

  this_case->geom_info_queried = false;

  this_case->modified = true;

  /* Return new case structure */

  return this_case;
}

/*----------------------------------------------------------------------------
 * Destroy a case file structure.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_t *
fvm_to_ensight_case_destroy(fvm_to_ensight_case_t  *this_case)
{
  int  i;

  /* Free names */

  BFT_FREE(this_case->name);
  BFT_FREE(this_case->case_file_name);
  BFT_FREE(this_case->file_name_prefix);

  BFT_FREE(this_case->geom_file_name);

  /* Free other members */

  for (i = 0 ; i < this_case->n_parts ; i++)
    BFT_FREE(this_case->part_name[i]);
  BFT_FREE(this_case->part_name);

  /* Free variable entries */

  _del_vars(this_case);

  /* Free time sets */

  for (i = 0 ; i < this_case->n_time_sets ; i++)
    _time_set_destroy(this_case->time_set[i]);
  BFT_FREE(this_case->time_set);

  /* Free structure and return */

  BFT_FREE(this_case);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

fvm_writer_time_dep_t
fvm_to_ensight_case_get_time_dep(fvm_to_ensight_case_t  *this_case)
{
  return  this_case->time_dependency;
}

/*----------------------------------------------------------------------------
 * Associate new time step with an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *   time_step  <-- time step number
 *   time_value <-- time_value number
 *
 * returns:
 *   0 if no time was added, 1 if a new time was added
 *----------------------------------------------------------------------------*/

int
fvm_to_ensight_case_set_geom_time(fvm_to_ensight_case_t  *const this_case,
                                  const int                     time_step,
                                  const double                  time_value)
{
  int retval = 0;

  if (this_case->geom_time_set == -1) {
    this_case->geom_time_set = this_case->n_time_sets;
    this_case->n_time_sets += 1;
    BFT_REALLOC(this_case->time_set,
                this_case->n_time_sets,
                fvm_to_ensight_case_time_t *);
    this_case->time_set[this_case->geom_time_set] = _time_set_create();
  }

  if (this_case->time_dependency != FVM_WRITER_FIXED_MESH)
    retval = _add_time(this_case->time_set[this_case->geom_time_set],
                       time_step,
                       time_value);

  if (retval > 0) {
    _update_geom_file_name(this_case);
    this_case->geom_info_queried = false;
    this_case->modified = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return current file name and "queried" indicator associated with an
 * EnSight geometry.
 *
 * The "queried" flag in the info structure is set to "false" the first
 * time this function returns a given file name, and to "true" all other
 * times.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   Info structure for geometry file
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_file_info_t
fvm_to_ensight_case_get_geom_file(fvm_to_ensight_case_t  *const this_case)
{
  fvm_to_ensight_case_file_info_t  retval = {NULL, false};

  retval.name    = this_case->geom_file_name;
  retval.queried = this_case->geom_info_queried;

  if (this_case->geom_info_queried == false)
    this_case->geom_info_queried = true;

  return retval;
}

/*----------------------------------------------------------------------------
 * Associate a part name with a case and return its number.
 * If the part was already associated, zero is returned.
 *
 * parameters:
 *   this_case  <-- case structure
 *   part_name  <-- part name
 *
 * returns:
 *   part number in case, or 0 if part already associated
 *----------------------------------------------------------------------------*/

int
fvm_to_ensight_case_add_part(fvm_to_ensight_case_t  *const this_case,
                             const char             *const part_name)
{
  int  i;

  assert(this_case != NULL);

  for (i = 0 ; i < this_case->n_parts ; i++) {
    if (strcmp(part_name, this_case->part_name[i]) == 0)
      break;
  }

  if (i < this_case->n_parts)
    i = 0;

  else {
    this_case->n_parts += 1;
    BFT_REALLOC(this_case->part_name, this_case->n_parts, char *);
    BFT_MALLOC(this_case->part_name[i], strlen(part_name) + 1, char);
    strcpy(this_case->part_name[i], part_name);
    i += 1;
  }

  return i;
}

/*----------------------------------------------------------------------------
 * Return the part number associated with a given part name, or 0.
 *
 * parameters:
 *   this_case  <-- case structure
 *   part_name  <-- part name
 *
 * returns:
 *   part number in case, or 0 if part name is not associated with this case
 *----------------------------------------------------------------------------*/

int
fvm_to_ensight_case_get_part_num(fvm_to_ensight_case_t  *const this_case,
                                 const char             *const part_name)
{
  int  i;

  assert(this_case != NULL);

  for (i = 0 ; i < this_case->n_parts ; i++) {
    if (strcmp(part_name, this_case->part_name[i]) == 0)
      break;
  }

  if (i == this_case->n_parts)
    i = 0;
  else
    i += 1;

  return i;
}

/*----------------------------------------------------------------------------
 * Return current file name and "queried" indicator associated with an
 * EnSight variable.
 *
 * The "queried" flag in the info structure is set to "false" the first
 * time this function returns a given file name, and to "true" all other
 * times.
 *
 * if the corresponding variable or physical time are not present in the
 * structure, the necessary elements are added.
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *   name       <-- variable name
 *   dimension  <-- variable dimension (0: constant, 1: scalar, 3: vector,
 *                  6: symmetrical tensor, 9: asymmetrical tensor)
 *   location   <-- variable definition location (nodes, elements, or particles)
 *   time_step  <-- number of time step to add
 *   time_value <-- associated time value
 *
 * returns:
 *   Info structure for file associated with the variable
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_file_info_t
fvm_to_ensight_case_get_var_file(fvm_to_ensight_case_t       *const this_case,
                                 const char                  *const name,
                                 const int                          dimension,
                                 const fvm_writer_var_loc_t         location,
                                 const int                          time_step,
                                 const double                       time_value)
{
  int i;
  int time_set = -1;
  int var_index = -1;

  fvm_to_ensight_case_var_t  *var = NULL;
  fvm_to_ensight_case_file_info_t  retval = {NULL, false};

  /* First, find if variable is already defined */

  for (i = 0 ; i < this_case->n_vars ; i++) {

    var = this_case->var[i];

    if (strcmp(var->name, name) == 0) {

      /* Variable is already defined, so check consistency */

      if (var->loc != location || var->dim != dimension)
        bft_error(__FILE__, __LINE__, 0,
                  _("A variable with the name \"%s\" has already been\n"
                    "defined, with dimension %d and location type %d,\n"
                    "which conflicts with the current definition with\n"
                    "dimension %d and location type %d.\n"),
                  name, var->dim, (int)(var->loc), dimension, (int)location);

      else if (var->time_set > -1 && time_step < 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("A transient variable with the name \"%s\" has already\n"
                    "been defined, which conflicts with the current constant"
                    " definition.\n"), name);

      else if (var->time_set < 0 && time_step > 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("A constant variable with the name \"%s\" has already\n"
                    "been defined, which conflicts with the current transient"
                    " definition.\n"), name);

      /* At this stage, variable is found, and seems consistent;
         only the time index should need to be updated in the transient case */

      break;

    }

  }

  /* Second, update time step */

  if (time_step > -1) {

    /* Variable based on geometry (and not particle/measured) */
    if (this_case->geom_time_set == -1) {
      this_case->geom_time_set = this_case->n_time_sets;
      this_case->n_time_sets += 1;
      BFT_REALLOC(this_case->time_set,
                  this_case->n_time_sets,
                  fvm_to_ensight_case_time_t *);
      this_case->time_set[this_case->geom_time_set] = _time_set_create();
    }

    time_set = this_case->geom_time_set;
    if (_add_time(this_case->time_set[time_set], time_step, time_value) > 0)
      this_case->modified = true;
  }

  /* If variable found, just update and return file name */

  if (i < this_case->n_vars) {
    retval.name = var->file_name;
    if (var->time_set > -1) {  /* if variable is time-dependant */
      var_index = (this_case->time_set[var->time_set])->n_time_values;
      sprintf(var->file_name + strlen(var->file_name) - 4, "%04d", var_index);
      if (var->last_time_step == time_step)
        retval.queried = true;
      else
        var->last_time_step = time_step;
    }
    else /* if variable is time-independant */
      retval.queried = true;
  }

  /* Otherwise, build variable entry */

  else {
    _add_var(this_case, name, dimension, location, time_set);
    var = this_case->var[this_case->n_vars - 1];
    if (time_step > -1)
      var->last_time_step = time_step;
    retval.name = var->file_name;
  }

 return retval;
}

/*----------------------------------------------------------------------------
 * Write an EnSight Gold case file.
 *
 * This function should only be called by one process in parallel mode.
 *
 * parameters:
 *   this_case  <-- case structure
 *----------------------------------------------------------------------------*/

void
fvm_to_ensight_case_write_case(fvm_to_ensight_case_t  *const this_case)
{
  int          i, j;
  bft_file_t  *f;
  _Bool        write_time_sets = false;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_case->modified == false)
    return;

  /* Open case file (overwrite it if present) */

  f = bft_file_open(this_case->case_file_name,
                    BFT_FILE_MODE_WRITE,
                    BFT_FILE_TYPE_TEXT);

  /* Output FORMAT */

  bft_file_printf(f,
                  "FORMAT\n"
                  "type: ensight gold\n");

  /* Output geometry */

  bft_file_printf(f,
                  "\n"
                  "GEOMETRY\n");

  if (this_case->time_dependency == FVM_WRITER_FIXED_MESH)
    bft_file_printf(f, "model: %s.geo\n",
                    this_case->file_name_prefix + this_case->dir_name_length);

  else if (this_case->time_dependency == FVM_WRITER_TRANSIENT_COORDS)
    bft_file_printf(f, "model: %d %s.geo.****  change_coords_only\n",
                    this_case->geom_time_set + 1,
                    this_case->file_name_prefix + this_case->dir_name_length);

  else
    bft_file_printf(f, "model: %d %s.geo.****\n",
                    this_case->geom_time_set + 1,
                    this_case->file_name_prefix + this_case->dir_name_length);

  /* Output variables */

  if (this_case->n_vars > 0) {

    bft_file_printf(f,
                    "\n"
                    "VARIABLE\n");

    for (i = 0 ; i < this_case->n_vars ; i++) {
      const fvm_to_ensight_case_var_t  *var = this_case->var[i];
      bft_file_printf(f, "%s\n", var->case_line);
    }

  }

  /* Output time section */

  if (this_case->n_time_sets > 0) {
    for (i = 0 ; i < this_case->n_time_sets ; i++) {
      if ((this_case->time_set[i])->n_time_values > 0) {
        write_time_sets = true;
        break;
      }
    }
  }

  if (write_time_sets == true) {

    bft_file_printf(f,
                    "\n"
                    "TIME\n");

    for (i = 0 ; i < this_case->n_time_sets ; i++) {

      const fvm_to_ensight_case_time_t  *ts = this_case->time_set[i];

      bft_file_printf(f, "time set:              %d\n", i+1);
      bft_file_printf(f, "number of steps:       %d\n", ts->n_time_values);
      bft_file_printf(f, "filename start number: 1\n");
      bft_file_printf(f, "filename increment:    1\n");
      bft_file_printf(f, "time values:\n");

      for (j = 0 ; j < ts->n_time_values ; j++)
        bft_file_printf(f, "            %g\n", ts->time_value[j]);

    }

  }

  /* Close case file */

  bft_file_free(f);

  this_case->modified = false;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
