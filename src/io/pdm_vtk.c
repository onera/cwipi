
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
#include "pdm_vtk.h"

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
static int _vtk_elt_type
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order
 )
{
  //assert (order <= 2);

  int vtk_elt_type;
  switch (elt_type) {
  case PDM_MESH_NODAL_POINT:
    vtk_elt_type = 1;
    break;
  case PDM_MESH_NODAL_BAR2:
    if (order == 1) {
      vtk_elt_type = 3;
      //} else if (order == 2) {
      // vtk_elt_type = 21;
    } else {
      vtk_elt_type = 68;
    }
    break;
  case PDM_MESH_NODAL_TRIA3:
    if (order == 1) {
      vtk_elt_type = 5;
      //} else if (order == 2) {
      //vtk_elt_type = 22;
    } else {
      vtk_elt_type = 69;
    }
    break;
  case PDM_MESH_NODAL_QUAD4:
    if (order == 1) {
      vtk_elt_type = 9;
      //} else if (order == 2) {
      //  vtk_elt_type = 23;
    } else {
      vtk_elt_type = 70;
    }
    break;
  case PDM_MESH_NODAL_TETRA4:
    if (order == 1) {
      vtk_elt_type = 10;
      //} else if (order == 2) {
      //vtk_elt_type = 24;
    } else {
      vtk_elt_type = 71;
    }
    break;
  case PDM_MESH_NODAL_PYRAMID5:
    if (order == 1) {
      vtk_elt_type = 14;
    } else if (order == 2) {
      vtk_elt_type = 27;
    } else {
      vtk_elt_type = 66;
    }
    break;
  case PDM_MESH_NODAL_PRISM6:
    if (order == 1) {
      vtk_elt_type = 13;
      //} else if (order == 2) {
      //vtk_elt_type = 26;
    } else {
      vtk_elt_type = 73;
    }
    break;
  case PDM_MESH_NODAL_HEXA8:
    if (order == 1) {
      vtk_elt_type = 12;
      //} else if (order == 2) {
      //vtk_elt_type = 25;
    } else {
      vtk_elt_type = 72;
    }
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "type %d is not a valid std elt type\n", elt_type);
  }

  return vtk_elt_type;
}



static void
_ijk_to_vtk
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order,
 int                        idx[]
 )
{
  switch (elt_type) {

  case PDM_MESH_NODAL_POINT:
    idx[0] = 0;
    break;

  case PDM_MESH_NODAL_BAR2:
    idx[0] = 0;
    idx[1] = order;
    for (int i = 1; i < order; i++) {
      idx[i+1] = i;
    }
    break;

  case PDM_MESH_NODAL_TRIA3:
    if (order == 1) {
      idx[2] = 2;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[2] = 5;
      idx[5] = 3; idx[4] = 4;
      idx[0] = 0; idx[3] = 1; idx[1] = 2;
    } else if (order == 3) {
      idx[2] = 9;
      idx[7] = 7; idx[6] = 8;
      idx[8] = 4; idx[9] = 5; idx[5] = 6;
      idx[0] = 0; idx[3] = 1; idx[4] = 2; idx[1] = 3;
    } else if (order == 4) {
      idx[ 2] = 14;
      idx[ 9] = 12; idx[ 8] = 13;
      idx[10] =  9; idx[14] = 10; idx[ 7] = 11;
      idx[11] =  5; idx[12] =  6; idx[13] =  7; idx[ 6] =  8;
      idx[ 0] =  0; idx[ 3] =  1; idx[ 4] =  2; idx[ 5] =  3; idx[ 1] =  4;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "TRIA VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_QUAD4:
    if (order == 1) {
      idx[3] = 2; idx[2] = 3;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[3] = 6; idx[6] = 7; idx[2] = 8;
      idx[7] = 3; idx[8] = 4; idx[5] = 5;
      idx[0] = 0; idx[4] = 1; idx[1] = 2;
    } else if (order == 3) {
      idx[ 3] = 12; idx[ 8] = 13; idx[ 9] = 14; idx[ 2] = 15;
      idx[11] =  8; idx[14] =  9; idx[15] = 10; idx[ 7] = 11;
      idx[10] =  4; idx[12] =  5; idx[13] =  6; idx[ 6] =  7;
      idx[ 0] =  0; idx[ 4] =  1; idx[ 5] =  2; idx[ 1] =  3;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "QUAD VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_TETRA4:
    if (order == 1) {
      idx[3] = 3;

      idx[2] = 2;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[3] = 9;

      idx[9] = 8;
      idx[7] = 6; idx[8] = 7;

      idx[2] = 5;
      idx[6] = 3; idx[5] = 4;
      idx[0] = 0; idx[4] = 1; idx[1] = 2;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "TETRA VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_PRISM6:
    if (order == 1) {
      idx[5] = 4;
      idx[3] = 3; idx[4] = 5;

      idx[2] = 1;
      idx[0] = 0; idx[1] = 2;
    } else if (order == 2) {
      idx[ 5] = 14;
      idx[11] = 13; idx[10] = 16;
      idx[ 3] = 12; idx[ 9] = 15; idx[ 4] = 17;

      idx[14] =  8;
      idx[17] =  7; idx[16] = 10;
      idx[12] =  6; idx[15] =  9; idx[13] = 11;

      idx[ 2] =  2;
      idx[ 8] =  1; idx[ 7] =  4;
      idx[ 0] =  0; idx[ 6] =  3; idx[ 1] =  5;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "PRISM VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_HEXA8:
    if (order == 1) {
      idx[7] = 6; idx[6] = 7;
      idx[4] = 4; idx[5] = 5;

      idx[3] = 2; idx[2] = 3;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[ 7] = 24; idx[14] = 25; idx[ 6] = 26;
      idx[15] = 21; idx[21] = 22; idx[13] = 23;
      idx[ 4] = 18; idx[12] = 19; idx[ 5] = 20;

      idx[19] = 15; idx[24] = 16; idx[18] = 17;
      idx[25] = 12; idx[26] = 13; idx[23] = 14;
      idx[16] =  9; idx[22] = 10; idx[17] = 11;

      idx[ 3] =  6; idx[10] =  7; idx[ 2] =  8;
      idx[11] =  3; idx[20] =  4; idx[ 9] =  5;
      idx[ 0] =  0; idx[ 8] =  1; idx[ 1] =  2;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "HEXA VTK ordering not implemented for order %d\n", order);
    }
    break;
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_vtk_write_boxes
(
 const char        *filename,
 const int          n_box,
 const double      *box_extents,
 const PDM_g_num_t *box_g_num
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nboxes\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_box);
  for (int i = 0; i < n_box; i++) {
    const double *e = box_extents + 6*i;
    fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_box, 9*n_box);
  for (int i = 0; i < n_box; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_box);
  for (int i = 0; i < n_box; i++) {
    fprintf(f, "12\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_box);
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int i = 0; i < n_box; i++) {
    fprintf(f, ""PDM_FMT_G_NUM"\n", box_g_num[i]);
  }

  fclose(f);
}

void
PDM_vtk_write_circles
(
 const char        *filename,
 const int          n_circles,
 const double      *center,
 const double      *radius,
 const PDM_g_num_t *g_num,
 const int         *color,
 const int          resolution
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "circles\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  double step = 2. * PDM_PI / (double) (resolution - 1);

  fprintf(f, "POINTS %d double\n", n_circles * resolution);
  for (int i = 0; i < n_circles; i++) {
    for (int j = 0; j < resolution; j++) {
      fprintf(f, "%f %f %f\n",
              center[3*i]     + radius[i]*cos(j*step),
              center[3*i + 1] + radius[i]*sin(j*step),
              center[3*i + 2]);
    }
  }

  fprintf(f, "CELLS %d %d\n", n_circles, n_circles * (resolution + 1));
  for (int i = 0; i < n_circles; i++) {
    fprintf(f, "%d \n", resolution);
    for (int j = 0; j < resolution-1; j++) {
      fprintf(f, "%d ", resolution*i + j);
    }
    fprintf(f, "%d\n", resolution*i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_circles);
  for (int i = 0; i < n_circles; i++) {
    fprintf(f, "4\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_circles);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  else if (color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_circles);
    fprintf(f, "SCALARS color long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  fclose(f);
}




void
PDM_vtk_write_polydata
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[],
 const int          face_color[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + face_vtx_idx[n_face]);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (face_g_num != NULL && face_color == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
     }
  } else if (face_color != NULL && face_g_num == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d\n", face_color[i]);
    }
  } else if (face_g_num != NULL && face_color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "FIELD field 2\n");
    fprintf(f, "face_gnum 1 %d long\n", n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM" ", face_g_num[i]);
    }
    fprintf(f, "\nface_color 1 %d int\n", n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d ", face_color[i]);
    }
  }




  fclose(f);
}




void
PDM_vtk_write_point_cloud
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          color[]
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "point cloud\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_vtx, 2*n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "1\n");
  }

  if (vtx_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  } else if (color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  fclose(f);
}



void
PDM_vtk_write_lines
(
 const char        *filename,
 const int          n_line,
 const double      *coord,
 const PDM_g_num_t *g_num,
 const int         *color
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "lines\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_line * 2);
  for (int i = 0; i < 2*n_line; i++) {
    fprintf(f, "%f %f %f\n",
            coord[3*i    ],
            coord[3*i + 1],
            coord[3*i + 2]);
  }

  fprintf(f, "CELLS %d %d\n", n_line, n_line * 3);
  for (int i = 0; i < n_line; i++) {
    fprintf(f, "%d %d %d\n", 2, 2*i, 2*i+1);
  }

  fprintf(f, "CELL_TYPES %d\n", n_line);
  for (int i = 0; i < n_line; i++) {
    fprintf(f, "3\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_line);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_line; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  else if (color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_line);
    fprintf(f, "SCALARS color long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_line; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  fclose(f);
}




void
PDM_vtk_write_std_elements
(
 const char                 *filename,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_ifield,
 const char                 *elt_ifield_name[],
 const int                  *elt_ifield[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (elt_type, 1);

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt * (1 + n_vtx_elt));
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d", n_vtx_elt);
    for (int j = 0; j < n_vtx_elt; j++) {
      fprintf(f, " %d", elt_vtx[n_vtx_elt*i + j] - 1);
    }
    fprintf(f, "\n");
  }


  int vtk_elt_type = _vtk_elt_type (elt_type, 1);
  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d\n", vtk_elt_type);
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (elt_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_elt);
    fprintf(f, "SCALARS elt_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", elt_g_num[i]);
     }
  }

  if (n_elt_ifield > 0) {
    assert (elt_ifield != NULL);

    if (elt_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_elt);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_ifield);
    for (int i = 0; i < n_elt_ifield; i++) {
      assert (elt_ifield[i] != NULL);
      assert (elt_ifield_name[i] != NULL);

      fprintf(f, "%s 1 %d int\n", elt_ifield_name[i], n_elt);
      for (int j = 0; j < n_elt; j++) {
        fprintf(f, "%d ", elt_ifield[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}


void
PDM_vtk_write_std_elements_double
(
 const char                 *filename,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_field,
 const char                 *elt_field_name[],
 const double               *elt_field[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (elt_type, 1);

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt * (1 + n_vtx_elt));
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d", n_vtx_elt);
    for (int j = 0; j < n_vtx_elt; j++) {
      fprintf(f, " %d", elt_vtx[n_vtx_elt*i + j] - 1);
    }
    fprintf(f, "\n");
  }


  int vtk_elt_type = _vtk_elt_type (elt_type, 1);
  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d\n", vtk_elt_type);
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (elt_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_elt);
    fprintf(f, "SCALARS elt_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", elt_g_num[i]);
     }
  }

  if (n_elt_field > 0) {
    assert (elt_field != NULL);

    if (elt_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_elt);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_field);
    for (int i = 0; i < n_elt_field; i++) {
      assert (elt_field[i] != NULL);
      assert (elt_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", elt_field_name[i], n_elt);
      for (int j = 0; j < n_elt; j++) {
        fprintf(f, "%lf ", elt_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}















void
PDM_vtk_write_std_elements_ho
(
 const char                 *filename,
 const int                   order,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_field,
 const char                 *elt_field_name[],
 const double               *elt_field[]
 )
{
  //assert (order < 3);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (elt_type, order);
  /*int *vtk_idx = malloc (sizeof(int) * n_vtx_elt);
  _ijk_to_vtk (elt_type,
                order,
                vtk_idx);*/

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt * (1 + n_vtx_elt));
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d", n_vtx_elt);
    for (int j = 0; j < n_vtx_elt; j++) {
      //fprintf(f, " %d", elt_vtx[n_vtx_elt*i + vtk_idx[j]] - 1);
      fprintf(f, " %d", elt_vtx[n_vtx_elt*i + j] - 1);
    }
    fprintf(f, "\n");
  }


  int vtk_elt_type = _vtk_elt_type (elt_type, order);
  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d\n", vtk_elt_type);
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (elt_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_elt);
    fprintf(f, "SCALARS elt_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", elt_g_num[i]);
     }
  }

  if (n_elt_field > 0) {
    assert (elt_field != NULL);

    if (elt_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_elt);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_field);
    for (int i = 0; i < n_elt_field; i++) {
      assert (elt_field[i] != NULL);
      assert (elt_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", elt_field_name[i], n_elt);
      for (int j = 0; j < n_elt; j++) {
        fprintf(f, "%lf ", elt_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}




void
PDM_vtk_write_ellipses
(
 const char        *filename,
 const int          n_ellipse,
 const double      *center,
 const double      *axes,
 const double      *radii,
 const PDM_g_num_t *g_num,
 const int         *color,
 const int          resolution
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "ellipses\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  double step = 2. * PDM_PI / (double) (resolution - 1);

  fprintf(f, "POINTS %d double\n", n_ellipse * resolution);
  for (int i = 0; i < n_ellipse; i++) {
    for (int j = 0; j < resolution; j++) {
      double x = radii[2*i  ] * cos(j*step);
      double y = radii[2*i+1] * sin(j*step);
      fprintf(f, "%f %f %f\n",
              center[3*i    ] + x*axes[6*i  ] + y*axes[6*i+3],
              center[3*i + 1] + x*axes[6*i+1] + y*axes[6*i+4],
              center[3*i + 2] + x*axes[6*i+2] + y*axes[6*i+5]);
    }
  }

  fprintf(f, "CELLS %d %d\n", n_ellipse, n_ellipse * (resolution + 1));
  for (int i = 0; i < n_ellipse; i++) {
    fprintf(f, "%d \n", resolution);
    for (int j = 0; j < resolution-1; j++) {
      fprintf(f, "%d ", resolution*i + j);
    }
    fprintf(f, "%d\n", resolution*i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_ellipse);
  for (int i = 0; i < n_ellipse; i++) {
    fprintf(f, "4\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_ellipse);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_ellipse; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  else if (color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_ellipse);
    fprintf(f, "SCALARS color long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_ellipse; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  fclose(f);
}
