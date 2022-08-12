#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2012-2017  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/


#include "cwp.h"
#include "pdm_mpi.h"
#include "pdm_writer.h"
#include "field.hxx"
#include <string>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class Visu visualization.hxx "visualization.hxx"
   * \brief Interface mesh
   *
   *  This class defines the interface mesh objects.
   *
   */
  class Mesh;
  class Visu {

  public:

   /**
    * \brief Visu constructor
    *
    * Construct the CWIPI visualization by using paradigm pario methods.
    *
    * \param [in] MPIComm Coupling Communicator.
    */

    Visu(const MPI_Comm &MPIComm, const CWP_Dynamic_mesh_t topology);

    /**
     * \brief Visu destructor
     *
     * Destroy the Visu object.
     *
     */

    ~Visu();


     /**
     * \brief Create the Visu
     *
     * \param [in]  freq             Output frequency
     * \param [in]  format           Output format to visualize exchanged fieldsDouble
     *                               on the coupled mesh. Choice between :
     *                               - "EnSight Gold"
     *                               - "MED_ficher"
     *                               - "CGNS"
     *                               .
     * \param [in]  format_option   Output options "opt1, opt2, ..." :
     *                         - text               output text files
     *                         - binary             output binary files (default)
     *                         - big_endian         force binary files
     *                                              to big-endian
     *                         - discard_polygons   do not output polygons
     *                                              or related values
     *                         - discard_polyhedra  do not output polyhedra
     *                                              or related values
     *                         - divide_polygons    tesselate polygons
     *                                              with triangles
     *                         - divide_polyhedra   tesselate polyhedra
     *                                              with tetrahedra and pyramids
     *                                              (adding a vertex near
     *                                               each polyhedron's center)
     *
     * \param [in] output_dir  Output directory
     * \param [in] output_name Output name
     *
     */

    void VisuCreate(const int          freq,
                    CWP_Visu_format_t  format,
                    const char        *format_option,
                    char              *output_dir,
                    char              *output_name);

   /**
    * \brief Create the visu geom
    *
    * \param [in] n_part       Number of mesh partition
    *
    */

    void GeomCreate(int n_part);

     /**
     * \brief Set the Visu geom coordinates
     *
     * \param [in] id_part       Index of the mesh partition
     * \param [in] n_pts        Number of partition vertices
     * \param [in] coords       Array containing the verticies coordinates
     * \param [in] global_num   Global numbering of the vertices
     *
     */

    void GeomCoordSet(int id_part,
                      int n_pts,
                      double *coords,
                      CWP_g_num_t *global_num);



     /**
     * \brief Addition of a block (set of cells) to a Visu geom partition
     *
     * This function add a block to a Visu geom partition
     *
     * \param [in] blockType  Type of the block addition
     *
     * \return block_id  Block Identifier
     */

    int  GeomBlockAdd(CWP_Block_t blockType);

   /**
    * \brief Write the Visu geom
    *
    */

    void GeomWrite(Mesh* mesh);

   /**
    * \brief Set a standard block to a Visu geom partition
    *
    * \param [in] id_part     Partition identifier
    * \param [in] id_block    Block identifier
    * \param [in] n_elt       Number of block elements
    * \param [in] connec      Vertices to elements connectivity
    * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
    *
    */

    void GeomBlockStdSet (int id_block,
                          int id_part,
                          int n_elt,
                          int *connec,
                          CWP_g_num_t *global_num);


   /**
    * \brief Set a standard block to a Visu geom partition
    *
    * \param [in] id_part     Partition identifier
    * \param [in] id_block    Block identifier
    * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
    *
    */

    void GeomBlockGNumMeshSet (int id_block,
                               int id_part,
                               CWP_g_num_t *global_num);

     /**
     * \brief Set a face polygon block to a Visu geom partition
     *
     * \param [in] id_part     Partition identifier
     * \param [in] id_block    Block identifier
     * \param [in] n_elt       Number of block elements
     * \param [in] connec_idx  Vertices to elements connectivity index
     * \param [in] connec      Vertices to elements connectivity
     * \param [in] global_num  Global numbering of the block vertices (or NULL)
     *
     */

     void GeomBlockPoly2D(int id_block,
                          int id_part,
                          int n_elt,
                          int *connec_idx,
                          int *connec,
                          CWP_g_num_t *global_num);

     /**
     * \brief Set a face polhedron block to a Visu geom partition
     *
     * \param [in] id_part           Partition identifier
     * \param [in] id_block          Block identifier
     * \param [in] n_elts            Number of block elements
     * \param [in] n_faces           Number of faces elements (or NULL)
     * \param [in] connec_faces_idx  Vertices to faces connectivity index
     * \param [in] connec_faces      Vertices to faces connectivity
     * \param [in] connec_cells_idx  Faces to cells connectivity index
     * \param [in] connec_cells      Faces to cells connectivity
     * \param [in] global_num        Global numbering of the block vertices (or NULL)
     *
     */

     void GeomBlockPoly3D(int              id_block,
                          int              id_part,
                          int              n_elts,
                          int              n_faces,
                          int              connec_faces_idx[],
                          int              connec_faces[],
                          int              connec_cells_idx[],
                          int              connec_cells[],
                          CWP_g_num_t      global_num[]
                        );

   /**
    * \brief Indicate the beginning of the writing step
    *
    * \param [in] physical_time Physical time corresponding to the Visu entity
    *
    */

    void WriterStepBegin(double physical_time, Mesh* mesh);

   /**
    * \brief Indicate the end of the writing step
    *
    */

    void WriterStepEnd();

    void fieldDataSet(Field* field, CWP_Field_map_t map_type, int i_part);

    void fieldDataFree(Field* field);

    void GeomFree();

    void WriterFieldCreate(Field* field);


    inline double physicalTimeGet();

    void WriterField(Field* field, int* n_ref_values, int **ref_values,  const CWP_Field_map_t  map_type);

   /**
    * \brief Convert a CWIPI block type CWP_Block_t to a PDM_writer block type
    *
    * \param [in] CWP_block_type  A CWIPI Block type
    *
    */

    PDM_writer_elt_geom_t PdmWriterBlockTypeFromCwpBlockType(CWP_Block_t CWP_block_type);

   /**
    * \brief Return True if the Visu has been created.
    *
    *
    */

    inline bool isCreated();

  private:
    PDM_writer_t*        _visu_id;           /*!< Visualization identifier */
    int                  _visu_mesh_id;      /*!< Visualization identifier */
    //int                  _freq;            /*!< Visualization frequency */
    char                *_output_dir;        /*!< Output directory */
    char                *_output_name;       /*!< Output Name */
    PDM_MPI_Comm         _pdmComm;           /*!< Paradigm communicator */
    PDM_writer_status_t  _divide_polygons;   /*!< Option to divide polygons */
    PDM_writer_status_t  _divide_polyhedra;  /*!< Option to divide polyhedra */
    bool                 _visuCreated;       /*!< True if the creation has be done */
    double               _physical_time;     /*!< Physical time for visualization */
    int                  _n_part;            /*!< Number of mesh partition */
    CWP_Dynamic_mesh_t   _topology;          /*!< Mesh topology */

    int _id_partitioning_field  ;
    std::vector<double*> _partitioning_field_data;

    int _id_ranking_field  ;
    std::vector<double*> _ranking_field_data;

    int _id_blocking_field  ;
    std::vector<double*> _blocking_field_data;

  };

  bool Visu::isCreated() {
    return _visuCreated;
  }


  double Visu::physicalTimeGet() {
    return _physical_time;
  }


}

/**
 * \endcond
 */

#endif //__VISUALIZATION_H__
