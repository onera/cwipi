.. _concepts:

This section underlines the structure of a coupling setup in this new version of CWIPI.
The concepts will be detailed working with the following coupling scheme:

.. image:: ./images/schema_basic_coupling.svg

Coupling
--------

To set up the coupling between `Solver1` and `Solver2`, create a Coupling instance and associate it with general information such as the geometric dimension of the coupling interface.

Mesh
----

Then specify the coupling interface geometry (see :ref:`Define mesh`).
In this case, we set a 2D triangle and quadrangle mesh for `Solver1` and a polygon mesh for `Solver2`.
After setting the mesh coordinates, a so called block of the mesh elements should be added.
This means that in the mesh instance a block for the given type of elements will be added.
After creating this block, the mesh element data can be given using the add function for standard elements (``CWP_Mesh_interf_block_std_set`` in C) for `Solver1` and for polygons (``CWP_Mesh_interf_f_poly_block_set`` in C) for `Solver2`.

For CWIPI to be able to do the internal geometric computations on the mesh, it must be "finalized" (``CWP_Mesh_interf_finalize`` in C).

Convention for standard elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. list-table:: Numbering convention used in CWIPI for standard elements
  :widths: 50 50

  * - .. figure:: ../../../images/pdm_mesh_nodal_point.svg
        :alt: PDM_MESH_NODAL_POINT

        ``PDM_MESH_NODAL_POINT``

    - .. figure:: ../../../images/pdm_mesh_nodal_bar2.svg
        :alt: PDM_MESH_NODAL_BAR2

        ``PDM_MESH_NODAL_BAR2``


  * - .. figure:: ../../../images/pdm_mesh_nodal_tria3.svg
        :alt: PDM_MESH_NODAL_TRIA3

        ``PDM_MESH_NODAL_TRIA3``


    - .. figure:: ../../../images/pdm_mesh_nodal_quad4.svg
        :alt: PDM_MESH_NODAL_QUAD4

        ``PDM_MESH_NODAL_QUAD4``


  * - .. figure:: ../../../images/pdm_mesh_nodal_tetra4.svg
        :alt: PDM_MESH_NODAL_TETRA4

        ``PDM_MESH_NODAL_TETRA4``


    - .. figure:: ../../../images/pdm_mesh_nodal_pyram5.svg
        :alt: PDM_MESH_NODAL_PYRAM5

        ``PDM_MESH_NODAL_PYRAM5``


  * - .. figure:: ../../../images/pdm_mesh_nodal_prism6.svg
        :alt: PDM_MESH_NODAL_PRISM6

        ``PDM_MESH_NODAL_PRISM6``


    - .. figure:: ../../../images/pdm_mesh_nodal_hexa8.svg
        :alt: PDM_MESH_NODAL_HEXA8

        ``PDM_MESH_NODAL_HEXA8``



Fields
------

A Field models a physical quantity with the geometrical support of a previously defined mesh.

The class provides the following main services:

- data storage description (type, interlacing, number of components)
- send, recv.

It is mandatory to define the interface mesh *before* creating field instances.
The degrees-of-freedom (dof) of a Field can either be located at mesh nodes, cell centers or user-defined points.
There can be no more than one user-defined point cloud per Coupling object.

For `Solver1` a field instance for sending the temperature will be created and another instance for receiving the pressure.
For `Solver2` the opposite will be done.

There are two ways to store field components:

* CWP_FIELD_STORAGE_INTERLACED   : The number of components is constant for each element. The field is stored according to this pattern :math:`(c_{1,1} ... c_{s,1} ... c_{1,n} ... c_{s,n})` , where :math:`s` is the number of components and :math:`n` the number of field elements;
* CWP_FIELD_STORAGE_INTERLEAVED  : The number of components is constant for each element. The field is stored according to this pattern :math:`(c_{1,1} ... c_{1,n} ... c_{s,1} ... c_{s,n})` , where :math:`s` is the number of components and :math:`n` the number of field elements. In this mode,


Control Parameters
------------------

*Not documented yet*
