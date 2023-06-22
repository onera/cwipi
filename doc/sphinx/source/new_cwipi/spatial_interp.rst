.. _spatial_interp:


From version 1.0 onwards, CWIPI offers multiple spatial interpolation methods.
A specific spatial interpolation method can be affected to each coupling object.


 (à détailler...)

Interpolation from nearest neighbors
------------------------------------

Weighted Least Square interpolation from nearest neighbors (à détailler...)

parameters:
  - ``n_neighbors``: number of neighbors to interpolate from
  - ``polyfit_degree``: degree of polynomial fit

``CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES``: each target dof is mapped to its ``n_neighbors`` nearest source dofs

``CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES``: each source dof is mapped to its ``n_neighbors`` nearest target dofs (note that some targets may be mapped to zero source)



Interpolation from mesh intersection
------------------------------------

.. .. math::
..    (a + b)^2  &=  (a + b)(a + b) \\
..               &=  a^2 + 2ab + b^2

``CWP_SPATIAL_INTERP_FROM_INTERSECTION``

parameters:
  - ``tolerance``: relative tolerance for bounding box inflation (especially useful for surface interface meshes)

Interpolation from localization
-------------------------------

The target degrees-of-freedom are localized in the source mesh elements.

``CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT``
``CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE``
``CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE``

parameters:
  - ``tolerance``: relative tolerance for bounding box inflation
