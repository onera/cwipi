.. _old_to_new:

Old to New
##########

Version 1.x introduces major changes in the CWIPI library API compared to version 0.x.
This API change enables much more advanced features.
However, users should be able to easily reproduce their current applications with this new version.

In order to ease this transition, we focus here on the equivalence between the features of versions 0.x and 1.x.

Initialize and Finalize
=======================

All function calls are still framed by a call to the initialization and finalization functions.
Version 1.x supports multiple codes running on a single MPI rank.
Besides, one can select for each code which MPI ranks will be available for CWIPI.
Future versions will support time interpolation.
In anticipation for this feature, the initial time of each local code must be provided. (?????)

In summary, from version 1.x onwards the following additional arguments are required at CWIPI initialization:
   - ``n_code``: the number of codes executed on current MPI rank ;
   - ``code_names``: the list of local code names ;
   - ``is_active_rank``: this array indicates whether current MPI rank will participate in the coupling for each local code ;
   - ``time_init``: the array of initial times of local codes.

.. For initialization the are some new arguments since the library can do more then before. Indeed, it handles the case that several codes
.. are running on the same MPI rank. That is the reason why one need to give the number of codes and a list of arguments at initialization.
.. Moreover, one can choose which MPI ranks they want to participate in the coupling. That is done by giving an array of
.. ``CWP_Status_t`` telling if the MPI rank is used or not for each code on the MPI rank. The ``time_init`` array is an argument in anticipation of
.. future developments, thus unused for now. The `intra_comms` array output argument gives for each code on the MPI rank
.. the MPI communicator to communicate between the MPI ranks of that code.

Create a coupling
=================

The key concepts used in CWIPI have been revisited in version 1.x.
Most importantly, the notion of *Field* as an object has been introduced (see section `Exchange fields`_).
The ``solver_type`` argument is now an attribute of the Field object.
Thus, a single Coupling object can now be used to exchange multiple fields with different degrees-of-freedom (mesh nodes, cell centers or user-defined target points).
All output-specfic arguments are now passed to ``CWP_Visu_set``.

.. Some of the coupling parameters used in version 0.x are now attributes of the Fields and are thus no longer specified at the Coupling creation.
.. Typically, the ``solver_type`` argument as well as all output-specfic arguments are now passed to separate functions.

On the other hand, some new arguments are now required.
First, as a coupling involves exactly two codes, the identifier of the local code involved in the coupling must be specified.
Next, since multiple spatial interpolation algorithms are now available, ``spatial_interp`` determines which one will be used for the coupling.
The specific properties for this algorithm (such as the geometric tolerance) are now defined via the function ``CWP_Spatial_interp_property_set``.
Since version 1.x supports multiple partitions (or subdomains) per MPI rank, ``n_part`` indicates the number of partitions of the coupling interface for the local code, on current MPI rank.
Finally, ``recv_freq_type`` will enable the user to choose a specific time interpolation scheme (not implemented yet).

The following table establishes the equivalence between the arguments that are essentially unchanged:

========================= =========================
**Version 0.x**           **Version 1.x**
========================= =========================
``coupling_name``         ``cpl_id``
``entitiesDim``           ``entities_dim``
``coupling_type``         ``comm_type``
``coupled_application``   ``coupled_code_name``
``mesh_type``             ``displacement``
========================= =========================


.. Let us focus on the arguments in the old version of CWIPI first.
.. They are now called in separate functions.
.. In the ``CWP_Cpl_create`` one still sets the lcoal code name, the coupled code name, the coupling type (``comm_type``) and the entities dimension.
.. Since more interpolation algorithms are available, their related property is not necessary the geometric tolerence.
.. Thus we set this information in the ``CWP_Spatial_interp_property_set`` function.
.. Almost all information about the mesh is now done using the ``CWP_Mesh_interf_*`` functions.
.. While creating the coupling, one only tells if the mesh is gonna be deformed/changed or not.
.. The ouput for visualisation is instrumented by ``CWP_Visu_set``.
.. So what's new?
.. Biggest change is the large amount of available spatial interpolation algorithms.
.. ``recv_freq_type`` is also an anticipation of future developments.


Exchange fields
================

The exchange of fields is no longer performed in one function call since this new API gives more flexibility.
The object-oriented aspect comes out more in this version.
The user creates a coupling between two codes each having a mesh on which several fields can be defined.
The non-blocking exchange functions were kept and are the way to go.

Summary
=======

.. image:: ./images/coupling.png
   :scale: 75%

.. image:: ./images/mesh.png
   :scale: 75%

.. image:: ./images/field.png
   :scale: 75%
