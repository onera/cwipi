.. _old_to_new:

Old to New
##########

Version 1.* brings major changes in the CWIPI library API compared to version 0.*. This new version is not only a API
change, it mainly is a feature enrichment. Thus user are still able to do what they did with the old version and even more.
Here we focus on the equivalent features between the two versions.

Initialize and Finalize
=======================

All function call are still framed by a call to the initialization function and the finalize function (the last didn't change).
For initialization the are some new arguments since the library can do more then before. Indeed, it handles the case that several codes
are running on the same MPI rank. That is the reason why one need to give the number of codes and a list of arguments at initialization.
Moreover, one can choose which MPI ranks they want to participate in the coupling. That is done by giving an array of
``CWP_Status_t`` telling if the MPI rank is used or not for each code on the MPI rank. The ``time_init`` array is an argument in anticipation of
future developments, thus unused for now. The `intra_comms` array output argument gives for each code on the MPI rank
the MPI communicator to communicate between the MPI ranks of that code.

Create a coupling
=================

Let us focus on the arguments in the old version of CWIPI first. They are now called in seperate functions. In the ``CWP_Cpl_create`` one
still sets the lcoal code name, the coupled code name, the coupling type (``comm_type``) and the entities dimension. Since more
interpolation algorithms are available, their related property is not necessary the geometric tolerence. Thus we set this information
in the ``CWP_Spatial_interp_property_set`` function. Almost all information about the mesh is now done using the ``CWP_Mesh_interf_*`` functions.
While creating the coupling, one only tells if the mesh is gonna be deformed/changed or not. The ouput for visualisation is
instrumented by ``CWP_Visu_set``. So what's new? Biggest change is the large amount of available spatial interpolation algorithms.
``recv_freq_type`` is also an anticipation of future developments.

Exchange a field
================

The exchange of fields is no more done in one function since this new API gives more flexibility. The object oriented aspect
comes out more in this version. The user creates a coupling between two codes each having a mesh on which several fields can be defined.
The non blocking exchange functions were kept and are the way to go.

Summary
=======

.. image:: ./images/coupling.png
   :scale: 75%

.. image:: ./images/mesh.png
   :scale: 75%

.. image:: ./images/field.png
   :scale: 75%
