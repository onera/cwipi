*********************
CWIPI's documentation
*********************

**CWIPI** (Coupling With Interpolation Parallel Interface) is a C++/C parallel coupling library under LGPL.
It allows fully parallel data exchanges based on distributed mesh definition. Those meshes are differently partitionned
on several processes. The interpolation is done between non coincident meshes on-the-fly. Arbitrary many codes can be coupled
using this library. The philosophy of CWIPI is to let the parallelism be fully transparent for the user.

The library does not rely on a central coupling instance. The coupling is minimally invasive because if suffices
to call primitives of the library in the codes. Still, such a central supervisor can be set up either by using
OpenPALM software (CERFACS) or by writing a Python supervisor script calling the Python interface of CWIPI.

Summary
-------


:ref:`Old CWIPI <old_cwipi>` is the documentation of the 0.* version of CWIPI. It describes most of the high level APIs provided by CWIPI.

:ref:`New CWIPI <new_cwipi>` is the documentation of the 1.* version of CWIPI. It describes most of the high level APIs provided by CWIPI.

:ref:`Client-Server <client_server_cwipi>` is the documentation of a client-server mode of CWIPI (based upon the 1.* API).

:ref:`FAQ <faq>` is compilation of frequently asked questions.

.. toctree::
  :hidden:
  :maxdepth: 1
  :caption: Reference

  installation
  faq
  old_cwipi/old_cwipi
  new_cwipi/new_cwipi
  client_server_cwipi/client_server_cwipi


