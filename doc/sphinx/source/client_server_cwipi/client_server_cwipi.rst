.. _client_server_cwipi:

Client-Server
#############

Client-Server CWIPI is a degenerate mode of the new version of CWIPI in terms of performance.
It has been developped to be able to do a coupling between parallel codes with a different version of MPI.
The API's of this mode and the new CWIPI version are broadly similar.

There is a Python interface for this mode. A Fortran might be developped upon request.

Server
------

The server ``cwp_server`` has to be launch on the same number of MPI ranks as the client coupling.

.. code:: bash


