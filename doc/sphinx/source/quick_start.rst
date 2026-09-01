.. _quick_start:

Quick Start
###########

Here is shown how to couple two codes using the 1.x version of CWIPI in the different
programming languages in which CWIPI is available. Both codes have the same basic
polygonal mesh (see figure). code1 sends a field on the nodes to code2 and the unmapped points are checked.

.. image:: ./images/mesh.png
   :scale: 50%

CWIPI C
-------

Code
~~~~

.. literalinclude:: ../../../tests/tutorial/c_new_api_polygon_sol.c
   :language: c

Compilation
~~~~~~~~~~~

.. code-block:: sh

  export CWIPI_INSTALL_DIR=<path>/<where>/<cwipi>/<is>/<installed>
  mpicc -I$CWIPI_INSTALL_DIR/include -L$CWIPI_INSTALL_DIR/lib -o <exec> <file>.c -lcwp

Execution
~~~~~~~~~

.. code-block:: sh

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CWIPI_INSTALL_DIR/lib/
  mpirun -np 2 ./<exec>


CWIPI Fortran
-------------

Code
~~~~

.. literalinclude:: ../../../tests/tutorial/fortran_new_api_polygon_sol.F90
   :language: fortran

Compilation
~~~~~~~~~~~

.. code-block:: sh

  export CWIPI_INSTALL_DIR=<path>/<where>/<cwipi>/<is>/<installed>
  mpif90 -I$CWIPI_INSTALL_DIR/include -L$CWIPI_INSTALL_DIR/lib -o <exec> <file>.f90 -lcwpf -lcwp -cpp

Execution
~~~~~~~~~

.. code-block:: sh

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CWIPI_INSTALL_DIR/lib/
  mpirun -np 2 ./<exec>


CWIPI Python
------------

Code
~~~~

.. literalinclude:: ../../../tests/tutorial/python_new_api_polygon_sol.py
   :language: python


Execution
~~~~~~~~~

.. code-block:: sh

  export CWIPI_INSTALL_DIR=<path>/<where>/<cwipi>/<is>/<installed>
  export PYTHONPATH=$CWIPI_INSTALL_DIR/lib/python<version>/site-packages:$PYTHONPATH
  mpirun -np 2 python <file>.py