.. _new_cwipi:

New CWIPI
#########

Released in 2023, this version relies on the parallel geometrical algortihm library ParaDiGM (Parallel Distributed Generalized
Mesh). ParaDiGM (LGPL) is developped by the CWIPI developpers team.

C API documentation
===================

.. doxygenfile:: cwp.h
   :project: cwipi

Python API documentation : pycwp
================================

.. currentmodule:: pycwp

.. automodule:: pycwp
   :imported-members:
   :members:
   :undoc-members:
   :show-inheritance:

Spatial interpolation methods
=============================

.. include:: spatial_interp.rst

Advanced functionalities
========================

For more precision, you may wish to adapt the spatial interpolation computed by CWIPI.
That is done by setting your own local interpolation function.
To write it you need local geometric data computed by CWIPI. Here is explained how to retrieve
within your own interpolation function.

Data getters for user interpolation functions in C
--------------------------------------------------

.. literalinclude:: ../../../../pattern/user_interp_function_c_pattern.txt
   :language: c

Data getters for user interpolation functions in Fortran
--------------------------------------------------------

.. literalinclude:: ../../../../pattern/user_interp_function_fortran_pattern.txt
   :language: fortran

Data getters for user interpolation functions in Python
-------------------------------------------------------

.. literalinclude:: ../../../../pattern/user_interp_function_python_pattern.txt
   :language: python


