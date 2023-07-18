.. _new_cwipi:

New CWIPI
#########

Released in 2023, this version relies on the parallel computational geometry library **ParaDiGM** which is developed by the same team as CWIPI.

General concepts
================

.. include:: concepts.rst

Spatial interpolation methods
=============================

.. include:: spatial_interp.rst

Advanced features
=================


User-defined spatial interpolation
----------------------------------

Users can perform customized spatial interpolation based on the geometric mapping computed by CWIPI.
To do so, one needs to provide CWIPI with a pointer to a user-written interpolation function (via ``CWP_Field_interp_function_set`` in C).
Each Field object can have a specific user-defined spatial interpolation.
The following sections show how to retrieve all the information required to write such a function.

Data getters for user interpolation functions in C
""""""""""""""""""""""""""""""""""""""""""""""""""

.. literalinclude:: ../../../../pattern/user_interp_function_c_pattern.txt
   :language: c

Data getters for user interpolation functions in Fortran
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. literalinclude:: ../../../../pattern/user_interp_function_fortran_pattern.txt
   :language: fortran

Data getters for user interpolation functions in Python
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. literalinclude:: ../../../../pattern/user_interp_function_python_pattern.txt
   :language: python


Additional exchange features
----------------------------

CWIPI v.1 offers two additional means to exchange data between the coupled codes.
As opposed to fields, these data need not be related to the interface meshes.

.. include:: global_data.rst

.. include:: part_data.rst



C API documentation
===================

.. doxygenfile:: cwp.h
   :project: cwipi

Fortran API documentation
=========================

*à remplir*

.. doxygenfile:: cwp_f.f90
   :project: cwipi

Python API documentation : pycwp
================================

.. currentmodule:: pycwp

.. automodule:: pycwp
   :imported-members:
   :members:
   :undoc-members:
   :show-inheritance:


