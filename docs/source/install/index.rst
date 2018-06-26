.. _installation:

Installing OpenFAST
===================

Required software
-----------------

In order to build OpenFAST using CMake, one needs the following minimum set of packages installed:

- Fortran compiler (GNU compiler version above 4.6.0 or Intel compiler version above 11)

- C/C++ compiler

- GNU Make (version 3.81 or later)

- CMake (version 2.8.12 or later)

- HDF5 (Optional, required for C++ API)

Build instructions
------------------

The following pages provide instructions for building OpenFAST and/or its modules from source code.  
The developer team is moving towards a CMake-only approach that well supports Window Visual Studio users, but at this time we provide a separate build path for those users.

.. toctree::
   :maxdepth: 1

   get_openfast.rst
   install_cmake_linux.rst
   install_spack.rst
   install_cmake_cygwin.rst
   install_vs_windows.rst

