.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.168586.svg
   :target: https://doi.org/10.5281/zenodo.168586
.. image:: https://travis-ci.org/scivision/zakharov.svg?branch=master
    :target: https://travis-ci.org/scivision/zakharov
.. image:: https://ci.appveyor.com/api/projects/status/n7cf0k2rwh5rggle?svg=true
    :target: https://ci.appveyor.com/project/scivision/zakharov



================================
Langmuir-1D Simulation
================================

Zakharov Simulation 

Originally authored by Hassanali Akbari as part of his PhD work.

Michael Hirsch converted to Fortran 2008, and works with `gfortran`, `ifort`, `flang` and other Fortran compilers on any operating system and computer.


:author: Hassanali Akbari, Michael Hirsch, Ph.D.

The procedure to use this program is as follows:

1. Build the Fortran code
2. run the simulation
3. plot with Matlab

.. contents::


Build
=====



Prereqs
-------

Linux
~~~~~
::

    apt install g++ cmake libboost-filesystem-dev libboost-program-options-dev
    
For `coarray` branch, additionally::

    apt install libcoarrays-dev
    
Mac
~~~
::

    brew install gcc boost


Compile
-------
::

    cd bin
    cmake ..
    make

Run Simulation
==============
arguments are:  output_directory simulation_end_time electron_beam_env(as many beams as you like)::

    ./zakhfort /tmp/test 1e-4 300


C++ usage
---------
Recommend using Fortran version instead::

    ./zakh --ev 300 -o /tmp/testcxx

--ev    beam energy
-o      output directory (will be created if it doesn't exist)

Plot Results
============
From GNU Octave or Matlab::

    Sim_v6_3_Linux(0, /tmp/test)
