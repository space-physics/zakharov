.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.168586.svg
   :target: https://doi.org/10.5281/zenodo.168586
.. image:: https://travis-ci.org/scivision/zakharov.svg?branch=master
    :target: https://travis-ci.org/scivision/zakharov

================================
Langmuir-1D Simulation
================================

Zakharov Simulation 

Originally authored by Hassanali Akbari as part of his PhD work.

Michael Hirsch converted to Fortran 2008, and should work with `gfortran`, `ifort`, `flang` and other Fortran compilers on any operating system and computer.

:author: Hassanali Akbari, Michael Hirsch, Ph.D.

.. contents::

Prereq
======

Linux / BSD / Windows Subsystem for Linux::

    apt install g++ cmake libboost-filesystem-dev libboost-program-options-dev
    
Mac::

    brew install gcc boost


Setup
=====
::

    cd bin
    cmake ..
    make

1-D ZakharovUsage
=================
arguments are:  output_directory simulation_end_time electron_beam_env(as many beams as you like)::

    ./zakhfort /tmp/testfort 100e-3 300


    

C++ usage
=========
Recommend using Fortran version instead::

    ./zakh --ev 300 -o /tmp/testcxx

--ev    beam energy
-o      output directory (will be created if it doesn't exist)
