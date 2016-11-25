.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.168586.svg
   :target: https://doi.org/10.5281/zenodo.168586
.. image:: https://travis-ci.org/scienceopen/zakharov.svg?branch=master
    :target: https://travis-ci.org/scienceopen/zakharov

================================
Langmuir-1D Simulation
================================

Zakharov Simulation 

Originally authored by Hassanali Akbari as part of his PhD work.

:author: Hassanali Akbari, Michael Hirsch

.. contents::

Prereq
======

Linux / BSD / Windows::

    sudo apt-get install g++ cmake libboost-filesystem-dev libboost-program-options-dev
    
Mac::

    brew install boost


Setup
=====
::

    cd bin
    cmake ..
    make

Usage
=====
::

    ./zahk
    
--ev    beam energy
-o      output directory (will be created if it doesn't exist)
