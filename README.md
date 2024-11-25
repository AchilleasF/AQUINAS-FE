# AQUINAS     `0.10.0`
__A Matlab based program for the analysis of axisymmetric shells__

## Authors

[Achilleas Filippidis](mailto:a.filippidis@imperial.ac.uk?subject=[GitHub]%20AQUINAS%20Query) and [Adam J. Sadowski](mailto:a.sadowski@imperial.ac.uk?subject=[GitHub]%20AQUINAS%20Query)

Department of Civil and Environmental Engineering, Imperial College London, UK

## Abstract

AQUINAS is a specialised MATLAB finite element (FE) toolbox for the nonlinear analysis of axisymmetric thin-walled shell systems, conceived by the Authors to recover a research functionality once widely available but now largely lost in the age of supposedly 'general' 3D FE software. AQUINAS exploits the powerful built-in sparse linear algebra routines of MATLAB with efficient matrix assembly in pre-compiled and parallelised C++ code via the MEX functionality, as well as MATLAB's extensive visualisation and quality-of-life features. In addition to making nonlinear bifurcation buckling calculations of axisymmetric shells more accessible to a wider research community, the toolbox is designed to facilitate the ongoing active development of the structural Eurocode EN 1993-1-6 governing the design of metal civil engineering shells through its capability to automate the accurate computation of nominal resistance capacity curves for reference axisymmetric shell systems.

## Requirements

### Usage
AQUINAS has been developed and tested on a MATLAB R2022a coding environment. Earlier versions of MATLAB may not be supported.
Only the base MATLAB installation is required, with no dependency on any of MATLAB's toolboxes.
Pre-compiled binaries `.mexw64` and `.mexa64` are included for Windows and Linux operating systems.

### Compilation/Development
MinGW-w64 6.3 (for OpenMP 4.5 support) is required for compilation of the C++ code into MATLAB's `MEX` binaries. This can be installed through MATLAB's Add-Ons manager.
Maple (2019 or later) can be used for automatic generation of the scripts (both for MATLAB and C++ assembly versions) responsible for computation of the stiffness coefficients. This is optional and can be done manually if a Maple license is unavailable.

## Contact
a.filippidis@imperial.ac.uk  |  a.sadowski@imperial.ac.uk
