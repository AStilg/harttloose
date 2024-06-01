# High-Aspect-Ratio T-Matrix Toolbox (HARTT): Loose file version
Whilst this work is under review for publication, it will be released under a AGPLv3.0 LICENSE (SEE LICENSE FILE).

## Getting started
HARTTloose is a package free repository of the HARTT functions. There is no installation for this software. This repository contains all of the needed files to produce the results showcased in **INSERT ARTICLE HERE**.

To run an example, simply enter the ``examples`` directory and run one of the examples, it will automatically add any paths you need to successfully run the code.

Personally, I suggest that the ``testing\run_test.m`` m-file be run first to see if there are some problems. If it fails the first time, run it twice more to see if it passes. It can fail in the tests checking wavefuncitons with a numerical calculation of curl because of numerical instability near poles. This is a far-future improvement for the toolbox unless it is something you would wish to look at!

## About the ``provided`` directory...
This directory contains several functions that can be found in other packages such as [ott](https://github.com/ilent2/ott/), in cases where you wish to use those functions instead you can interchange them without issue.

## Summary of toolbox functions
All main toolbox functions come with some sort of help message that describes what goes in and comes out of them. For more details use the MATLAB command line: ``>> help <function name>``.
- There are 4 coordinate transforms:
  - ``xyz2xietaphi``
  - ``rtp2xietaphi``
  - ``xietaphi2xyz``
  - ``xietaphi2rtp``
- There are 4 vector-field transforms that also call their matching coordinate transform:
  - ``xyzv2xietaphiv``
  - ``rtpv2xietaphiv``
  - ``xietaphiv2xyzv``
  - ``xietaphiv2rtpv``
- Three base wave functions:
  - ``spheroidalS1`` - Angular spheroidal wavefunctions of the first kind
  - ``spheroidalR1`` - Radial spheroidal wavefunctions of the first kind
  - ``spheroidalR2`` - Radial spheroidal wavefunctions of the second kind
- Two sets of vector wavefunctions
  - ``spheroidalvwf`` - Outputs outgoing, incoming, and regular vector wavefunctions
  - ``spheroidalvwf_farfield`` - Outputs only outgoing and incoming vector wavefunctions that are accurate in double precision only for $c\xi>10^{7.5}$.
- sT-Matrix functions that all start with ``stmatrix``:
  - ``stmatrix_spheroid_ebcm``
  - ``stmatrix_spheroid_pm``
  - ``stmatrix_cylinder_ebcm``
  - ``stmatrix_cylinder_pm``
- Beam shape coefficient generators that start with ``bsc`` and end with ``spheroidal``
  - ``bsc_plane_spheroidal``
  - ``bsc_bessel_spheroidal``
- Mode transformation functions:
  - ``spheroidal_u_coefficients`` - The fundamental routine that allows the toolbox to function. It is responsible for solving the spheroidal eigenvalue problem and generating the basis weights.
  - ``spheroidal_expansion`` - Generates the vector mode transformation between spherical-and-spheroidal wavefunctions.
  - ``spherical_to_spheroidal`` - Enables *ott* T-matrix and beam shape coefficients (the data) to be converted to spheroidal coefficients.
  - ``spheroidal_to_spherical`` - Transforms spheroidal sT-matrix and beam shape coefficients to be converted to data for import into *ott*.
- Other needed functions:
  - ``aspect_ratio_to_conk`` - Estimates the interfocal radius needed for an expansion in spheroidal wavefunctions.
  - ``change_nmax_spheroidal``
  - ``mode_couplings_spheroidal`` - For sT-matrix that use rotational and mirror symmetry reduction
  - ``ndotxietaphicross`` - For EBCM
  - ``neg_spheroidal_continued_fractions`` - used by ``spheroidal_u_coefficients``
  - ``spheroidal_continued_fractions`` - used by ``spheroidal_u_coefficients``
  - ``spheroidal_scale_factors`` - (pretty optional though is used somewhere).
