# TMM-structure-functions
Structural analysis functions required for implementing Transfer Matrix Method with applied axial loading.

> Please cite as: Chatterjee, P., & Bryant, M. (2019). Analysis of Tension-Tunable Clamped-Clamped Piezoelectric Beams for Harvesting Energy from Wind and Vibration, Journal of Intelligent Material Systems and Structures. https://doi.org/10.1177/1045389X19862390

This repository contains the MATLAB functions to run a structural modal analysis with bending-torsional degrees of freedom [Dr. Punnag Chatterjee](https://sites.google.com/view/punnagchatterjee/home).

**Why TMM?**
In short, TMM provides a very methodical way to include geometrical and material discontinuities in a structure. As opposed to FEA, this is an exact method. It means, you can extract the exact symbolic/mathematical mode shape functions of the structure.
Also, this is faster compared to FEA as the matrix size is limited to 4x4 (or 6x6), unlike FEA where the more nodes result in larger and mlarger matrices that are solved. In TMM, the idea is to cleverly dicrstize a structure into multiple segments such that geometrical and material properties are homogenous within the segment. One limitation that immediately comes to my mind is it is more suited for analyzing simple structures. 

## Short instruction - to get started
run the code "MAIN_bending_torsion_modal_analysis.m"

## Code Structure

**Structural Module**

a. reading geometry

b. reading material

c. structural simulation options

d. geometry discretization

e. structural module - evaluates mode shape and natural frequencies for bending/torsion

## Coming up
1. Detailed instructions to change the simulation for different geometries
2. Electromechanical module
3. Aerodynamic module
4. Couplings to connect all the modules together

## How to contribute?

Please let me know if you find bugs or wish to contribute.

## Copyright and License

All content is under Creative Commons Attribution [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode.txt), and all [code is under BSD-3 clause](https://github.com/engineersCode/EngComp/blob/master/LICENSE). 

I am fine if you re-use the content in any way!

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)


