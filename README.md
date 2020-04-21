# TMM-structure-functions
Structural analysis functions required for implementing Transfer Matrix Method with applied axial loading.

> Please cite as: Chatterjee, P., & Bryant, M. (2019). Analysis of Tension-Tunable Clamped-Clamped Piezoelectric Beams for Harvesting Energy from Wind and Vibration, Journal of Intelligent Material Systems and Structures. https://doi.org/10.1177/1045389X19862390

This repository is part of many others used to successully run the primary MATLAB code used during graduate studies by [Dr. Punnag Chatterjee](https://sites.google.com/view/punnagchatterjee/home).

The structure of the MAIN (driver) code looks like this:

**Step 1: Structural**

a. reading geometry

b. reading material

c. structural simulation options

d. geometry discretization

e. structural module - evaluates mode shape and natural frequencies for bending/torsion

**Step 2: Electromechanical**

electromechanical module - evaluates the piezo-mechanical coupling parameters

**Step 3: Aerodynamic**

a. aerodynamic simulation options

b. creating aerodynamic matrices

**Step 4: Coupled - solution**

a. ode solver 

**Step 5: Post-processing**

## How to contribute?

Please let me know if you find bugs or wish to contribute.

## Copyright and License

All content is under Creative Commons Attribution [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode.txt), and all [code is under BSD-3 clause](https://github.com/engineersCode/EngComp/blob/master/LICENSE). 

I am fine if you re-use the content in any way!

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
