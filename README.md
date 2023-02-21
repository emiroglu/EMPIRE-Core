# EMPIRE-Core
The core of the code for the EMPIRE-Multiphysics software developed in Chair of Structural Analysis at TUM

This repository contains the developments related to the following research:
- [Emiroglu, A. (2019)](https://mediatum.ub.tum.de/1473366). Multiphysics Simulation and CAD Integrated Shape Optimization in Fluid-Structure Interaction. (Doctoral dissertation, Technische Universität München)
- [Apostolatos, A. (2019)](https://mediatum.ub.tum.de/1453663). Isogeometric Analysis of Thin-Walled Structures on Multipatch Surfaces in Fluid-Structure Interaction. (Doctoral dissertation, Technische Universität München)
- [Wang, T. (2016)](https://mediatum.ub.tum.de/1281102). Development of Co-Simulation Environment and Mapping Algorithms. (Doctoral dissertation, Technische Universität München)
- [Sicklinger, S. (2014)](https://mediatum.ub.tum.de/1223319). Stabilized Co-Simulation of Coupled Problems Including Fields and Signals. (Doctoral dissertation, Technische Universität München)

# Building
```
export CC=icc
export CXX=icpc
cd build
cmake ..
```
Following make targets are available:
```
make          # compilation and linking)
make clean    # remove object files and executable including all folders)
make doc      # generates documentation) html main file is  /EMPEROR/doc/html/index.html
make cleandoc # removes documentation)
```
