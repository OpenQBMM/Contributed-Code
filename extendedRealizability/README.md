# Extended Moment Inversion Through Realizability of Transformed Moments

Implementation of the algorithm by Pigou et al. (2018) for extended moment inversion 
using the realizability of transformed moments.

**Reference**

Maxime Pigou, Jérôme Morchain, Pascal Fede, Marie-Isabelle Penet, Geoffrey Laronze, 
New developments of the Extended Quadrature Method of Moments to solve Population Balance Equations, 
Journal of Computational Physics, 
Volume 365, 2018, Pages 243-268, ISSN 0021-9991, doi.org/10.1016/j.jcp.2018.03.027. 

## Prerequisites

The following libraries are required:

* OpenFOAM-6
* OpenQBMM (development version, tested with 8fbff4c)

## Installing

The library and test script may be installed through:

```
./Allwmake
```

## Test script

The test script performs the inversion of a given set of moments. It is executed with:

```
cd Test-extendedRealizability
Test-ExtendedMomentInversionRealizability
```

To change the moment set, edit the file *Test-ExtendedMomentInversion.C* and recompile.

## Authors

* **Gabriel Gonçalves** (Núcleo Interdisciplinar de Dinâmica dos Fluidos - Universidade Federal do Rio de Janeiro) - *Initial work*
