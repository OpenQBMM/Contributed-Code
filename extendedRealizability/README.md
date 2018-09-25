# Extended Moment Inversion Through Realizability of Degenerate Moments

Implementation of the algorithm by Pigou et al. (2018) for extended moment inversion 
using the realizability of degenerate moments.

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

1. Replace the folder *src/quadratureMethods/momentInversion/univariate/extended* 
for the one provided here.

2. In *src/quadratureMethods/momentInversion/Make/files*, replace :

```
univariate/extended/extendedMomentInversion/extendedMomentInversion.C
univariate/extended/extendedMomentInversion/newExtendedMomentInversion.C
univariate/extended/gamma/gammaEQMOM.C
univariate/extended/lognormal/lognormalEQMOM.C
univariate/extended/beta/betaEQMOM.C
```
for:
```
univariate/extended/extendedMomentInversions/extendedMomentInversion/extendedMomentInversion.C
univariate/extended/extendedMomentInversions/extendedMomentInversion/newExtendedMomentInversion.C
univariate/extended/extendedMomentInversions/highestMomentReconstruction/highestMomentReconstruction.C
univariate/extended/extendedMomentInversions/mStarRealizability/mStarRealizability.C

univariate/extended/kernelDensityFunctions/kernelDensityFunction/kernelDensityFunction.C
univariate/extended/kernelDensityFunctions/kernelDensityFunction/newKernelDensityFunction.C
univariate/extended/kernelDensityFunctions/gamma/gammaEQMOM.C
univariate/extended/kernelDensityFunctions/lognormal/lognormalEQMOM.C
univariate/extended/kernelDensityFunctions/beta/betaEQMOM.C
```

3. Apply the patch provided by *canonical_moments.patch*.

4. Recompile OpenQBMM.

## How to Use

The provided implementation should be backwards compatible with the existing
cases, with no syntax changes required. The previously implemented method 
*highestMomentReconstruction* is considered the default one.

To use the new inversion algorithm, add to the *extendedMomentInversion* dict in
*quadratureProperties*:
```
extendedMomentInversionMethod	mStarRealizability;
```

Three convergence criteria are available; *sigmaTolRel* and *targetFunctionTol*
follow the original paper proposal, and *sigmaTol* was added for easier 
comparison to the previous method. One can use, for example:
```
sigmaTol           0.0;
sigmaTolRel        1.0e-11;
targetFunctionTol  1.0e-11;
```

## Authors

* **Gabriel Gonçalves** (Núcleo Interdisciplinar de Dinâmica dos Fluidos - Universidade Federal do Rio de Janeiro) - *Initial work*
