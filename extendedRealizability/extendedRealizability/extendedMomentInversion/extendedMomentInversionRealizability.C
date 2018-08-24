/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTAbiLITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "extendedMomentInversionRealizability.H"
#include "eigenSolver.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedMomentInversionRealizability, 0);
    defineRunTimeSelectionTable(extendedMomentInversionRealizability, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedMomentInversionRealizability::extendedMomentInversionRealizability
(
    const dictionary& dict,
    const label nMoments,
    const label nSecondaryNodes
)
:
    momentInverter_
    (
        univariateMomentInversion::New(dict.subDict("basicQuadrature"))
    ),
    nMoments_(nMoments),
    nPrimaryNodes_((nMoments_ - 1)/2),
    nSecondaryNodes_(nSecondaryNodes),
    primaryWeights_(nPrimaryNodes_, 0.0),
    primaryAbscissae_(nPrimaryNodes_, 0.0),
    sigma_(0.0),
    secondaryWeights_(nPrimaryNodes_, nSecondaryNodes_),
    secondaryAbscissae_(nPrimaryNodes_, nSecondaryNodes_),
    minMean_(dict.lookupOrDefault("minMean", 1.0e-8)),
    minVariance_(dict.lookupOrDefault("minVariance", 1.0e-8)),
    maxSigmaIter_(dict.lookupOrDefault<label>("maxSigmaIter", 1000)),
    momentsTol_(dict.lookupOrDefault("momentsTol", 1.0e-12)),
    sigmaMin_(dict.lookupOrDefault("sigmaMin", 1.0e-6)),
    sigmaTol_(dict.lookupOrDefault("sigmaTol", 1.0e-12)),
    targetFunctionTol_(dict.lookupOrDefault("targetFunctionTol", 1.0e-12)),
    foundUnrealizableSigma_(false),
    nullSigma_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedMomentInversionRealizability::~extendedMomentInversionRealizability()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::extendedMomentInversionRealizability::invert(const univariateMomentSet& moments)
{
    univariateMomentSet m(moments);

    reset();

    // Terminate execution if negative number density is encountered
    if (m[0] < 0.0)
    {
        FatalErrorInFunction
            << "The zero-order moment is negative." << nl
            << "    Moment set: " << m
            << abort(FatalError);
    }

    // Exclude cases where the zero-order moment is very small to avoid
    // problems in the inversion due to round-off error
    if (m[0] < SMALL)
    {
        sigma_ = 0.0;
        nullSigma_ = true;

        return;
    }

    label nRealizableMoments = m.nRealizableMoments();

    // If the moment set is on the boundary of the moment space, the
    // distribution will be reconstructed by a summation of Dirac delta,
    // and no attempt to use the extended quadrature method of moments is made.
    if (m.isOnMomentSpaceBoundary())
    {
        sigma_ = 0.0;
        nullSigma_ = true;
        momentInverter_().invert(m);

        secondaryQuadrature
        (
            momentInverter_().weights(),
            momentInverter_().abscissae()
        );

        return;
    }

    if (nRealizableMoments % 2 == 0)
    {
        // If the number of realizable moments is even, we apply the standard
        // QMOM directly to maximize the number of preserved moments.
        sigma_ = 0.0;
        nullSigma_ = true;
        momentInverter_().invert(m);

        secondaryQuadrature
        (
            momentInverter_().weights(),
            momentInverter_().abscissae()
        );
    }
    else
    {
        // Do not attempt the EQMOM reconstruction if mean or variance of the
        //  moment set are small to avoid numerical problems. These problems are
        // particularly acute in the calculation of the recurrence relationship
        // of the Jacobi orthogonal polynomials used for the beta kernel density
        // function.
        if (m[1]/m[0] < minMean_ || (m[2]/m[0] - sqr(m[1]/m[0])) < minVariance_)
        {
            sigma_ = 0.0;
            nullSigma_ = true;
            momentInverter_().invert(m);

            secondaryQuadrature
            (
                momentInverter_().weights(),
                momentInverter_().abscissae()
            );

            return;
        }

        // Resizing the moment set to avoid copying again
        m.resize(nRealizableMoments);

        // Local set of starred moments
        univariateMomentSet mStarNew(nRealizableMoments, m.support());
        univariateMomentSet mStar(nRealizableMoments, m.support());
        univariateMomentSet mStarHigh(nRealizableMoments, m.support());
        univariateMomentSet mStarMid(nRealizableMoments, m.support());

        // Compute target function for sigma = 0
        sigma_ = 0.0;
        momentsToMomentsStar(sigma_, m, mStar);
        // Realizability parameter (beta, canonical moments or zeta) that determines
        // whether the mStar is on moment space boundary
        scalarList realizParam = realizabilityParameter(mStar);
        label targetParameter = realizParam.size() - 1;
        scalar fZero = targetFunction(mStar, targetParameter);

        // Check if sigma = 0 leads to mStar on moment space boundary
        if (mStar.isOnMomentSpaceBoundary())
        {
            sigma_ = 0.0;
            nullSigma_ = true;
            momentInverter_().invert(m);

            secondaryQuadrature
            (
                momentInverter_().weights(),
                momentInverter_().abscissae()
            );

            return;
        }

        // Compute moments star for sigma = sigmaMax
        scalar sigMax = sigmaMax(m);
        scalar sigmaHigh = sigMax;
        momentsToMomentsStar(sigmaHigh, m, mStarHigh);

        // Check if sigmaHigh is solution
        if (mStarHigh.isOnMomentSpaceBoundary() && mStarHigh.nRealizableMoments() == nRealizableMoments)
        {
            sigma_ = sigmaHigh;
            momentInverter_().invert(mStarHigh);

            secondaryQuadrature
            (
                momentInverter_().weights(),
                momentInverter_().abscissae()
            );

            return;
        }

        // Apply Ridder's algorithm to find sigma
        for (label iter = 0; iter < maxSigmaIter_; iter++)
        {
            scalar sigmaMid = (sigma_ + sigmaHigh)/2.0;
            momentsToMomentsStar(sigmaMid, m, mStarMid);

            // Check if sigmaMid is solution
            if (mStarMid.isOnMomentSpaceBoundary() && mStarMid.nRealizableMoments() == nRealizableMoments)
            {
                sigma_ = sigmaMid;
                momentInverter_().invert(mStarMid);

                secondaryQuadrature
                (
                    momentInverter_().weights(),
                    momentInverter_().abscissae()
                );

                return;
            }

            // Update position of the first non-realizable parameter 
            // (beta, canonical moments or zeta)
            label realizabilityParameterI = firstNonRealizableParameter(mStarHigh);

            scalar fLow = targetFunction(mStar, realizabilityParameterI);
            scalar fHigh = targetFunction(mStarHigh, realizabilityParameterI);
            scalar fMid = targetFunction(mStarMid, realizabilityParameterI);

            scalar s = sqrt(sqr(fMid) - fLow*fHigh);

            if (s == 0.0)
            {
                FatalErrorInFunction
                    << "Singular value encountered searching for root." << nl
                    << "    Moment set = " << m << nl
                    << "    sigma = " << sigma_ << nl
                    << "    fLow = " << fLow << nl
                    << "    fMid = " << fMid << nl
                    << "    fHigh = " << fHigh
                    << abort(FatalError);
            }

            scalar sigmaNew = sigmaMid + (sigmaMid - sigma_)*sign(fLow - fHigh)*fMid/s;
            momentsToMomentsStar(sigmaNew, m, mStarNew);

            // Check if sigmaNew is solution
            if (mStarNew.isOnMomentSpaceBoundary() && mStarNew.nRealizableMoments() == nRealizableMoments)
            {
                sigma_ = sigmaNew;
                momentInverter_().invert(mStarNew);

                secondaryQuadrature
                (
                    momentInverter_().weights(),
                    momentInverter_().abscissae()
                );

                return;
            }

            // Refine search
            // Get the highest value value between sigma_, sigmaMid and sigmaNew 
            // that generates a realizable mStar
            if (mStarMid.isFullyRealizable())
            {
                if (mStarNew.isFullyRealizable() && sigmaNew > sigmaMid)
                {
                    sigma_ = sigmaNew;
                    mStar = mStarNew;
                }
                else
                {
                    sigma_ = sigmaMid;
                    mStar = mStarMid;
                }
            }
            // Get the lowest value value between sigmaHigh, sigmaMid and sigmaNew
            // that generates a non-realizable mStar
            if (!mStarMid.isFullyRealizable())
            {
                if (!mStarNew.isFullyRealizable() && sigmaNew < sigmaMid)
                {
                    sigmaHigh = sigmaNew;
                    mStarHigh = mStarNew;
                }
                else
                {
                    sigmaHigh = sigmaMid;
                    mStarHigh = mStarMid;
                }
            }

            // Check for convergence
            scalar dSigma = sigmaHigh - sigma_;
            scalar fTarget = targetFunction(mStar, targetParameter);
            if (mag(fTarget) <= targetFunctionTol_ * mag(fZero) || mag(dSigma) <= sigmaTol_ * sigMax)
            {
                // Root finding converged

                // If sigma_ is small, use QMOM
                if (mag(sigma_) < sigmaMin_)
                {
                    sigma_ = 0.0;
                    nullSigma_ = true;
                    momentInverter_().invert(m);

                    secondaryQuadrature
                    (
                        momentInverter_().weights(),
                        momentInverter_().abscissae()
                    );

                    return;
                }

                // Found a value of sigma that preserves all the moments
                momentInverter_().invert(mStar);

                secondaryQuadrature  // Secondary quadrature from mStar
                (
                    momentInverter_().weights(),
                    momentInverter_().abscissae()
                );

                return;

            }
        }

        FatalErrorInFunction
            << "Number of iterations exceeded." << nl
            << "    Max allowed iterations = " << maxSigmaIter_
            << abort(FatalError);
    }
}

void Foam::extendedMomentInversionRealizability::reset()
{
    foundUnrealizableSigma_ = false;
    nullSigma_ = false;

    forAll(primaryWeights_, pNodei)
    {
        primaryWeights_[pNodei] = 0.0;
        primaryAbscissae_[pNodei] = 0.0;

        for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
        {
            secondaryWeights_[pNodei][sNodei] = 0.0;
            secondaryAbscissae_[pNodei][sNodei] = 0.0;
        }
    }
}

void Foam::extendedMomentInversionRealizability::secondaryQuadrature
(
    const scalarList& pWeights,
    const scalarList& pAbscissae
)
{
    // Copy primary weights and abscissae
    forAll(pWeights, pNodei)
    {
        primaryWeights_[pNodei] = pWeights[pNodei];
        primaryAbscissae_[pNodei] = pAbscissae[pNodei];
    }

    if (!nullSigma_)
    {
        // Coefficients of the recurrence relation
        scalarDiagonalMatrix a(nSecondaryNodes_, 0.0);
        scalarDiagonalMatrix b(nSecondaryNodes_, 0.0);

        forAll(pWeights, pNodei)
        {
            // Compute coefficients of the recurrence relation
            recurrenceRelation(a, b, primaryAbscissae_[pNodei], sigma_);

            // Define the Jacobi matrix
            scalarSquareMatrix J(nSecondaryNodes_, 0.0);

            // Fill diagonal of Jacobi matrix
            forAll(a, ai)
            {
                J[ai][ai] = a[ai];
            }

            // Fill off-diagonal terms of the Jacobi matrix
            for (label bi = 0; bi < nSecondaryNodes_ - 1; bi++)
            {
                J[bi][bi + 1] = Foam::sqrt(b[bi + 1]);
                J[bi + 1][bi] = J[bi][bi + 1];
            }

            // Compute Gaussian quadrature used to find secondary quadrature
            eigenSolver JEig(J, true);

            const scalarDiagonalMatrix& JEigenvaluesRe(JEig.eigenvaluesRe());

            // Compute secondary weights before normalization and calculate sum
            for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei]
                    = sqr(JEig.eigenvectors()[0][sNodei]);

                secondaryAbscissae_[pNodei][sNodei] =
                    secondaryAbscissa(primaryAbscissae_[pNodei],
                        JEigenvaluesRe[sNodei], sigma_);
            }
        }

        // Set weights and abscissae of unused nodes to zero
        for (label pNodei = pWeights.size(); pNodei < nPrimaryNodes_; pNodei++)
        {
            primaryWeights_[pNodei] = 0.0;
            primaryAbscissae_[pNodei] = 0.0;

            for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei] = 0.0;
                secondaryAbscissae_[pNodei][sNodei] = 0.0;
            }
        }
    }
    else
    {
        // Manage case with null sigma to avoid redefining source terms
        forAll(pWeights, pNodei)
        {
            secondaryWeights_[pNodei][0] = 1.0;
            secondaryAbscissae_[pNodei][0] = primaryAbscissae_[pNodei];

            for (label sNodei = 1; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei] = 0.0;
                secondaryAbscissae_[pNodei][sNodei] = 0.0;
            }
        }

        // Set weights and abscissae of unused nodes to zero
        for (label pNodei = pWeights.size(); pNodei < nPrimaryNodes_; pNodei++)
        {
            primaryWeights_[pNodei] = 0.0;
            primaryAbscissae_[pNodei] = 0.0;

            for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei] = 0.0;
                secondaryAbscissae_[pNodei][sNodei] = 0.0;
            }
        }
    }
}

Foam::scalarList Foam::extendedMomentInversionRealizability::realizabilityParameter
(
    univariateMomentSet& momentsStar
)
{
    return momentsStar.zetas();
}

Foam::label Foam::extendedMomentInversionRealizability::firstNonRealizableParameter
(
    univariateMomentSet& momentsStar
)
{
    return momentsStar.negativeZeta() - 1; // The -1 converts from "mathematical index" to "array index"
}

Foam::scalar Foam::extendedMomentInversionRealizability::targetFunction
(
    univariateMomentSet& momentsStar,
    const label realizabilityParameterI
)
{
    scalarList zetas = momentsStar.zetas();
    return zetas[realizabilityParameterI];
}

// ************************************************************************* //
