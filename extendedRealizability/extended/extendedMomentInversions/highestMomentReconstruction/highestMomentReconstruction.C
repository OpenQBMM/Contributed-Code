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

#include "highestMomentReconstruction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(highestMomentReconstruction, 0);

    addToRunTimeSelectionTable
    (
        extendedMomentInversion,
        highestMomentReconstruction,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highestMomentReconstruction::highestMomentReconstruction
(
    const dictionary& dict,
    const label nMoments,
    const label nSecondaryNodes
)
:
    extendedMomentInversion(dict, nMoments, nSecondaryNodes),
    momentsTol_(dict.lookupOrDefault("momentsTol", 1.0e-12)),
    sigmaTol_(dict.lookupOrDefault("sigmaTol", 1.0e-12)),
    targetFunctionTol_(dict.lookupOrDefault("targetFunctionTol", 1.0e-12))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highestMomentReconstruction::~highestMomentReconstruction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::highestMomentReconstruction::invert(const univariateMomentSet& moments)
{
    univariateMomentSet m(moments);

    reset();
	
	invertSingular(m);
	if (nullSigma_) return;
	
	label nRealizableMoments = m.nRealizableMoments();
	
	// Resizing the moment set to avoid copying again
	m.resize(nRealizableMoments);

	// Local set of starred moments
	univariateMomentSet mStar(nRealizableMoments, m.support());

	// Compute target function for sigma = 0
	scalar sigmaLow = 0.0;
	scalar fLow = targetFunction(sigmaLow, m, mStar);
	sigma_ = sigmaLow;

	// Check if sigma = 0 is root
	if (mag(fLow) <= targetFunctionTol_)
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

	// Compute target function for sigma = sigmaMax
	scalar sigMax = kernel_().sigmaMax(m);
	scalar sigmaHigh = sigMax;
	scalar fHigh = targetFunction(sigmaHigh, m, mStar);

	// This should not happen with the new algorithm
	if (fLow*fHigh > 0)
	{
		// Root not found. Minimize target function in [0, sigma_]
		sigma_ = minimizeTargetFunction(0, sigmaHigh, m, mStar);

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

		targetFunction(sigma_, m, mStar);
		secondaryQuadrature  // secondary quadrature from mStar
		(
			momentInverter_().weights(),
			momentInverter_().abscissae()
		);

		return;
	}

	// Apply Ridder's algorithm to find sigma
	for (label iter = 0; iter < maxSigmaIter_; iter++)
	{
		scalar fMid, sigmaMid;

		sigmaMid = (sigmaLow + sigmaHigh)/2.0;
		fMid = targetFunction(sigmaMid, m, mStar);

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

		sigma_ = sigmaMid + (sigmaMid - sigmaLow)*sign(fLow - fHigh)*fMid/s;

		kernel_().momentsToMomentsStar(sigma_, m, mStar);

		scalar fNew = targetFunction(sigma_, m, mStar);
		scalar dSigma = (sigmaHigh - sigmaLow)/2.0;

		// Check for convergence
		if (mag(fNew) <= targetFunctionTol_ || mag(dSigma) <= sigmaTol_)
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

			scalar momentError = normalizedMomentError(sigma_, m, mStar);

			if
			(
				momentError < momentsTol_
			)
			{
				// Found a value of sigma that preserves all the moments
				secondaryQuadrature  // Secondary quadrature from mStar
				(
					momentInverter_().weights(),
					momentInverter_().abscissae()
				);

				return;
			}
			else
			{
				// Root not found. Minimize target function in [0, sigma_]
				sigma_ = minimizeTargetFunction(0, sigma_, m, mStar);

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

				targetFunction(sigma_, m, mStar);

				secondaryQuadrature // Secondary quadrature from  mStar
				(
					momentInverter_().weights(),
					momentInverter_().abscissae()
				);

				return;
			}
		}
		else
		{
			// Root finding did not converge. Refine search.

			if (fNew*fMid < 0 && sigma_ < sigmaMid)
			{
				sigmaLow = sigma_;
				fLow = fNew;
				sigmaHigh = sigmaMid;
				fHigh = fMid;
			}
			else if (fNew*fMid < 0 && sigma_ > sigmaMid)
			{
				sigmaLow = sigmaMid;
				fLow = fMid;
				sigmaHigh = sigma_;
				fHigh = fNew;
			}
			else if (fNew*fLow < 0)
			{
				sigmaHigh = sigma_;
				fHigh = fNew;
			}
			else if (fNew*fHigh < 0)
			{
				sigmaLow = sigma_;
				fLow = fNew;
			}
		}
	}

	FatalErrorInFunction
		<< "Number of iterations exceeded." << nl
		<< "    Max allowed iterations = " << maxSigmaIter_
		<< abort(FatalError);
}

Foam::scalar Foam::highestMomentReconstruction::minimizeTargetFunction
(
    scalar sigmaLow,
    scalar sigmaHigh,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    const scalar goldenRatio = (sqrt(5.0) - 1.0)/2.0;

    scalar a = sigmaLow;
    scalar b = sigmaHigh;
    scalar x = b - goldenRatio*(b - a);
    scalar y = a + goldenRatio*(b - a);

    label iter = 0;

    while (mag (x - y) > sigmaTol_ && iter < maxSigmaIter_)
    {
        // Square the target function to find closest value to zero,
        // independently from the sign of the function
        scalar fx = sqr(targetFunction(x, moments, momentsStar));
        scalar fy = sqr(targetFunction(y, moments, momentsStar));

        if (fx < fy)
        {
            b = y;
            y = x;
            x = b - goldenRatio*(b - a);
        }
        else
        {
            a = x;
            x = y;
            y = a + goldenRatio*(b - a);
        }

        iter++;
    }

    if (iter > maxSigmaIter_)
    {
        FatalErrorInFunction
            << "Number of iterations exceeded." << nl
            << "    Max allowed iterations = " << maxSigmaIter_
            << abort(FatalError);
    }

    return (a + b)/2.0;
}

Foam::scalar Foam::highestMomentReconstruction::normalizedMomentError
(
    scalar sigma,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    scalar norm = 0.0;

    targetFunction(sigma, moments, momentsStar);

    univariateMomentSet approximatedMoments(moments.size(), moments.support());

    kernel_().momentsStarToMoments(sigma, approximatedMoments, momentsStar);

    for (label momenti = 0; momenti < moments.size(); momenti++)
    {
        norm += mag(1.0 - approximatedMoments[momenti]/moments[momenti]);
    }

    return sqrt(norm);
}

Foam::scalar Foam::highestMomentReconstruction::targetFunction
(
    scalar sigma,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    kernel_().momentsToMomentsStar(sigma, moments, momentsStar);

    momentInverter_().invert(momentsStar);
    momentsStar.update
    (
        momentInverter_().weights(),
        momentInverter_().abscissae()
    );

    scalar highestMoment = moments.last();

    return (highestMoment - kernel_().m2N(sigma, momentsStar))/highestMoment;
}

// ************************************************************************* //
