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

#include "mStarRealizability.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mStarRealizability, 0);

    addToRunTimeSelectionTable
    (
        extendedMomentInversion,
        mStarRealizability,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mStarRealizability::mStarRealizability
(
    const dictionary& dict,
    const label nMoments,
    const label nSecondaryNodes
)
:
    extendedMomentInversion(dict, nMoments, nSecondaryNodes),
    sigmaTol_(dict.lookupOrDefault("sigmaTol", 1e-12)),
    sigmaTolRel_(dict.lookupOrDefault("sigmaTolRel", 1.0e-10)),
    targetFunctionTol_(dict.lookupOrDefault("targetFunctionTol", 1.0e-10))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mStarRealizability::~mStarRealizability()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::mStarRealizability::invert(const univariateMomentSet& moments)
{
    univariateMomentSet m(moments);

    reset();

	invertSingular(m);
	if (nullSigma_) return;
	
	label nRealizableMoments = m.nRealizableMoments();

	// Resizing the moment set to avoid copying again
	m.resize(nRealizableMoments);

	// Local set of starred moments
	univariateMomentSet mStarNew(nRealizableMoments, m.support());
	univariateMomentSet mStar(nRealizableMoments, m.support());
	univariateMomentSet mStarHigh(nRealizableMoments, m.support());
	univariateMomentSet mStarMid(nRealizableMoments, m.support());

	// Compute target function for sigma = 0
	sigma_ = 0.0;
	kernel_().momentsToMomentsStar(sigma_, m, mStar);
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
	scalar sigMax = kernel_().sigmaMax(m);
	scalar sigmaHigh = sigMax;
	kernel_().momentsToMomentsStar(sigmaHigh, m, mStarHigh);

  scalar fTol = targetFunctionTol_ * mag(fZero);
  scalar sTol = max(sigmaTolRel_ * sigMax, sigmaTol_);

	// Check if sigmaHigh is solution
	if (mStarHigh.isOnMomentSpaceBoundary() && 
      mStarHigh.nRealizableMoments() == nRealizableMoments)
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
		kernel_().momentsToMomentsStar(sigmaMid, m, mStarMid);

		// Check if sigmaMid is solution
		if (mStarMid.isOnMomentSpaceBoundary() && 
        mStarMid.nRealizableMoments() == nRealizableMoments)
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

    if (m.support() == "01" &&
        realizabilityParameterI < m.nRealizableMoments() &&
        fHigh > 1.0)
    {
      fLow = 1.0 - fLow;
      fHigh = 1.0 - fHigh;
      fMid = 1.0 - fMid;
    }

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
		kernel_().momentsToMomentsStar(sigmaNew, m, mStarNew);

		// Check if sigmaNew is solution
		if (mStarNew.isOnMomentSpaceBoundary() && 
        mStarNew.nRealizableMoments() == nRealizableMoments)
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
		if (mag(fTarget) <= fTol || mag(dSigma) <= sTol)
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

Foam::scalarList Foam::mStarRealizability::realizabilityParameter
(
    univariateMomentSet& momentsStar
)
{
	scalarList parameterList;
  if (momentsStar.support() == "RPlus")
  {
    parameterList = momentsStar.zetas();
  } 
  else if (momentsStar.support() == "01")
  {
    parameterList = momentsStar.canonicalMoments();
  }
  else if (momentsStar.support() == "R")
  {
    NotImplemented;
  }
  else
  {
    FatalErrorInFunction
        << "The specified support is invalid." << nl
        << "    Valid supports are: R, RPlus and 01."
        << abort(FatalError);
  }
  return parameterList;
}

Foam::label Foam::mStarRealizability::firstNonRealizableParameter
(
    univariateMomentSet& momentsStar
)
{
	label idx = 0;
  if (momentsStar.support() == "RPlus")
  {
    idx = momentsStar.negativeZeta() - 1;
  } 
  else if (momentsStar.support() == "01")
  {
    idx = momentsStar.nRealizableMoments() - 1;
  }
  else if (momentsStar.support() == "R")
  {
    NotImplemented;
  }
  else
  {
    FatalErrorInFunction
        << "The specified support is invalid." << nl
        << "    Valid supports are: R, RPlus and 01."
        << abort(FatalError);
  }
  return idx;
}

Foam::scalar Foam::mStarRealizability::targetFunction
(
    univariateMomentSet& momentsStar,
    const label realizabilityParameterI
)
{
  scalarList realizParams = realizabilityParameter(momentsStar);
  return realizParams[realizabilityParameterI];
}

// ************************************************************************* //
