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

#include "extendedMomentInversion.H"
#include "eigenSolver.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedMomentInversion, 0);
    defineRunTimeSelectionTable(extendedMomentInversion, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedMomentInversion::extendedMomentInversion
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
    kernel_
    (
        kernelDensityFunction::New(dict)
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
    sigmaMin_(dict.lookupOrDefault("sigmaMin", 1.0e-6)),
    foundUnrealizableSigma_(false),
    nullSigma_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedMomentInversion::~extendedMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::extendedMomentInversion::reset()
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

void Foam::extendedMomentInversion::secondaryQuadrature
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
            kernel_().recurrenceRelation(a, b, primaryAbscissae_[pNodei], sigma_);

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
                    kernel_().secondaryAbscissa(primaryAbscissae_[pNodei],
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

void Foam::extendedMomentInversion::invertSingular(univariateMomentSet& m)
{
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

    if (m.nRealizableMoments() % 2 == 0)
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
	}
}

Foam::tmp<Foam::scalarField> Foam::extendedMomentInversion::f(const scalarField& x) const
{
    tmp<scalarField> tmpY
    (
        new scalarField(x.size(), 0.0)
    );
    scalarField& y = tmpY.ref();

    for (label pNodei = 0; pNodei < nPrimaryNodes_; pNodei++)
    {
        y += kernel_().f(x, primaryAbscissae_[pNodei], sigma_)
           *primaryWeights_[pNodei];
    }

    return tmpY;
}

// ************************************************************************* //
