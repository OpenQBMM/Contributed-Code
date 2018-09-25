/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mStarRealizability

Description
    Class for the optimized algorithm for the extended quadrature method of 
	moments proposed by Pigou. The value of sigma is chosen such as to 
	have the transformed moments mStar be on the boundary of the realizability
	space.

    References
    \verbatim
        "New developments of the Extended Quadrature Method of Moments to 
        solve Population Balance Equations"
        M Pigou, J Morchain, P Fede, M-I Penet, G Laronze
        Journal of Computational Physics
        Volume 365, Pages 243-268, 2018
    \endverbatim

SourceFiles
    mStarRealizability.C

\*---------------------------------------------------------------------------*/

#ifndef mStarRealizability_H
#define mStarRealizability_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "extendedMomentInversion.H"
#include "runTimeSelectionTables.H"
#include "univariateMomentSet.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mStarRealizability Declaration
\*---------------------------------------------------------------------------*/

class mStarRealizability
:
    public extendedMomentInversion
{
	// Private data

        //- Tolerance for the change in the sigma parameter
        scalar sigmaTol_;

        //- Relative tolerance for the change in sigma parameter
        scalar sigmaTolRel_;

        //- Tolerance for the target function
        scalar targetFunctionTol_;


    // Private member functions

        //- Return list of parameters for moment set realizability
        //  (beta, canonical moments or zeta)
        scalarList realizabilityParameter
        (
            univariateMomentSet& momentsStar
        );

        //- Return index of first invalid realizability parameter in moment set
        label firstNonRealizableParameter
        (
            univariateMomentSet& momentsStar
        );

        //- Compute target function whose root is sigma_
        scalar targetFunction
        (
            univariateMomentSet& momentsStar,
            const label realizabilityParameterI
        );

        //- Disallow default bitwise copy construct
        mStarRealizability(const mStarRealizability&);

        //- Disallow default bitwise assignment
        void operator=(const mStarRealizability&);


public:

    //- Runtime type information
    TypeName("mStarRealizability");


    // Constructors

        // Construct from dictionary and label
        mStarRealizability
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
    virtual ~mStarRealizability();


    // Member Functions

        //- Invert moments to find weight, abscissae and sigma
        void invert(const univariateMomentSet& moments);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
