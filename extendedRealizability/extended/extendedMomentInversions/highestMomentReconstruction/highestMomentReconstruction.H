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
    Foam::highestMomentReconstruction

Description
    Class for the classical algorithm for the extended quadrature method of 
	moments. The value of sigma is chosen such as to reproduce as closely
	as possible the highest order moment.

    References
    \verbatim
        "An extended quadrature method of moments for population balance
        equations"
        C Yuan, F Laurent, R O Fox
        Journal of Aerosol Science
        Volume 51, Pages 1-23, 2012
    \endverbatim
	
    \verbatim
        "The theory of canonical moments with applications in Statistics,
        Probability and Analysis"
        H Dette, W J Studden
        Wiley & Sons, 1997
    \endverbatim

SourceFiles
    highestMomentReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef highestMomentReconstruction_H
#define highestMomentReconstruction_H

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
                     Class highestMomentReconstruction Declaration
\*---------------------------------------------------------------------------*/

class highestMomentReconstruction
:
    public extendedMomentInversion
{
	// Private data
	
        //- Tolerance on moment conservation
        scalar momentsTol_;

        //- Tolerance for the change in the sigma parameter
        scalar sigmaTol_;

        //- Tolerance for the target function
        scalar targetFunctionTol_;
	
    // Private member functions

        //- Attempt to find a min or max of the target function on an interval
        scalar minimizeTargetFunction
        (
            scalar sigmaLow,
            scalar sigmaHigh,
            const univariateMomentSet& moments,
            univariateMomentSet& momentsStar
        );

        //- Compute target function whose root is sigma_
        scalar targetFunction
        (
            scalar sigma,
            const univariateMomentSet& moments,
            univariateMomentSet& momentsStar
        );

        //- Compute the normalized moment error
        scalar normalizedMomentError
        (
            scalar sigma,
            const univariateMomentSet& moments,
            univariateMomentSet& momentsStar
        );

        //- Disallow default bitwise copy construct
        highestMomentReconstruction(const highestMomentReconstruction&);

        //- Disallow default bitwise assignment
        void operator=(const highestMomentReconstruction&);


public:

    //- Runtime type information
    TypeName("highestMomentReconstruction");


    // Constructors

        // Construct from dictionary and label
        highestMomentReconstruction
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
    virtual ~highestMomentReconstruction();


    // Member Functions

        //- Invert moments to find weight, abscissae and sigma
        void invert(const univariateMomentSet& moments);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
