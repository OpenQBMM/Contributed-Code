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
    Foam::extendedMomentInversion

Description
    Abstract class to allow run-time selection of the algorithm to be used
    for the extended quadrature method of moments.

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

    \verbatim
        "Orthogonal Polynomials: Computation and Approximation"
        W Gautschi
        Oxford University Press, 2004
    \endverbatim

SourceFiles
    extendedMomentInversion.C
    newExtendedMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef extendedMomentInversion_H
#define extendedMomentInversion_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "kernelDensityFunction.H"
#include "runTimeSelectionTables.H"
#include "univariateMomentSet.H"
#include "univariateMomentInversion.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class extendedMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class extendedMomentInversion
{
    // Private member functions

        //- Disallow default bitwise copy construct
        extendedMomentInversion(const extendedMomentInversion&);

        //- Disallow default bitwise assignment
        void operator=(const extendedMomentInversion&);


protected:

    // Protected data

        //- Univariate moment inversion method
        autoPtr<univariateMomentInversion> momentInverter_;

        //- Kernel density function
        autoPtr<kernelDensityFunction> kernel_;

        //- Number of moments
        label nMoments_;

        //- Number of primary quadrature nodes
        label nPrimaryNodes_;

        //- Number of secondary quadrature nodes
        label nSecondaryNodes_;

        //- Primary quadrature weights
        scalarList primaryWeights_;

        //- Primary quadrature abscissae
        scalarList primaryAbscissae_;

        //- Parameter sigma of the kernel density function
        scalar sigma_;

        //- Secondary quadrature weights
        scalarRectangularMatrix secondaryWeights_;

        //- Secondary quadrature abscissae
        scalarRectangularMatrix secondaryAbscissae_;

        //- Minimum mean to perform EQMOM reconstruction
        scalar minMean_;

        //- Minimum variance to perform EQMOM reconstruction
        scalar minVariance_;

        //- Maximum number of iterations allowed to find sigma_
        label maxSigmaIter_;

        //- Minimum allowed value of the parameter sigma of the kernel density
        //  function. Below this value, QMOM is used.
        scalar sigmaMin_;

        //- Bool to track values of sigma_ that lead to unrealizable moments
        bool foundUnrealizableSigma_;

        //- Bool to track if sigma = 0 is root
        bool nullSigma_;

    // Protected member functions

        //- Reset inverter
        void reset();

		//- Invert moments if the distribution is singular
        void invertSingular(univariateMomentSet& moments);

        //- Compute secondary weigths and abscissae
        void secondaryQuadrature
        (
            const scalarList& pWeights,
            const scalarList& pAbscissae
        );


public:

    //- Runtime type information
    TypeName("extendedMomentInversion");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            extendedMomentInversion,
            dictionary,
            (
                const dictionary& dict,
                const label nMoments,
                const label nSecondaryNodes
            ),
            (dict, nMoments, nSecondaryNodes)
        );

    // Constructors

        // Construct from dictionary and label
        extendedMomentInversion
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
    virtual ~extendedMomentInversion();


    // Selectors

        static autoPtr<extendedMomentInversion> New
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    // Member Functions

        //- Invert moments to find weight, abscissae and sigma
        virtual void invert(const univariateMomentSet& moments) = 0;

        //- Return number of moments
        inline label nMoments()
        {
            return nMoments_;
        }

        //- Return number of primary quadrature nodes
        inline label nPrimaryNodes()
        {
            return nPrimaryNodes_;
        }

        //- Return number of secondary quadrature nodes
        inline label nSecondaryNodes()
        {
            return nSecondaryNodes_;
        }

        //- Return the value of the sigma parameter
        inline scalar sigma() const
        {
            return sigma_;
        }

        //- Return primary quadrature weights
        inline const scalarList& primaryWeights() const
        {
            return primaryWeights_;
        }

        //- Return primary quadrature abscissae
        inline const scalarList& primaryAbscissae() const
        {
            return primaryAbscissae_;
        }

        //- Return secondary quadrature weights
        inline const scalarRectangularMatrix& secondaryWeights() const
        {
            return secondaryWeights_;
        }

        //- Return secondary quadrature abscissae
        inline const scalarRectangularMatrix& secondaryAbscissae() const
        {
            return secondaryAbscissae_;
        }
		
        //- Return the sum of the kernel density functions at a field of points
        tmp<scalarField> f(const scalarField& x) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
