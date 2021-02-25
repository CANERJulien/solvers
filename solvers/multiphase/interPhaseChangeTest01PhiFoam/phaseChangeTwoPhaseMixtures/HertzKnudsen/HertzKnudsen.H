/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::phaseChangeTwoPhaseMixtures::HertzKnudsen

Description
    HertzKnudsen cavitation model.

    Reference:
    \verbatim
        C. L. HertzKnudsen, J. Feng, and P. E. O. Buelow,
        "Computational modeling of the dynamics of sheet cavitation",
        in Proceedings Third International Symposium on Cavitation
        Grenoble, France 1998.
   \endverbatim

SourceFiles
    HertzKnudsen.C

\*---------------------------------------------------------------------------*/

#ifndef HertzKnudsen_H
#define HertzKnudsen_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*---------------------------------------------------------------------------*\
                                 Class HertzKnudsen
\*---------------------------------------------------------------------------*/

class HertzKnudsen
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

        
        dimensionedScalar fc_;
        dimensionedScalar fv_;
        dimensionedScalar M_;
        dimensionedScalar lv_;
        dimensionedScalar T0_;
        dimensionedScalar P0_;
        
        dimensionedScalar Rgas_;
         
       dimensionedScalar unitLength_;

        dimensionedScalar mcCoeff_;
        dimensionedScalar mvCoeff_;


public:

    //- Runtime type information
    TypeName("HertzKnudsen");


    // Constructors

        //- Construct from components
        HertzKnudsen
        (
            const volVectorField& U,
            const surfaceScalarField& phi
);

    //- Destructor
    virtual ~HertzKnudsen()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const;

        //- Correct the HertzKnudsen phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
