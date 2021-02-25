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

\*---------------------------------------------------------------------------*/

#include "HertzKnudsen.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(HertzKnudsen, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, HertzKnudsen, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::HertzKnudsen::HertzKnudsen
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    
    fc_("fc", dimless, phaseChangeTwoPhaseMixtureCoeffs_),
    fv_("fv", dimless, phaseChangeTwoPhaseMixtureCoeffs_),
    M_("M", dimensionSet (1, 0, 0, 0, -1, 0, 0), phaseChangeTwoPhaseMixtureCoeffs_),

    lv_("lv", dimensionSet (0, 2, -2, 0, 0, 0, 0), phaseChangeTwoPhaseMixtureCoeffs_),
    T0_("T0", dimensionSet (0, 0, 0, 1, 0, 0, 0), phaseChangeTwoPhaseMixtureCoeffs_),
    P0_("P0", dimensionSet (1, -1, -2, 0, 0, 0, 0), phaseChangeTwoPhaseMixtureCoeffs_),  
    
    Rgas_("Rgas", dimensionSet (1, 2, -2, -1, -1, 0, 0), phaseChangeTwoPhaseMixtureCoeffs_),

    unitLength_("unitLength", dimensionSet (0, 1, 0, 0, 0, 0, 0), phaseChangeTwoPhaseMixtureCoeffs_),
    
    
    mcCoeff_(scalar(2)*fc_/(scalar(2)-fc_)*sqrt(M_/(scalar(2*3.141459265358932384)*Rgas_*T0_))/unitLength_),
    mvCoeff_(scalar(2)*fv_/(scalar(2)-fv_)*sqrt(M_/(scalar(2*3.141459265358932384)*Rgas_*T0_))/unitLength_)
    
    
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::HertzKnudsen::mDotAlphal() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pSatu = alpha1_.db().lookupObject<volScalarField>("pSatu");    
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");


    dimensionedScalar zeroPres
(
"zeroPres",
dimensionSet(1,-1,-2,0,0,0,0),
scalar(0)
);

    forAll(pSatu, cellI)
{
pSatu[cellI]=P0_.value()*exp((M_.value()*lv_.value()/Rgas_.value())*(scalar(1)/(T0_.value())-scalar(1)/T[cellI]));

} 
    
   // pSatu = P0_*exp((M_*lv_/Rgas_)*(scalar(1)/T0_-scalar(1)/T));

    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*max((p-pSatu),zeroPres),
        mvCoeff_*min((p- pSatu),zeroPres)
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::HertzKnudsen::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    
    volScalarField pSatu = alpha1_.db().lookupObject<volScalarField>("pSatu");
    
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    
 

dimensionedScalar Rgas
(
"Rgas",
dimensionSet(1, 2, -2, -1, -1, 0, 0),
scalar(8.31447)
);

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    
    forAll(pSatu, cellI)
{
pSatu[cellI]=P0_.value()*exp((M_.value()*lv_.value()/Rgas_.value())*(scalar(1)/(T0_.value())-scalar(1)/T[cellI]));
}
   // pSatu = P0_*exp((M_*lv_/Rgas_)*(scalar(1)/T0_-scalar(1)/T));

    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*(1.0 - limitedAlpha1)*pos0(p-pSatu),
        (-mvCoeff_)*limitedAlpha1*neg(p-pSatu)      
       
    );
}


void Foam::phaseChangeTwoPhaseMixtures::HertzKnudsen::correct()
{}


bool Foam::phaseChangeTwoPhaseMixtures::HertzKnudsen::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");

        
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("fc") >> fc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("fv") >> fv_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("M") >> M_;
        //phaseChangeTwoPhaseMixtureCoeffs_.lookup("lv") >> lv_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("T0") >> T0_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("P0") >> P0_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Rgas") >> Rgas_;
       
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("unitLength") >> unitLength_;

        
        
        mcCoeff_ = scalar(2)*fc_/(scalar(2)-fc_)*sqrt(M_/(scalar(2*3.141459265358932384)*Rgas_*T0_))/unitLength_;
        mvCoeff_ = scalar(2)*fv_/(scalar(2)-fv_)*sqrt(M_/(scalar(2*3.141459265358932384)*Rgas_*T0_))/unitLength_;
        
        
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
