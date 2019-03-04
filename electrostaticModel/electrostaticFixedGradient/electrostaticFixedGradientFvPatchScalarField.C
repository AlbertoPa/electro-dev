/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "electrostaticFixedGradientFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrostaticFixedGradientFvPatchScalarField::electrostaticFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
    //electrostaticGradient_()
{}


/*Foam::electrostaticFixedGradientFvPatchScalarField::electrostaticFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const Field<scalar>& fld
)
:
    fixedGradientFvPatchScalarField(p, iF, fld)
    //electrostaticGradient_()
{}*/


Foam::electrostaticFixedGradientFvPatchScalarField::electrostaticFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
    //electrostaticGradient_(Function1<Type>::New("electrostaticGradient", dict))
{
    this->evaluate();
}


Foam::electrostaticFixedGradientFvPatchScalarField::electrostaticFixedGradientFvPatchScalarField
(
    const electrostaticFixedGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
    //electrostaticGradient_(ptf.electrostaticGradient_, false)
{}


Foam::electrostaticFixedGradientFvPatchScalarField::electrostaticFixedGradientFvPatchScalarField
(
    const electrostaticFixedGradientFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf)
    //electrostaticGradient_(ptf.electrostaticGradient_, false)
{}


Foam::electrostaticFixedGradientFvPatchScalarField::electrostaticFixedGradientFvPatchScalarField
(
    const electrostaticFixedGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
    //electrostaticGradient_(ptf.electrostaticGradient_, false)
{
     //Evaluate the profile if defined
    //if (ptf.electrostaticGradient_.valid())
    //{
    //    this->evaluate();
    //}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electrostaticFixedGradientFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    
    const objectRegistry& mesh(this->db());
    const dictionary& electrostaticProperties(mesh.lookupObject<dictionary>("electrostaticProperties"));
             
    //const dimensionedScalar& initialRhoq = electrostaticProperties.lookup("initialRhoq");
    const dimensionedScalar& saturatedRhoq = electrostaticProperties.lookup("saturatedRhoq");
    const dimensionedScalar& poissonsRatio = electrostaticProperties.lookup("poissonsRatio");
    const dimensionedScalar& youngsModulus = electrostaticProperties.lookup("youngsModulus");
    const dimensionedScalar& chargingEfficiency = electrostaticProperties.lookup("chargingEfficiency");
    const dimensionedScalar& potentialDifference = electrostaticProperties.lookup("potentialDifference");
    const dimensionedScalar& vacuumPermittivity = electrostaticProperties.lookup("vacuumPermittivity");
    const dimensionedScalar& particleRelPermittivity = electrostaticProperties.lookup("particleRelPermittivity");
    const dimensionedScalar& particleDiameter = electrostaticProperties.lookup("particleDiameter");
    const dimensionedScalar& particleDensity = electrostaticProperties.lookup("particleDensity");
    const dimensionedScalar& alphaMax = electrostaticProperties.lookup("alphaMax");
          
    //const volScalarField& alphas = db().lookupObject<volScalarField>("alpha.particles");
    //const volScalarField& thetas = db().lookupObject<volScalarField>("Theta.particles");
    //const volScalarField& rhoq = db().lookupObject<volScalarField>("rhoq"); 
        
    const fvPatch& patch(this->patch());
       
    const fvPatchScalarField& alphasPatchField(patch.lookupPatchField<volScalarField, scalar>("alpha.particles"));
    
    const fvPatchScalarField& thetasPatchField(patch.lookupPatchField<volScalarField, scalar>("Theta.particles"));
    
    //const fvPatchScalarField& rhoqPatchField(patch.lookupPatchField<volScalarField, scalar>("rhoq"));
    const scalarField rhoqPatchField(patchInternalField());
    
    //calculate Vp, g0, ks, z0, kq, C1, qnew = qinf + (qold -qinf)*exp(-C1/qinf*delt)
    
    scalarField N = max(alphasPatchField/(constant::mathematical::pi*pow(particleDiameter.value(),3)/6.0),1e-10);
        
    scalarField g0 = 1.0/(1.0-(pow((alphasPatchField/alphaMax.value()),1.0/3.0)));

    scalar ks = 1.364*(pow(particleDiameter.value(),2.0))*(pow((particleDensity.value()*(1.0-poissonsRatio.value())/youngsModulus.value()),2.0/5.0));
    
    scalar z0 = particleDiameter.value()/2000.0;
    
    scalarField kq = potentialDifference.value()*chargingEfficiency.value()*vacuumPermittivity.value()*particleRelPermittivity.value()/z0*(1.0-rhoqPatchField/N/saturatedRhoq.value());
    
    this->gradient() = 3.0/5.0*g0*ks*kq*(pow(2.0,19.0/10.0))*(tgamma(2.0/5.0))/(pow(constant::mathematical::pi,3.0/2.0))*alphasPatchField*(pow(thetasPatchField,9.0/10.0))/(pow(particleDiameter.value(),3.0))*alphasPatchField*particleDensity.value();

    fixedGradientFvPatchScalarField::updateCoeffs();    
}


void Foam::electrostaticFixedGradientFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    //electrostaticGradient_->writeData(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        electrostaticFixedGradientFvPatchScalarField
    );
}
// ************************************************************************* //
