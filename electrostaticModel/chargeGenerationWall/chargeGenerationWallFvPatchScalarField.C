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

#include "chargeGenerationWallFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    ),

    vacuumPermittivity_
    (
        electrostaticProperties_.lookup("vacuumPermittivity")
    ),
    particleRelPermittivity_
    (
        electrostaticProperties_.lookup("particleRelPermittivity")
    ),
    fluidRelPermittivity_
    (
        electrostaticProperties_.lookup("fluidRelPermittivity")
    ),
    saturatedRhoq_
    (
        electrostaticProperties_.lookup("saturatedRhoq")
    ),
    poissonsRatio_
    (
        electrostaticProperties_.lookup("poissonsRatio")
    ),
    youngsModulus_
    (
        electrostaticProperties_.lookup("youngsModulus")
    ),
    chargingEfficiencyWall_
    (
        electrostaticProperties_.lookup("chargingEfficiencyWall")
    ),
    potentialDifference_
    (
        electrostaticProperties_.lookup("potentialDifference")
    ),
    particleWorkFunction_
    (
        electrostaticProperties_.lookup("particleWorkFunction")
    ),
    particleDiameter_
    (
        electrostaticProperties_.lookup("particleDiameter")
    ),
    particleDensity_
    (
        electrostaticProperties_.lookup("particleDensity")
    ),
    alphaMax_
    (
        electrostaticProperties_.lookup("alphaMax")
    ),
    z0_
    (
        electrostaticProperties_.lookup("z0")
    ),
    ks_
    (
        1.364*(pow(particleDiameter_,2.0))*(pow((particleDensity_*(1.0-poissonsRatio_)/youngsModulus_),2.0/5.0))
    ),
    kqs_
    (
        ((5.0/7.0)*tgamma(12.0/5.0)/pow(2.0,1.0/10.0)/
         pow(constant::mathematical::pi,0.5))
        *ks_*chargingEfficiencyWall_*vacuumPermittivity_*
        potentialDifference_/z0_/saturatedRhoq_
    ),
    kqe_
    (
        ((5.0/7.0)*tgamma(12.0/5.0)/pow(2.0,1.0/10.0)/
         pow(constant::mathematical::pi,0.5))
        *ks_*chargingEfficiencyWall_*vacuumPermittivity_
    )

{
    this->evaluate();
}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const chargeGenerationWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const chargeGenerationWallFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const chargeGenerationWallFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chargeGenerationWallFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatch& patch(this->patch());

    const scalarField rhoqCell(patchInternalField());

    // Particle volume fraction at the wall boundary
    const fvPatchScalarField& alphasPatchField
    (
        patch.lookupPatchField<volScalarField, scalar>("alpha.particles")
    );

    // Particle theta at the wall boundary
    const fvPatchScalarField& thetasPatchField
    (
        patch.lookupPatchField<volScalarField, scalar>("Theta.particles")
    );

    // Electric field at the wall boundary
    const fvPatchVectorField& EqsPatchField
    (
        patch.lookupPatchField<volVectorField, vector>("Eq")
    );

    const fvPatchScalarField& diffusivityPatchField
    (
        patch.lookupPatchField<volScalarField, scalar>("RhoqDiffusivity")
    );

    const fvPatchScalarField& fieldDiffusivityPatchField
    (
        patch.lookupPatchField<volScalarField, scalar>("FieldDiffusivity")
    );

    // Radial distribution function at the wall boundary
    scalarField g0
    (
        1.0/(1.0-(pow((alphasPatchField/alphaMax_.value()),1.0/3.0)))
    );

    scalarField Kw1
    (
        max(alphasPatchField, scalar(1.0e-6))*particleDensity_.value()*
        pow(max(thetasPatchField, scalar(1.0e-10)),(9.0/10.0))*
        g0*kqs_.value()

    );

    scalarField Kw2
    (
        max(alphasPatchField, scalar(1.0e-6))*particleDensity_.value()*
        pow(max(thetasPatchField, scalar(1.0e-10)),(9.0/10.0))*
        g0*kqe_.value()

    );

    scalarField cRhoq
    (
        diffusivityPatchField*patch.deltaCoeffs()
    );

    scalarField cField
    (
        fieldDiffusivityPatchField
    );

    // Effect of electric field on wall charge BC
    scalarField EqWallNormal
    (
        (EqsPatchField & patch.Sf())/patch.magSf()
    );

    // Compute particle charge at wall boundary
    scalarField qWall
    (
        ((Kw1*saturatedRhoq_.value())+(Kw2*EqWallNormal)+(cField*EqWallNormal)-(cRhoq*rhoqCell))/(Kw1 - cRhoq)
    );

    operator==(qWall);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::chargeGenerationWallFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        chargeGenerationWallFvPatchScalarField
    );
}
// ************************************************************************* //
