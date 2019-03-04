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

Class
    Foam::chargeGenerationWallFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an updated charge on the wall due to generation.

Usage
    \table
        Property     | Description             | Required    | Default value
        chargeGenerationWall | electrostatic charge generation at wall     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            chargeGenerationWall;

    }
    \endverbatim

Note


See also
    Foam::fixedGradientFvPatchField

SourceFiles
    chargeGenerationWallFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chargeGenerationWallFvPatchScalarField_H
#define chargeGenerationWallFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class chargeGenerationWallFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

//template<class Type>
class chargeGenerationWallFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Gradient
        //autoPtr<Function1<Type>> electrostaticGradient_;

    const objectRegistry& mesh_;

    const dictionary& electrostaticProperties_;

    //- Vacuum electrical permittivityparticleRelPermittivitymix
    const dimensionedScalar vacuumPermittivity_;

    //- Particle relative electrical permittivity
    const dimensionedScalar particleRelPermittivity_;

    //- Fluid relative electrical permittivity
    const dimensionedScalar fluidRelPermittivity_;

    //- Saturated charge density of the particle
    const dimensionedScalar saturatedRhoq_;

    //- Poisson ration for the particle material
    const dimensionedScalar poissonsRatio_;

    //- Young's modulus for the particle material
    const dimensionedScalar youngsModulus_;

    //- Charging efficiency for the particle-particle charging model
    const dimensionedScalar chargingEfficiencyWall_;

    //- Work function of the particle material
    const dimensionedScalar potentialDifference_;

     //- Work function of the particle material
    const dimensionedScalar particleWorkFunction_;

    //- Diameter of the particle
    const dimensionedScalar particleDiameter_;

    //- Particle material density
    const dimensionedScalar particleDensity_;

    //- Maximum packing volume fraction
    const dimensionedScalar alphaMax_;

    //- Electrostatic cutoff distance
    const dimensionedScalar z0_;

    //- Hertzian deformation area model constant
    const dimensionedScalar ks_;

    //- Particle-wall charge diffusion constants
    const dimensionedScalar kqs_;
    const dimensionedScalar kqe_;

    // Electric field at the wall boundary
    //const fvPatchVectorField& EqsPatchField;



public:

    //- Runtime type information
    TypeName("chargeGenerationWall");


    // Constructors

        //- Construct from patch and internal field
        chargeGenerationWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        /*chargeGenerationWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const Field<scalar>& fld
        );*/

        //- Construct from patch, internal field and dictionary
        chargeGenerationWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chargeGenerationWallFvPatchScalarField
        //  onto a new patch
        chargeGenerationWallFvPatchScalarField
        (
            const chargeGenerationWallFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        chargeGenerationWallFvPatchScalarField
        (
            const chargeGenerationWallFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chargeGenerationWallFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        chargeGenerationWallFvPatchScalarField
        (
            const chargeGenerationWallFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new chargeGenerationWallFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
