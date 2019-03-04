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
    Foam::electrostaticFixedGradientFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a electrostatic fixed gradient condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        electrostaticGradient | electrostatic gradient     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            electrostaticFixedGradient;
        
    }
    \endverbatim

Note


See also
    Foam::Function1Types
    Foam::fixedGradientFvPatchField

SourceFiles
    electrostaticFixedGradientFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef electrostaticFixedGradientFvPatchScalarField_H
#define electrostaticFixedGradientFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class electrostaticFixedGradientFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

//template<class Type>
class electrostaticFixedGradientFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Gradient
        //autoPtr<Function1<Type>> electrostaticGradient_;


public:

    //- Runtime type information
    TypeName("electrostaticFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        electrostaticFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        /*electrostaticFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const Field<scalar>& fld
        );*/

        //- Construct from patch, internal field and dictionary
        electrostaticFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given electrostaticFixedGradientFvPatchScalarField
        //  onto a new patch
        electrostaticFixedGradientFvPatchScalarField
        (
            const electrostaticFixedGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electrostaticFixedGradientFvPatchScalarField
        (
            const electrostaticFixedGradientFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new electrostaticFixedGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        electrostaticFixedGradientFvPatchScalarField
        (
            const electrostaticFixedGradientFvPatchScalarField&,
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
                new electrostaticFixedGradientFvPatchScalarField(*this, iF)
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
