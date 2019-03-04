/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2016 Alberto Passalacqua
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
    Foam::mixturePermittivityModel

Description
    Mixture electrical permittivity models for gas-particle flows

SourceFiles
    mixturePermittivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef mixturePermittivityModel_H
#define mixturePermittivityModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mixturePermittivityModel Declaration
\*---------------------------------------------------------------------------*/

class mixturePermittivityModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        mixturePermittivityModel(const mixturePermittivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const mixturePermittivityModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("mixturePermittivityModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        mixturePermittivityModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        mixturePermittivityModel(const dictionary& dict);


    // Selectors

        static autoPtr<mixturePermittivityModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~mixturePermittivityModel();


    // Member Functions

        //- Mixture permittivity
        virtual tmp<volScalarField> mixturePermittivity
        (
            const volScalarField& alphaParticles,
            const dimensionedScalar& particleRelPermittivity,
            const dimensionedScalar& fluidRelPermittivity
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
