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
    Foam::LooyengaPermittivity

Description
    Mixture permittivity model of Looyenga.

    Reference
    \verbatim
        Looyenga, H. (1965).
        Dielectric constants of heterogeneous mixtures.
        Physica, 31(3), 401–406.
    \endverbatim

SourceFiles
    LooyengaPermittivity.C

\*---------------------------------------------------------------------------*/

#ifndef LooyengaPermittivity_H
#define LooyengaPermittivity_H

#include "mixturePermittivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LooyengaPermittivity Declaration
\*---------------------------------------------------------------------------*/

class LooyengaPermittivity
:
    public mixturePermittivityModel
{

public:

    //- Runtime type information
    TypeName("Looyenga");


    // Constructors

        //- Construct from components
        LooyengaPermittivity(const dictionary& dict);


    //- Destructor
    virtual ~LooyengaPermittivity();


    // Member Functions

        tmp<volScalarField> mixturePermittivity
        (
            const volScalarField& alphaParticles,
            const dimensionedScalar& particleRelPermittivity,
            const dimensionedScalar& fluidRelPermittivity
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
