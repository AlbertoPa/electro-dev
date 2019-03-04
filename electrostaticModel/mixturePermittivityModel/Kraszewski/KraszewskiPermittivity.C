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

\*---------------------------------------------------------------------------*/

#include "KraszewskiPermittivity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(KraszewskiPermittivity, 0);

    addToRunTimeSelectionTable
    (
        mixturePermittivityModel,
        KraszewskiPermittivity,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KraszewskiPermittivity::KraszewskiPermittivity(const dictionary& dict)
:
    mixturePermittivityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::KraszewskiPermittivity::~KraszewskiPermittivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::KraszewskiPermittivity::mixturePermittivity
(
    const volScalarField& alphaParticles,
    const dimensionedScalar& particleRelPermittivity,
    const dimensionedScalar& fluidRelPermittivity
) const
{

    return sqr(alphaParticles*sqrt(particleRelPermittivity)
        + (scalar(1) - alphaParticles)*sqrt(fluidRelPermittivity));
}

// ************************************************************************* //
