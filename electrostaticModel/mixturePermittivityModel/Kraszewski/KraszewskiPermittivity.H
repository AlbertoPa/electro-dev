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
    Foam::KraszewskiPermittivity

Description
    Mixture permittivity model of Looyenga.

    Reference
    \verbatim
        Kraszewski, A., Kulinski, S., & Matuszewski, M. (1976).
        Dielectric properties and a model of biphase water suspension at 9.4 GHz.
        Journal of Applied Physics, 47(4), 1275–1277.
    \endverbatim

SourceFiles
    KraszewskiPermittivity.C

\*---------------------------------------------------------------------------*/

#ifndef KraszewskiPermittivity_H
#define KraszewskiPermittivity_H

#include "mixturePermittivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class KraszewskiPermittivity Declaration
\*---------------------------------------------------------------------------*/

class KraszewskiPermittivity
:
    public mixturePermittivityModel
{

public:

    //- Runtime type information
    TypeName("Kraszewski");


    // Constructors

        //- Construct from components
        KraszewskiPermittivity(const dictionary& dict);


    //- Destructor
    virtual ~KraszewskiPermittivity();


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
