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
    Foam::electrostaticModel

Description
    Electrostatic model to describe electrostatic phenomena in fluidized
    systems.

SourceFiles
    electrostaticModel.C

\*---------------------------------------------------------------------------*/

#ifndef electrostaticModel_H
#define electrostaticModel_H

#include "mixturePermittivityModel.H"
#include "autoPtr.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class electrostaticModel Declaration
\*---------------------------------------------------------------------------*/

class electrostaticModel
{
    // Private data

        //- Particle phase subject to electrostatic forces
        const twoPhaseSystem& phaseSystem_;

        const phaseModel& phase_;

        //- Volume fraction field of the phase subject to electrostatic effects
        const volScalarField& alpha_;

        //- Density field of the phase subject to electrostatic effects
        const volScalarField& rho_;

        //- Flux of the phase subject to electrostatic effects
        const surfaceScalarField& alphaRhoPhi_;

        //- Surface flux of velocity
        const surfaceScalarField& phiParticle_;

        //- m2/s2
        const volScalarField& theta_;

        //- dictionary holding the modeling info
        IOdictionary electrostaticProperties_;

        //- dictionary holding the modeling info
        IOdictionary turbulenceParticleProperties_;

        //- flag to check if electrostatic model is activated
        Switch electrostatics_;

        autoPtr<mixturePermittivityModel> mixturePermittivityModel_;

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
        const dimensionedScalar chargingEfficiencyParticle_;

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

        //- Particle-particle charge diffusion constant
        const dimensionedScalar kq_;

        //Coefficient of restitution
        const dimensionedScalar e_;

        const uniformDimensionedVectorField g_;

        //- Electric potential
        volScalarField Vq_;

        //- Electric field
        volVectorField Eq_;

        // Diffusivity of charge
        volScalarField diffusivity;

        // Diffusivity of charge due to electric field
        volScalarField fieldDiffusivity;

        //- Charge density
        volScalarField rhoq_;

        //- Electrostatic force acting on the particle phase
        volVectorField Fq_;

        //- Electric potential energy density
        volScalarField PE_;

        //- Mixture relative electrical permittivity
        volScalarField mixtureRelPermittivity_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        electrostaticModel(const electrostaticModel&);

        //- Disallow default bitwise assignment
        void operator=(const electrostaticModel&);


public:

    // Constructors

        //- Construct from components
        electrostaticModel
        (
            const twoPhaseSystem& phaseSystem
        );


    //- Destructor
	virtual ~electrostaticModel();

    // Member Functions

        //- Solve for charge density and electric potential
        void correct();

        //- Return status of electrostatic model
        bool on() const
        {
            return electrostatics_;
        }

        //- Return electric potential
        const volScalarField& Vq() const
        {
            return Vq_;
        }

        //- Return electric field
        const volVectorField& Eq() const
        {
            return Eq_;
        }

        //- Return electric force
        const volVectorField& Fq() const
        {
            return Fq_;
        }

        //- Return electric potential energy density
        const volScalarField& PE() const
        {
            return PE_;
        }

        //- Return relative permittivity of the multiphase mixture
        const volScalarField& mixtureRelPermittivity() const
        {
            return mixtureRelPermittivity_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
