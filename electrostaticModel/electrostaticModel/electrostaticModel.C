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

#include "electrostaticModel.H"
#include "fvCFD.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrostaticModel::electrostaticModel
(
    const Foam::twoPhaseSystem& phaseSystem
)
:
    phaseSystem_(phaseSystem),
    phase_(phaseSystem.phase1()),
    alpha_(phase_),
    rho_(phase_.rho()),
    alphaRhoPhi_(phase_.alphaRhoPhi()),
    phiParticle_(phase_.phi()),
    theta_
    (
        phaseSystem.mesh().lookupObject<volScalarField>
        (IOobject::groupName("Theta", phase_.name()))
    ),

    electrostaticProperties_
    (
        IOobject
        (
            "electrostaticProperties",
            alpha_.time().constant(),
            alpha_.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    turbulenceParticleProperties_
    (
        IOobject
        (
            "turbulenceProperties.particles",
            alpha_.time().constant(),
            alpha_.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    electrostatics_
    (
        electrostaticProperties_.lookup("electrostatics")
    ),
    mixturePermittivityModel_
    (
        mixturePermittivityModel::New
        (
            electrostaticProperties_
        )
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
    chargingEfficiencyParticle_
    (
        electrostaticProperties_.lookup("chargingEfficiencyParticle")
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
        0.783*(pow(particleDiameter_,2.0))*
        (pow((particleDensity_*(1.0-poissonsRatio_)/youngsModulus_),2.0/5.0))
    ),
    kq_
    (
        chargingEfficiencyParticle_*particleDiameter_/6.0       // Laurentie
    ),
    e_
    (
        "e", dimless,
        turbulenceParticleProperties_.subDict("RAS").
        subDict("kineticTheoryCoeffs").lookup("e")
    ),
    g_
    (
        alpha_.mesh().lookupObject<uniformDimensionedVectorField>("g")
    ),

    Vq_
    (
        IOobject
        (
            "Vq",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh()
    ),
    Eq_
    (
        IOobject
        (
            "Eq",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        ),
        -fvc::grad(Vq_)

    ),

    q_flux_
    (
        IOobject
        (
            "Qflux",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedVector("dimqflux", dimensionSet(1,-5,0,0,0,1,0), Zero)
    ),

    q_flux_no_rhoq
    (
        IOobject
        (
            "QNoRhoq",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedVector("dimqflux", dimensionSet(1,-5,0,0,0,1,0), Zero)
    ),

    diffusivity
    (
        IOobject
        (
            "RhoqDiffusivity",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("dimRhoqDiffusivity", dimensionSet(1,-1,-1,0,0,0,0), Zero)
    ),

    rhoq_
    (
        IOobject
        (
            "rhoq",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh()
    ),

    varRhoq_
    (
        IOobject
        (
            "varRhoq",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("dimVarRhoq", dimensionSet(0,-6,2,0,0,2,0), Zero)
    ),

    cq_
    (
        IOobject
        (
            "CQ",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedVector("dimCQ", dimensionSet(0,-2,0,0,0,1,0), Zero)
    ),

    rhoPearson_
    (
        IOobject
        (
            "PearsonCorrelation",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimless, 1.0)
    ),

    Fq_
    (
        IOobject
        (
            "Fq",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhoq_*alpha_*Eq_
    ),

    PE_
    (
        IOobject
        (
            "PE",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhoq_*Vq_
    ),

    mixtureRelPermittivity_
    (
        IOobject
        (
            "mixtureRelPermittivity",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimless, 1.0)
    )

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electrostaticModel::~electrostaticModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electrostaticModel::correct()
{
    // Radial distribution function
    volScalarField g0
    (
        1.0/(1.0-(pow((alpha_/alphaMax_),1.0/3.0)))
    );

    // Limiting minimum value of alpha to a non-zero positive number
    volScalarField alpha
    (
        max(alpha_, scalar(1.0e-6))
    );

    // Setting minimum non-zero value for theta
    dimensionedScalar thetaMin
    (
        "thetaMin",dimensionSet(0,2,-2,0,0,0,0),scalar(1.0e-10)
    );

    volScalarField theta
    (
        max(theta_,thetaMin)
    );

    // Compute mixture relative permittivity
    mixtureRelPermittivity_ =
        mixturePermittivityModel_->mixturePermittivity
        (
            alpha_,
            particleRelPermittivity_,
            fluidRelPermittivity_
        );

    // Compute mixture absolute permittivity
    volScalarField mixtureEffPermittivity
    (
        vacuumPermittivity_*mixtureRelPermittivity_
    );

    // Solve for electric potential
    fvScalarMatrix VqEqn
    (
        fvm::laplacian(mixtureEffPermittivity, Vq_) + alpha_*rhoq_
    );

    //VqEqn.relax();
    VqEqn.solve();

    // Compute electrical potential energy density
    PE_ = rhoq_*Vq_;

    // Compute electric field
    Eq_ = -fvc::grad(Vq_);

    // * * * * *  Compute coefficients for charge diffusion  * * * * * * //

    dimensionedScalar cq_k1
    (
        ks_*kq_             //MR (10/17/2018) Changed the particle charging model to Laurentie
    );

    dimensionedScalar cq_k2
    (
        ks_*chargingEfficiencyParticle_*vacuumPermittivity_
    );

    volScalarField thetaPow0_9
    (
        pow(theta,(9.0/10.0))
    );

    volScalarField thetaPow0_4
    (
        pow(theta,(2.0/5.0))
    );

    volScalarField thetaPow1_3
    (
        pow(theta,(13.0/10.0))
    );

    volScalarField thetaPow0_8
    (
        pow(theta,(4.0/5.0))
    );

    volScalarField cq_A
    (
        "cq_A",
        ((1.0+e_)*pow(constant::mathematical::pi*theta,0.5))
        +(pow(2.0,(14.0/5.0))*tgamma(12.0/5.0)
        *pow(constant::mathematical::pi,0.5)*cq_k1
        *thetaPow0_9*(5.0/7.0 - (5.0*(1.0+e_)/12.0)))
    );

    volScalarField cq1
    (
        "cq1",
        (1.0/cq_A)*pow(2.0,(19.0/5.0))*tgamma(29.0/10.0)*
        pow(constant::mathematical::pi,0.5)*cq_k2*theta*
        thetaPow0_4*((5.0/57.0)-(23*(1.0+e_)/264.0))
    );

    volScalarField cq2
    (
        "cq2",
        (1.0/cq_A)*5.0/108.0*pow(2.0,14.0/5.0)*tgamma(12.0/5.0)
        *particleDiameter_*pow(constant::mathematical::pi,0.5)*cq_k2
        *thetaPow0_9*alpha/mixtureEffPermittivity
    );

    volScalarField cq3
    (
        "cq3",
        (1.0/cq_A)*
        (((constant::mathematical::pi*particleDiameter_)/(g0*6.0*alpha))
         +(0.5*(1.0+e_)*particleDiameter_))*theta
    );

    // MR (10/17/2018) - Updated cq4,5,6 to correct missing alpha
    volScalarField cq4
    (
        "cq4",
        g0*(60.0/19.0)*pow(2.0, 9.0/5.0)*tgamma(19.0/10.0)*
        (1.0/pow(constant::mathematical::pi,0.5))
        *cq_k1*thetaPow0_4*alpha_
    );

    volScalarField cq5
    (
        "cq5",
        g0*(5.0/9.0)*pow(2.0, 9.0/5.0)*tgamma(12.0/5.0)*
        (1.0/pow(constant::mathematical::pi,0.5))*cq_k2*
        thetaPow0_9*alpha_
    );

    volScalarField cq6
    (
        "cq6",
        g0*(5.0/21.0)*pow(2.0, 9.0/5.0)*tgamma(12.0/5.0)*
        (1.0/pow(constant::mathematical::pi,0.5))*cq_k1*
        thetaPow0_9*alpha_*
        particleDiameter_
    );

    volScalarField cq7
    (
        (1.0/cq_A)*
        ((constant::mathematical::pi*particleDiameter_)/
        (g0*6.0*alpha))*(rhoq_*rhoq_/particleDensity_)
    );

    volScalarField cq_g
    (
        "cq_g",
        (1.0/cq_A)*
        ((constant::mathematical::pi*particleDiameter_)/
        (g0*6.0*alpha))
    );

    volScalarField qq_A
    (
        "qq_A",
        g0*alpha*alpha*rho_*(1.0/particleDiameter_)*cq_k1*
        (cq_k1*(5.0/3.0)*(1.0/pow(constant::mathematical::pi,0.5))
        *pow(2.0,(28.0/5.0))*tgamma(14.0/5.0)*thetaPow1_3
        - (15.0/7.0)*(1.0/pow(constant::mathematical::pi,0.5))
        *pow(2.0,(24.0/5.0))*tgamma(12.0/5.0)*thetaPow0_9)
    );

    volScalarField qq_B
    (
        "qq_B",
        g0*alpha*alpha*rho_*cq_k2*cq_k2*(1.0/particleDiameter_)
        *(39.0/35.0)*(1.0/pow(constant::mathematical::pi,0.5))
        *pow(2.0,(18.0/5.0))*tgamma(14.0/5.0)*thetaPow1_3
    );

    volScalarField qq_C
    (
        "qq_C",
        g0*alpha*alpha*rho_*(1.0/particleDiameter_)*cq_k2*(cq_k1
        *(15.0/23.0)*(1.0/pow(constant::mathematical::pi,0.5))
        *pow(2.0,(28.0/5.0))*tgamma(23.0/10.0)*thetaPow0_8
        - (15.0/19.0)*(1.0/pow(constant::mathematical::pi,0.5))
        *pow(2.0,(24.0/5.0))*tgamma(19.0/5.0)*thetaPow0_4)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Compute self diffusion flux
    cq_ = - (cq1+cq7)*Eq_
          - cq3*fvc::grad(rhoq_)
          + cq_g*rhoq_*g_;

    // Compute variance of charge
    varRhoq_ = (-1.0/qq_A)*(qq_B*(magSqr(Eq_)) + qq_C*(Eq_ & cq_));

    // Setting minimum non-zero value for variance
    dimensionedScalar varRhoqMin
    (
        "varRhoqMin",dimensionSet(0,-6,2,0,0,2,0),scalar(1.0e-4)
    );

    varRhoq_ = max(varRhoq_, 0.0*varRhoqMin);

    volScalarField varRhoq
    (
        max(varRhoqMin, varRhoq_)
    );

    // Compute flux of charge without contribution of gradient of rhoq
    q_flux_no_rhoq =
            (fvc::reconstruct(alphaRhoPhi_))*rhoq_
            + (alpha_*rho_*(cq_g+(cq_g*cq4)))*rhoq_*g_
            - (alpha_*rho_*(-cq1-(cq1*cq4)+cq5+cq7+(cq7*cq4)))*fvc::grad(Vq_)
            + (alpha_*rho_*cq6*(rhoq_/varRhoq)*fvc::grad(varRhoq_));

    // Compute diffusivity of rhoq alone
    diffusivity = alpha_*rho_*((cq3+(cq3*cq4))+cq6);

    // Solve for charge density
    fvScalarMatrix rhoqEqn
    (
        // AP Time derivative in non-conservative form to avoid singularity
        //    when alpha -> 0.

        fvm::ddt(alpha*rho_, rhoq_)
      + rho_*rhoq_*fvc::ddt(alpha_)
      + fvm::div(alphaRhoPhi_, rhoq_, "div(" + alphaRhoPhi_.name() + ",rhoq)")
      + (g_ & fvc::grad(alpha_*rho_*(cq_g+(cq_g*cq4))*rhoq_))
      - fvc::laplacian(alpha_*rho_*(-cq1-(cq1*cq4)+cq5+cq7+(cq7*cq4)), Vq_)
      - fvm::laplacian(alpha*rho_*((cq3+(cq3*cq4))+cq6),rhoq_)
      + fvc::laplacian(alpha*rho_*cq6*rhoq_/varRhoq, varRhoq_)
      - fvm::Sp(phase_.continuityError(), rhoq_)
    );

    rhoqEqn.relax();
    rhoqEqn.solve();

    Info << "Charge: min = " << min(rhoq_).value()
         << " - max = " << max(rhoq_).value()
         << " - average = " << sum(alpha_*rhoq_).value()/sum(alpha_).value()
         << endl;

    // Compute electric force acting on the particle phase (per unit volume)
    Fq_ = rhoq_*alpha_*Eq_;

    Info << "Electric potential: min(Vq) = " << min(Vq_).value()
         << " - max(Vq) = " << max(Vq_).value() << endl;

    Info << "Electrostatic field magnitude: " << max(mag(Eq_)).value() << endl;

    Info << "Electrostatic force magnitude: " << max(mag(Fq_)).value() << endl;

}

// ************************************************************************* //
