/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "quasiIncompressibleTwoPhaseMixture.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quasiIncompressibleTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quasiIncompressibleTwoPhaseMixture::quasiIncompressibleTwoPhaseMixture
(
            const IOdictionary &dict,
            const volScalarField &volumeFraction1,
            const volScalarField &volumeFraction2,
            const volScalarField &p,
            const volScalarField &T
)
:
    transportModel(),
    volumeFraction1_(volumeFraction1),
    volumeFraction2_(volumeFraction2),
    p_(p),
    T_(T),
    R1_
    (
        dimensioned<scalar>("R1", dict)
    ),
    R2_
    (
        dimensioned<scalar>("R2", dict)
    ),
    molM1_
    (
        dimensioned<scalar>("molM1", dict)
    ),
    molM2_
    (
        dimensioned<scalar>("molM2", dict)
    ),
    rho01_
    (
        dimensioned<scalar>("rho01", dict)
    ),
    rho02_
    (
        dimensioned<scalar>("rho02", dict)
    ),
    rho1Min_
    (
        dimensioned<scalar>("rho1Min", dict)
    ),
    rho2Min_
    (
        dimensioned<scalar>("rho2Min", dict)
    ),
    Cp1_
    (
        dimensioned<scalar>("Cp1", dict)
    ),
    Cp2_
    (
        dimensioned<scalar>("Cp2", dict)
    ),
    gamma1_
    (
        dimensioned<scalar>
        (
            "gamma1",
            Cp1_/(Cp1_ - (Foam::constant::physicoChemical::R/molM1_))
        )
    ),
    gamma2_
    (
        dimensioned<scalar>
        (
            "gamma2",
            Cp2_/(Cp2_ - (Foam::constant::physicoChemical::R/molM2_))
        )
    ),
    mu1_
    (
        dimensioned<scalar>("mu1", dict)
    ),
    mu2_
    (
        dimensioned<scalar>("mu2", dict)
    ),
    Pr1_
    (
        dimensioned<scalar>("Pr1", dict)
    ),
    Pr2_
    (
        dimensioned<scalar>("Pr2", dict)
    ),
    alpha1_
    (
        "alpha1",
        mu1_/(Pr1_)
    ),
    alpha2_
    (
        "alpha2",
        mu2_/(Pr2_)
    ),
    Prt_
    (
        dimensioned<scalar>("Prt", dict)
    ),
    psi1_
    (
        "psi1",
        1/(R1_*T_)
    ),
    psi2_
    (
        "psi2",
        1/(R2_*T_)
    ),
    rho1_
    (
        IOobject
        (
            "rho1",
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
         psi1_*p_
    ),
    rho2_
    (
        IOobject
        (
            "rho2",
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        psi2_*p_
    ),
    rho0_
    (
        "rho0",
        volumeFraction1_*rho01_
        +
        volumeFraction2_*rho02_
    ),
    rho_
    (
        "rho",
        volumeFraction1_*rho1_
        +
        volumeFraction2_*rho2_
    ),
    mu_
    (
        "mu",
        volumeFraction1_*mu1_
        +
        volumeFraction2_*mu2_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::quasiIncompressibleTwoPhaseMixture::~quasiIncompressibleTwoPhaseMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::quasiIncompressibleTwoPhaseMixture::nu() const
{
    tmp<volScalarField> nu_mixture = mu_ / rho_;
    nu_mixture.ref().rename("nu");
    return nu_mixture;
}


Foam::tmp<Foam::scalarField>
Foam::quasiIncompressibleTwoPhaseMixture::nu(const label patchi) const
{
    tmp<scalarField> nu_mixture_patch
    (
        mu_.boundaryField()[patchi]/rho_.boundaryField()[patchi]
    );
    return nu_mixture_patch;
}


void Foam::quasiIncompressibleTwoPhaseMixture::correct()
{
    rho_ = volumeFraction1_*rho1_ + volumeFraction2_*rho2_;
    rho0_ = volumeFraction1_*rho01_ + volumeFraction2_*rho02_;
    mu_ = volumeFraction1_*mu1_ + volumeFraction2_*mu2_;
}

void Foam::quasiIncompressibleTwoPhaseMixture::compressibility()
{
    psi1_ = 1.0/(R1_*T_);
    psi2_ = 1.0/(R2_*T_);
}

void Foam::quasiIncompressibleTwoPhaseMixture::density()
{
    rho1_ = psi1_*p_ + rho01_;
    rho2_ = psi2_*p_ + rho02_;
    rho1_ = max(rho1_,rho1Min_);
    rho2_ = max(rho2_,rho2Min_);
    rho1_.correctBoundaryConditions();
    rho2_.correctBoundaryConditions();
}

void Foam::quasiIncompressibleTwoPhaseMixture::saveOldTime()
{
    rho1_.oldTime();
    rho2_.oldTime();
    rho_.oldTime();
    psi1_.oldTime();
    psi2_.oldTime();
}

void Foam::quasiIncompressibleTwoPhaseMixture::limitRho1()
{
    rho1_ = max(rho1_, rho1Min_);
}

void Foam::quasiIncompressibleTwoPhaseMixture::limitRho2()
{
    rho2_ = max(rho2_, rho2Min_);
}

// bool Foam::quasiIncompressibleTwoPhaseMixture::read()
// {
//     if (regIOobject::read())
//     {
//         return viscosityModelPtr_->read(*this);
//     }

//     return false;
// }


// ************************************************************************* //
