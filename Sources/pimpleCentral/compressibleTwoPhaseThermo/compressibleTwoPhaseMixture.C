/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
       hybridCentralSolvers | Copyright (C) 2016-2018 ISP RAS (www.unicfd.ru)
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

#include "compressibleTwoPhaseMixture.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixture::compressibleTwoPhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    liqPhaseName_(dict.lookupOrDefault<word>("LiqPhaseName", "Liq")),
    gasPhaseName_(dict.lookupOrDefault<word>("GasPhaseName", "Gas")),

    YLiq_
    (
        IOobject
        (
            ("Y." + liqPhaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    YGas_
    (
        IOobject
        (
            ("Y." + gasPhaseName_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - YLiq_
    ),

    YbarLiq_
    (
        IOobject
        (
            ("Ybar." + liqPhaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    YbarGas_
    (
        IOobject
        (
            ("Ybar." + gasPhaseName_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - YbarLiq_
    ),

    rhoEff_
    (
        IOobject
        (
            "thermo:rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    )
{
    YGas_ = 1.0 - YLiq_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibleTwoPhaseMixture::updateVolFrac(const volScalarField& rhoLiq, const volScalarField& rhoGas)
{
    YbarLiq_ = (YLiq_ * rhoGas / rhoLiq) / (1.0 - YLiq_ + YLiq_ * rhoGas / rhoLiq);
    YbarGas_ = 1.0 - YbarLiq_;
}


// ************************************************************************* //
