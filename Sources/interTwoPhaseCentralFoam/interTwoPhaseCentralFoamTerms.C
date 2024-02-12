/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
       hybridCentralSolvers | Copyright (C) 2016-2021 ISP RAS (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.
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
#include "interTwoPhaseCentralFoam.H"
#include "gaussGrad.H"

void Foam::interTwoPhaseCentralFoam::pressureGradient()
{
    surfaceScalarField p_own
    (
        fvc::interpolate(p_rgh_, own_, "reconstruct(p)")
    );
    surfaceScalarField p_nei
    (
        fvc::interpolate(p_rgh_, nei_, "reconstruct(p)")
    );
    surfaceScalarField pf
    (
        linearInterpolate(p_rgh_)
    );

    surfaceVectorField phase1_coeffs
    (
        (
            kappa_*(alpha1_own_ *p_own + alpha1_nei_*p_nei)
            +
            onemkappa_*pf
        ) * p_rgh_.mesh().Sf()
    );

    surfaceVectorField phase2_coeffs
    (
        (
            kappa_*(alpha2_own_ *p_own + alpha2_nei_*p_nei)
            +
            onemkappa_*pf
        ) * p_rgh_.mesh().Sf()
    );

    gradp_ =
        fvc::div(phase1_coeffs)*volumeFraction1_
        +
        fvc::div(phase2_coeffs)*volumeFraction2_;
}

//* * * * * * * * * * * * * * * * * Viscosity * * * * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::divDevRhoReff()
{
    mu_ = volumeFraction1_*mu1_ + volumeFraction2_*mu2_;
    divDevRhoReff_ =
    (
        - fvm::laplacian(mu_, U_)
        - fvc::div((mu_)*dev2(Foam::T(fvc::grad(U_))))
    );
}

//
//END-OF-FILE
//
