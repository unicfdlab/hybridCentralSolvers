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

void Foam::interTwoPhaseCentralFoam::pressureGradient()
{
    surfaceScalarField p_own = fvc::interpolate(p_, own_, "reconstruct(p)");
    surfaceScalarField p_nei = fvc::interpolate(p_, nei_, "reconstruct(p)");
    // gradp_ = fvc::div((alpha_own_ *p_own + alpha_nei_*p_nei)*U_.mesh().Sf());

//    gradp_ = fvc::div((alpha1_own_ *p_own + alpha1_nei_*p_nei)*U_.mesh().Sf());
//    gradp_ = fvc::div((alpha2_own_ *p_own + alpha2_nei_*p_nei)*U_.mesh().Sf());
//    gradp_ = fvc::grad(p_);

    gradp_ = volumeFraction1_*
        fvc::div((alpha1_own_ *p_own + alpha1_nei_*p_nei)*U_.mesh().Sf());
    gradp_ += volumeFraction2_*
        fvc::div((alpha2_own_ *p_own + alpha2_nei_*p_nei)*U_.mesh().Sf());
}

//* * * * * * * * * * * * * * * * * Viscosity * * * * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::divDevRhoReff()
{
    divDevRhoReff1_ =
    (
        - fvm::laplacian(mu1_, U_)
        - fvc::div((mu1_)*dev2(Foam::T(fvc::grad(U_))))
    );

    divDevRhoReff2_ =
    (
        - fvm::laplacian(mu2_, U_)
        - fvc::div((mu2_)*dev2(Foam::T(fvc::grad(U_))))
    );
}


void Foam::interTwoPhaseCentralFoam::viscosityTEqn()
{
    Tviscosity1 = - fvm::laplacian(alpha1_*Cp1_, T_);

    Tviscosity2 =- fvm::laplacian(alpha2_*Cp2_, T_);
}


void Foam::interTwoPhaseCentralFoam::devRhoReff()
{
    devRhoReff1_ = (-(alpha1_)*dev(twoSymm(fvc::grad(U_))));
    devRhoReff2_ = (-(alpha2_)*dev(twoSymm(fvc::grad(U_))));
    // Rename TSourse
    TSource1_ =
    fvc::div((linearInterpolate((-devRhoReff1_) & U_) & U_.mesh().Sf())());

    TSource2_ =
    fvc::div((linearInterpolate((-devRhoReff2_) & U_) & U_.mesh().Sf())());
}

//* * * * * * * * * * * * * * * * * * Others * * * * * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::divU()
{

    surfaceScalarField rbyAf = fvc::interpolate(rbyA_);
/*
    Foam::CorrectPhi
    (
        U_,
        phi_,
        p_,
        rbyAf,
        divU_,
        pimple_
    );
*/
}

//
//END-OF-FILE
//

